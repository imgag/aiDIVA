import argparse
import networkx as nx
import numpy as np
import os
import pandas as pd
import multiprocessing as mp
import pickle
import re
from itertools import combinations
from operator import itemgetter
from scipy.stats import poisson

if not __name__=="__main__":
    from . import get_HPO_similarity_score as gs


coding_variants = ["splice_acceptor_variant",
                   "splice_donor_variant",
                   "stop_gained",
                   "frameshift_variant",
                   "stop_lost",
                   "start_lost",
                   "inframe_insertion",
                   "inframe_deletion",
                   "missense_variant",
                   "protein_altering_variant",
                   "splice_region_variant",
                   "incomplete_terminal_codon_variant",
                   "start_retained_variant",
                   "stop_retained_variant",
                   "synonymous_variant",
                   "coding_sequence_variant",
                   "5_prime_UTR_variant",
                   "3_prime_UTR_variant"]

supported_coding_variants = ["stop_gained",
                             "stop_lost",
                             "start_lost",
                             "inframe_insertion",
                             "inframe_deletion",
                             "missense_variant",
                             "protein_altering_variant",
                             "incomplete_terminal_codon_variant",
                             "start_retained_variant",
                             "stop_retained_variant",
                             "coding_sequence_variant"]

cadd_identifier = "CADD_PHRED"
duplication_identifier = "segmentDuplication"
repeat_identifier = "simpleRepeat"
family = None
family_type = "SINGLE"
genes2exclude = None
gene_2_HPO = None
hgnc_2_gene = None
gene_2_interacting = None
HPO_graph = None
HPO_query = None
HPO_query_distances = 0


def prioritize_variants(variant_data, hpo_resources_folder, num_cores, family_file=None, fam_type="SINGLE", hpo_list=None, gene_exclusion_list=None):
    #load HPO resources
    gene_2_HPO_f = hpo_resources_folder + "gene2hpo.pkl"
    hgnc_2_gene_f = hpo_resources_folder + "hgnc2gene.pkl"
    gene_2_interacting_f = hpo_resources_folder + "gene2interacting.pkl"
    HPO_graph_file = hpo_resources_folder + "hpo_graph.pkl"
    hpo_list_file = hpo_list
    gene_exclusion_file = gene_exclusion_list

    global hgnc_2_gene
    hgnc_2_gene = pickle.load(open(hgnc_2_gene_f, "rb"))

    global gene_2_interacting
    gene_2_interacting = pickle.load(open(gene_2_interacting_f, "rb"))

    global gene_2_HPO
    gene_2_HPO = pickle.load(open(gene_2_HPO_f, "rb"))
    hpo_nodes, hpo_edges = pickle.load(open(HPO_graph_file, "rb"))

    global HPO_graph
    HPO_graph = nx.Graph()
    HPO_graph.add_nodes_from(hpo_nodes)
    HPO_graph.add_edges_from(hpo_edges)

    global genes2exclude
    genes2exclude = set()
    if gene_exclusion_file:
        if os.path.isfile(gene_exclusion_file):
            with open(gene_exclusion_file, "r") as exclusion_file:
                for line in exclusion_file:
                    if line.startswith("#"):
                        continue
                    if line == "\n":
                        continue
                    gene = line.rstrip().split("\t")[0]
                    genes2exclude.add(gene)
        else:
            print("The specified gene exclusion list %s is not a valid file" % (gene_exclusion_file))
            print("No genes are excluded during filtering!")

    global HPO_query
    global HPO_query_distances
    HPO_query = set()
    if hpo_list_file:
        if os.path.isfile(hpo_list_file):
            with open(hpo_list_file, "r") as hpo_file:
                for line in hpo_file:
                    HPO_term = line.rstrip()
                    HPO_query.add(HPO_term)
            HPO_query = list(HPO_query) # removes duplicate entries in the list
            HPO_query.sort() # makes sure that the gene symbols are ordered (could lead to problems otherwise)
            HPO_query_distances = gs.precompute_query_distances(HPO_graph, HPO_query)
        else:
            print("The specified HPO list %s is not a valid file" % (hpo_list_file))
            print("Skip HPO score finalization!")

    # read family relationship (from PED file)
    global family
    family = dict()
    if family_file:
        if os.path.isfile(family_file):
            with open(family_file, "r") as fam_file:
                for line in fam_file:
                    line = line.rstrip()
                    splitline = line.split("\t")

                    if splitline[5] == "2":
                       family[splitline[1]] = 1
                    elif splitline[5] == "1":
                       family[splitline[1]] = 0
                    else:
                       print("ERROR: There is a problem with the given PED file describing the family.")
        else:
            print("The specified family file %s is not a valid file" % (family_file))
            print("Skip inheritance assessment!")

    global family_type
    family_type = fam_type

    variant_data = parallelize_dataframe_processing(variant_data, parallelized_variant_processing, num_cores)
    variant_data = variant_data.sort_values(["FINAL_AIDIVA_SCORE"], ascending=[False])
    variant_data = variant_data.reset_index(drop=True)

    return variant_data


def parallelize_dataframe_processing(variant_data, function, num_cores):
    num_partitions = num_cores * 2

    if len(variant_data) <= num_partitions:
        dataframe_splitted = np.array_split(variant_data, 1)
    else:
        dataframe_splitted = np.array_split(variant_data, num_partitions)

    pool = mp.Pool(num_cores)
    variant_data = pd.concat(pool.map(function, dataframe_splitted))
    pool.close()
    pool.join()

    return variant_data


def parallelized_variant_processing(variant_data):
    variant_data = check_inheritance(variant_data, family_type, family)
    variant_data[["HPO_RELATEDNESS", "HPO_RELATEDNESS_INTERACTING", "FINAL_AIDIVA_SCORE"]] = variant_data.apply(lambda variant: pd.Series(compute_hpo_relatedness_and_final_score(variant)), axis=1)
    variant_data[["FILTER_PASSED", "FILTER_COMMENT"]] = variant_data.apply(lambda variant: pd.Series(check_filters(variant)), axis=1)

    return variant_data


def compute_hpo_relatedness_and_final_score(variant):
    if HPO_query:
        if np.isnan(variant["AIDIVA_SCORE"]):
            final_score = np.nan
            hpo_relatedness = np.nan
            hpo_relatedness_interacting = np.nan
        else:
            variant_gene = str(variant["SYMBOL"])
            hgnc_id = str(variant["HGNC_ID"])
            gene_distances = []
            gene_distances_interacting = []
            processed_HPO_genes = dict()

            if variant_gene in gene_2_interacting.keys():
                interacting_genes = [gene_interaction["interacting_gene"] for gene_interaction in gene_2_interacting[variant_gene]]
            else:
                interacting_genes = []

            # we use the hgnc ID to prevent problems if a given gene symbol isn't used anymore           
            if variant_gene not in genes2exclude:
                if variant_gene in gene_2_HPO.keys():
                    gene_HPO_list = gene_2_HPO.get(variant_gene, [])
                else:
                    if (hgnc_id != "nan") and (hgnc_id in hgnc_2_gene.keys()):
                        gene_symbol = hgnc_2_gene[hgnc_id]
                        gene_HPO_list = gene_2_HPO.get(gene_symbol, [])
                    else:
                        print("WARNING: Given gene is not covered!")
                        gene_HPO_list = []

                g_dist = gs.list_distance(HPO_graph, HPO_query, gene_HPO_list, HPO_query_distances)
                gene_distances.append(g_dist)

                ## TODO: look at interacting genes only if HPO relation missing or zero???
                for interacting_gene in interacting_genes:
                    if interacting_gene in genes2exclude:
                        continue
                    if interacting_gene in gene_2_HPO.keys():
                        gene_HPO_list = gene_2_HPO.get(interacting_gene, [])
                    else:
                        print("WARNING: Given gene is not covered!")
                        gene_HPO_list = []

                    g_dist = gs.list_distance(HPO_graph, HPO_query, gene_HPO_list, HPO_query_distances)
                    gene_distances_interacting.append(g_dist)

                if gene_distances or gene_distances_interacting:
                    # take only the maximum HPO relatedness to prevent downvoting of genes if no HPO relation in interacting genes is observed
                    hpo_relatedness = max(gene_distances, default=0.0)
                    hpo_relatedness_interacting = max(gene_distances_interacting, default=0.0)
                    ## TODO: try different weighting of AIDIVA_SCORE and HPO_RELATEDNESS and HPO_RELATEDNESS_INTERACTING (eg 0.6 and 0.3 and 0.1)
                    # predicted pathogenicity has a higher weight than the HPO relatedness
                    final_score = (float(variant["AIDIVA_SCORE"]) * 0.67 + float(max(hpo_relatedness, hpo_relatedness_interacting)) * 0.33) #+ float(hpo_relatedness_interacting) * 0.1)# / 3
                    #final_score = (float(variant["AIDIVA_SCORE"]) + float(hpo_relatedness) + float(hpo_relatedness_interacting)) / 3
            else:
                final_score = np.nan
                hpo_relatedness = np.nan
                hpo_relatedness_interacting = np.nan
    else:
        final_score = variant["AIDIVA_SCORE"]
        hpo_relatedness = np.nan
        hpo_relatedness_interacting = np.nan

    return [hpo_relatedness, hpo_relatedness_interacting, final_score]


def check_inheritance(variant_data, family_type="SINGLE", family=None):
    variant_columns = variant_data.columns
    if family_type == "SINGLE":
        variant_data["COMPOUND"] = 0
        variant_data["RECESSIVE"] = variant_data.apply(lambda variant: check_recessive_single(variant, variant_columns), axis=1)

        variant_data_grouped = [group for key, group in variant_data.groupby("SYMBOL")]

        for group in variant_data_grouped:
            check_compound_single(group, variant_columns)
        variant_data = pd.concat(variant_data_grouped)

    elif family_type == "TRIO":
        if not family is None:
            variant_data["COMPOUND"] = 0
            variant_data["DOMINANT_DENOVO"] = variant_data.apply(lambda variant: check_denovo(variant, family), axis=1)

            ## TODO: do we need a check for affected family members?
            variant_data["DOMINANT"] = variant_data.apply(lambda variant: check_dominant(variant, family), axis=1)
            variant_data["XLINKED"] = variant_data.apply(lambda variant: check_xlinked(variant, family), axis=1)
            variant_data["RECESSIVE"] = variant_data.apply(lambda variant: check_recessive(variant, family, family_type), axis=1)

            variant_data_grouped = [group for key, group in variant_data.groupby("SYMBOL")]

            affected_child = ""
            parent_1 = ""
            parent_2 = ""

            ## TODO: handle the case if affected status is given as unknown
            ## TODO: change to be less strict and to not rely on unaffected family members
            for name in family.keys():
                if family[name] == 1:
                    affected_child = name
                elif family[name] == 0:
                    if not parent_1:
                        parent_1 = name
                        continue
                    elif not parent_2:
                        parent_2 = name
                        continue
                    else:
                        print("Something went wrong!")
            if  affected_child and parent_1 and parent_2:
                for group in variant_data_grouped:
                    check_compound(group, affected_child, parent_1, parent_2)
            variant_data = pd.concat(variant_data_grouped)
        else:
            print("ERROR: If family type (TRIO) is used a proper PED file defining the family relations is required")
        
    elif (family_type == "FAMILY") and (family is not None):
        pass
        if family is not None:
            variant_data["DOMINANT"] = variant_data.apply(lambda variant: check_dominant(variant, family), axis=1)
            variant_data["XLINKED"] = variant_data.apply(lambda variant: check_xlinked(variant, family), axis=1)
            variant_data["RECESSIVE"] = variant_data.apply(lambda variant: check_recessive(variant, family, family_type), axis=1)
        else:
            print("ERROR: If family type (FAMILY) is used a proper PED file defining the family relations is required")
        
    else:
        print("ERROR: Unsupported family type (%s) was given!" % family_type)

    variant_columns = variant_data.columns
    variant_data["INHERITANCE"] = variant_data.apply(lambda variant: add_inheritance_mode(variant, variant_columns), axis=1)

    return variant_data


def add_inheritance_mode(variant, variant_columns):
    inheritance_list = []

    if "DOMINANT" in variant_columns:
        if variant["DOMINANT"] == 1:
            inheritance_list.append("DOMINANT")

    if "DOMINANT_DENOVO" in variant_columns:
        if variant["DOMINANT_DENOVO"] == 1:
            inheritance_list.append("DOMINANT_DENOVO")
    
    if "RECESSIVE" in variant_columns:
        if variant["RECESSIVE"] == 1:
            inheritance_list.append("RECESSIVE")
    
    if "COMPOUND" in variant_columns:
        if variant["COMPOUND"] == 1:
            inheritance_list.append("COMPOUND")

    if "XLINKED" in variant_columns:
        if variant["XLINKED"] == 1:
            inheritance_list.append("XLINKED")
    
    inheritance_mode = "&".join(inheritance_list)

    return inheritance_mode


def check_filters(variant):
    variant_genes = re.sub("\(.*?\)", "", str(variant["SYMBOL"]))
    genenames = set(variant_genes.split(";"))

    consequences = str(variant["Consequence"])
    found_consequences = consequences.split("&")
    seg_dup = float(variant[duplication_identifier])
    repeat = str(variant[repeat_identifier])
    cadd = float(variant[cadd_identifier])
    filter_comment = ""

    try:
        maf = variant["MAX_AF"]
    except Exception as e:
        print("Allele frequency could not be identified, use 0.0 instead")
        maf = 0.0

    # exclude gene, if it is on the exclusion list
    if len(genes2exclude & genenames) > 0:
        for gene in genenames:
            if gene in genes2exclude:
                filter_passed = 0 # gene in exclusion list
                filter_comment = "gene exclusion"
                return filter_passed, filter_comment

    if (repeat != "NA") and (repeat != "") and (repeat != "nan"):
        filter_passed = 0 # tandem repeat
        filter_comment = "tandem repeat"
        return filter_passed, filter_comment

    ## TODO: change frequency based on inheritance mode (hom/het)
    if float(maf) <= 0.01:
        if any(term for term in coding_variants if term in found_consequences):
            if "synonymous_variant" not in found_consequences:
                if ("frameshift_variant" not in found_consequences) and (any(term for term in supported_coding_variants if term in found_consequences)):
                    if not np.isnan(variant["AIDIVA_SCORE"]):
                        if len(HPO_query) >= 1:
                            if float(variant["HPO_RELATEDNESS"]) > 0.0:
                                filter_passed = 1
                                filter_comment = "passed all"
                            elif float(variant["HPO_RELATEDNESS_INTERACTING"]) > 0.0:
                                filter_passed = 1
                                filter_comment = "HPO related to interacting genes"
                            else:
                                filter_passed = 0 # no relation to reported HPO terms
                                filter_comment = "no HPO relation"
                        else:
                            filter_passed = 1 # skip hpo filter if no terms are present
                            filter_comment = "no HPO terms given"
                    else:
                        filter_passed = 0 # no prediction present (eg. variant type not covered by the used ML models)
                        filter_comment = "missing AIDIVA_SCORE"
                else:
                    filter_passed = 0 # splicing, intron and frameshift variants are not covered by the used ML models (the predictions cannot be trusted), skip variants that have no covered consequence
                    filter_comment = "variant type not covered"
            else:
                filter_passed = 0 # synonymous variant
                filter_comment = "synonymous"
        else:
            filter_passed = 0 # not coding
            filter_comment = "not coding"
    else:
        filter_passed = 0 # allele frequency to high
        filter_comment = "high MAF"

    return filter_passed, filter_comment


def check_compound(gene_variants, affected_child, parent_1, parent_2):
    num_variant_candidates = gene_variants.shape[0]

    if num_variant_candidates >= 2:
        candidate_indices = [x for x in combinations(gene_variants.index.tolist(), 2)]
        for index_pair in candidate_indices:
            affected_child_zygosity_a = gene_variants.loc[index_pair[0], "GT." + affected_child].replace("|", "/")
            affected_child_zygosity_b = gene_variants.loc[index_pair[1], "GT." + affected_child].replace("|", "/")

            parent_1_zygosity_a = gene_variants.loc[index_pair[0], "GT." + parent_1].replace("|", "/")
            parent_1_zygosity_b = gene_variants.loc[index_pair[1], "GT." + parent_1].replace("|", "/")

            parent_2_zygosity_a = gene_variants.loc[index_pair[0], "GT." + parent_2].replace("|", "/")
            parent_2_zygosity_b = gene_variants.loc[index_pair[1], "GT." + parent_2].replace("|", "/")

            if (affected_child_zygosity_a == "0/1") and (affected_child_zygosity_b == "0/1"):
                if ((parent_1_zygosity_a == "0/0") and (parent_2_zygosity_a == "0/1")) and ((parent_1_zygosity_b == "0/1") and (parent_2_zygosity_b == "0/0")):
                    gene_variants.loc[index_pair[0], "COMPOUND"] = 1
                    gene_variants.loc[index_pair[1], "COMPOUND"] = 1

                elif ((parent_1_zygosity_a == "0/1") and (parent_2_zygosity_a == "0/0")) and ((parent_1_zygosity_b == "0/0") and (parent_2_zygosity_b == "0/1")):
                    gene_variants.loc[index_pair[0], "COMPOUND"] = 1
                    gene_variants.loc[index_pair[1], "COMPOUND"] = 1

            if ("." in affected_child_zygosity_a) or ("." in affected_child_zygosity_b):
                print("WARNING: Skip variant pair, uncalled genotype in affected sample!")



def check_compound_single(gene_variants, variant_columns):
    genotype_column = [column for column in variant_columns if column.startswith("GT.")][0]
    num_variant_candidates = gene_variants.shape[0]

    if num_variant_candidates >= 2:
        candidate_indices = [x for x in combinations(gene_variants.index.tolist(), 2)]
        for index_pair in candidate_indices:
            genotype_variant_a = gene_variants.loc[index_pair[0], genotype_column].replace("|", "/")
            genotype_variant_b = gene_variants.loc[index_pair[1], genotype_column].replace("|", "/")

            if (((genotype_variant_a == "0/1") and (genotype_variant_b == "0/1"))):
                gene_variants.loc[index_pair[0], "COMPOUND"] = 1
                gene_variants.loc[index_pair[1], "COMPOUND"] = 1


def check_denovo(variant, family):
    judgement = 0
    check_samples = dict()

    # create data structure for completeness check
    for name in family.keys():
        check_samples[name] = 0

        zygosity = variant["GT." + name].replace("|", "/")

        # check if sample is found in pedigree
        # sample info complete?
        ## TODO: do we need this check???
        if name in check_samples:
            check_samples[name] = 1

        # heterozygous in affected individual - good
        if zygosity == "0/1" and family[name] == 1:
            judgement = 1
            continue

        # hom ref, not affected - good
        elif zygosity == "0/0" and family[name] == 0 :
            judgement = 1
            continue

        # heterozygous in non-affected - bad
        elif zygosity == "0/1" and family[name] == 0:
            judgement = 0
            break

        # hom ref in affected - bad
        elif zygosity == "0/0" and family[name] == 1:
            judgement = 0
            break

        # homozygous can"t be denovo
        elif zygosity == "1/1":
            judgement = 0
            break

        else:
            # reject missing genotype here or beforehand
            pass
        
    for vals in check_samples.values():
       if vals == 0:
            judgement = 0
            break

    return judgement


def check_dominant(variant, family):
    judgement = 0
    check_samples = dict()

    # create data structure for completeness check
    for name in family.keys():
        check_samples[name] = 0

        zygosity = variant["GT." + name].replace("|", "/")

        if name in check_samples:
            check_samples[name] = 1

        # affected family members should have the mutation (hom ref not allowed)
        if zygosity == "0/0" and family[name] == 1:
            judgement = 0
            break

        # affected family members might be het
        elif zygosity == "0/1" and family[name] == 1:
            judgement = 1
            continue

        # affected family members might be hom alt
        # that"s the major difference to de novo...
        elif zygosity == "1/1" and family[name] == 1:
            judgement = 1
            continue

        # non-affected family members must not have the mutation - hom ref is OK
        elif zygosity == "0/0" and family[name] == 0:
            judgement = 1
            continue

        # non-affected family members must not have the mutation - het is bad
        elif zygosity == "0/1" and family[name] == 0:
            judgement = 0
            break

        # non-affected family members must not have the mutation - hom alt is worst
        elif zygosity == "1/1" and family[name] == 0:
            judgement = 0
            break
        
        else:
            # reject missing genotype here or beforehand
            pass

    for vals in check_samples.values():
       if vals == 0:
            judgement = 0
            break

    return judgement


def check_dominant_single(variant, variant_columns):
    genotype_column = [column for column in variant_columns if column.startswith("GT.")][0]
    judgement = 0

    if (variant[genotype_column] == "0/1"):
                judgement = 1

    return judgement


def check_recessive(variant, family, family_type):
    judgement = 0
    check_samples = dict()

    # create data structure for completeness check
    for name in family.keys():
        check_samples[name] = 0

        zygosity = variant["GT." + name].replace("|", "/")

        if name in check_samples:
            check_samples[name] = 1

        # affected individuals have to be homozygous
        if zygosity == "1/1" and family[name] == 1:
            judgement = 1
            continue

        # affected individuals should not be hom ref or het
        elif ( zygosity == "0/0" or zygosity == "0/1" ) and family[name] == 1:
            judgement = 0
            break

        # non-affected individuals might be het
        elif zygosity == "0/1" and family[name] == 0:
            judgement = 1
            continue

        # non-affected individuals might be hom ref, if a family is interrogated
        elif zygosity == "0/0" and family[name] == 0 and family_type == "FAMILY":
            judgement = 1
            continue

        # non-affected individuals in a trio are the parents and have to be het
        elif zygosity == "0/0" and family[name] == 0 and family_type == "TRIO":
            judgement = 0
            break

        # non-affected individuals must not be hom alt
        elif zygosity == "1/1" and family[name] == 0:
            judgement = 0
            break

        else:
            # reject missing genotype here or beforehand
            pass
        
    for vals in check_samples.values():
        if vals == 0:
            judgement = 0
            break

    return judgement


def check_recessive_single(variant, variant_columns):
    genotype_column = [column for column in variant_columns if column.startswith("GT.")][0]
    variant_genotype = variant[genotype_column].replace("|", "/")
    is_recessive = 0

    # recessive variants have to be homozygous in affected samples
    if (variant_genotype == "1/1"):
        is_recessive = 1
    
    return is_recessive


def check_xlinked(variant, family):
    judgement = 0
    check_samples = dict()
    inheritance_logic = dict()

    if not ((variant["CHROM"] == "X") or (variant["CHROM"] == "x") or (variant["CHROM"] == "chrX") or (variant["CHROM"] == "chrx") or (variant["CHROM"] == "23")):
        return 0

    # create data structure for completeness check
    for name in family.keys():
        check_samples[name] = 0

        zygosity = variant["GT." + name].replace("|", "/")

        if name in check_samples:
            check_samples[name] = 1

        if family[name] == 0:
            inheritance_logic[name] = zygosity

        # affected individuals have to be homozygous
        if zygosity == "1/1" and family[name] == 1:
            judgement = 1
            continue

        # affected individuals should not be hom ref or het
        elif ( zygosity == "0/0" or zygosity == "0/1" ) and family[name] == 1:
            judgement = 0
            break

        # non-affected individuals might be het
        elif zygosity == "0/1" and family[name] == 0:
            judgement = 1
            continue

        # non-affected individuals might be hom ref
        elif zygosity == "0/0" and family[name] == 0:
            judgement = 1
            continue

        # non-affected individuals must not be hom alt
        elif zygosity == "1/1" and family[name] == 0:
            judgement = 0
            break

        else:
            # reject missing genotype here or beforehand
            pass
        
    # sanity check
    het_checker = 0
    hom_checker = 0

    for values in inheritance_logic.values():
        # mother
        if values == "0/1":
            het_checker = 1
        
        # father
        if values == "0/0":
            hom_checker = 1

    if het_checker == 1 and hom_checker == 1:
        judgement = 1
    else:
        judgement = 0

    # another sanity check
    for vals in check_samples.values():
        if vals == 0:
            judgement = 0
            break

    return judgement


if __name__=="__main__":
    import get_HPO_similarity_score as gs

    parser = argparse.ArgumentParser(description = "Filter variants and finalize the AIDIVA_SCORE based on the given HPO terms (if this information is present)")
    parser.add_argument("--in_file", type=str, dest="in_file", required=True, help="Tab separated input annotated and scored file [required]")
    parser.add_argument("--out_file", type=str, dest="out_filename", required=True, help="Name to save the results [required]")
    parser.add_argument("--family", type=str, dest="family", required=False, help="Tab separated list of samples annotated with affection status. [required]")
    parser.add_argument("--family_type", type=str, choices=["TRIO", "FAMILY", "SINGLE"], dest="family_type", required=False, help="Choose if the data you provide is a trio or a larger family [required]")
    parser.add_argument("--gene_exclusion", type=str, dest="gene_exclusion_list", required=False, help="List of genes that should be excluded in the prioritization")
    parser.add_argument("--hpo_list", type=str, dest="hpo_list", default=None, required=False, help="List of HPO terms that are observed in the patient. These terms are used to adjust the AIDIVA_SCORE\n")
    parser.add_argument("--hpo_resources", type=str, dest="hpo_resources", default="../../data/", required=True, help="Folder where the HPO resources (HPO_graph,...) are found\n")
    args = parser.parse_args()

    input_data = pd.read_csv(args.in_file, sep="\t", low_memory=False)

    prioritized_variants = prioritize_variants(input_data, args.hpo_resources, args.family, args.family_type, args.hpo_list, args.gene_exclusion_list)
    prioritized_variants.to_csv(args.out_filename, sep="\t", index=False)
