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
HPO_graph = None
HPO_query = None
query_dist = 0
num_partitions = 10
num_cores = 5


def prioritize_variants(variant_data, hpo_resources_folder, family_file=None, fam_type="SINGLE", hpo_list=None, gene_exclusion_list=None, n_cores=1):
    #load HPO resources
    gene_2_HPO_f = hpo_resources_folder + "gene2hpo.pkl"
    hgnc_2_gene_f = hpo_resources_folder + "hgnc2gene.pkl"
    HPO_graph_file = hpo_resources_folder + "hpo_graph.pkl"
    hpo_list_file = hpo_list
    gene_exclusion_file = gene_exclusion_list

    global hgnc_2_gene
    hgnc_2_gene = pickle.load(open(hgnc_2_gene_f, "rb"))

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
                for gene in exclusion_file:
                    gene = gene.rstrip()
                    genes2exclude.add(gene)
        else:
            print("The specified gene exclusion list %s is not a valid file" % (gene_exclusion_file))
            print("No genes are excluded during filtering!")

    global HPO_query
    global query_dist
    HPO_query = set()
    if hpo_list_file:
        if os.path.isfile(hpo_list_file):
            with open(hpo_list_file, "r") as hpo_file:
                for line in hpo_file:
                    HPO_term = line.rstrip()
                    HPO_query.add(HPO_term)
            HPO_query = list(HPO_query) # removes duplicate entries in the list
            HPO_query.sort() # makes sure that the gene symbols are ordered (could lead to problems otherwise)
            ## TODO: the following method is obsolete, "list_distance" can be used instead
            query_dist = gs.precompute_query_distances(HPO_graph, HPO_query, 0)
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
                    if line.startswith("sample"):
                        continue
                    line = line.rstrip()
                    splitline = line.split("\t")

                    if splitline[5] == 2:
                       family[splitline[1]] = 1
                    elif splitline[5] == 1:
                       family[splitline[1]] = 0
                    else:
                       print("ERROR: There is a problem with the given PED file describing the family.")
        else:
            print("The specified family file %s is not a valid file" % (family_file))
            print("Skip inheritance assessment!")

    global family_type
    family_type = fam_type

    variant_data = parallelize_dataframe_processing(variant_data, parallelized_variant_processing, n_cores)

    variant_data.sort_values(["FINAL_AIDIVA_SCORE"], ascending=[False], inplace=True)
    variant_data.reset_index(inplace=True, drop=True)

    return variant_data


def parallelize_dataframe_processing(variant_data, function, n_cores):
    if n_cores is None:
        num_cores = 1
    else:
        num_cores = n_cores

    global num_partitions
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
    variant_data[["HPO_RELATEDNESS", "FINAL_AIDIVA_SCORE"]] = variant_data.apply(lambda variant: pd.Series(compute_hpo_relatedness_and_final_score(variant)), axis=1)
    variant_data[["FILTER_PASSED", "FILTER_COMMENT"]] = variant_data.apply(lambda variant: pd.Series(check_filters(variant)), axis=1)

    # Perform inheritance check only if the family information is passed
    if family and family_type != "SINGLE":
        variant_data = variant_data.apply(lambda variant: pd.Series(check_inheritance(variant_data)), axis=1)

    return variant_data


def compute_hpo_relatedness_and_final_score(variant):
    if HPO_query:
        if np.isnan(variant["AIDIVA_SCORE"]):
            final_score = np.nan
            hpo_relatedness = np.nan
        else:
            variant_genes = re.sub("\(.*?\)", "", str(variant["SYMBOL"])) # shouldn't be needed
            hgnc_id = str(variant["HGNC_ID"])
            genenames = set(variant_genes.split(";"))
            gene_distances = []
            processed_HPO_genes = dict()

            # we use the hgnc ID to prevent problems if a given gene symbol isn't used anymore
            for gene_id in genenames:
                if gene_id in genes2exclude:
                    continue
                if gene_id in processed_HPO_genes.keys():
                    gene_distances.append(processed_HPO_genes[gene_id])
                else:
                    if gene_id in gene_2_HPO.keys():
                        gene_HPO_list = gs.extract_HPO_related_to_gene(gene_2_HPO, gene_id)
                    else:
                        if (hgnc_id != "nan") & (hgnc_id in hgnc_2_gene.keys()):
                            gene_symbol = hgnc_2_gene[hgnc_id]
                            gene_HPO_list = gs.extract_HPO_related_to_gene(gene_2_HPO, gene_symbol)
                        else:
                            print("WARNING: Given gene is not covered!")
                            gene_HPO_list = []

                    # do we need to update query_dist here???
                    global query_dist
                    (g_dist, query_distance) = gs.list_distance(HPO_graph, HPO_query, gene_HPO_list, query_dist)
                    query_dist = query_distance

                    gene_distances.append(g_dist)
                    processed_HPO_genes[gene_id] = g_dist

            if gene_distances:
                hpo_relatedness = str(max(gene_distances, default=0.0))
                ## TODO: try different weighting of AIDIVA_SCORE and HPO_RELATEDNESS (eg 0.7 and 0.3)
                # final_score = str((float(variant["AIDIVA_SCORE"]) * 0.7) + (float(hpo_relatedness) * 0.3))
                final_score = str((float(variant["AIDIVA_SCORE"]) + float(hpo_relatedness)) / 2)
            else:
                final_score = variant["AIDIVA_SCORE"]
                hpo_relatedness = np.nan
    else:
        final_score = variant["AIDIVA_SCORE"]
        hpo_relatedness = np.nan

    return [hpo_relatedness, final_score]


def check_inheritance(variant):
    variant_data["COMPOUND"] = 0
    variant_data["DOMINANT_DENOVO"] = variant_data.apply(lambda variant: check_denovo(variant, family), axis=1)

    ## TODO: do we need a check for affected family members?
    variant_data["DOMINANT_INHERITED"] = variant_data.apply(lambda variant: check_dominant(variant, family), axis=1)
    variant_data["XLINKED"] = variant_data.apply(lambda variant: check_xlinked(variant, family), axis=1)
    variant_data["RECESSIVE"] = variant_data.apply(lambda variant: check_recessive(variant, family, family_type), axis=1)

    if family_type == "TRIO":
        variant_data_grouped = [group for key, group in variant_data.groupby("SYMBOL")]

        affected_child = ""
        parent_1 = ""
        parent_2 = ""

        for name in family.keys():
            if family[name] == "1":
                affected_child = name
            elif family[name] == "0":
                if not parent_1:
                    parent_1 = name
                    continue
                if not parent_2:
                    parent_2 = name
                    continue
                else:
                    print("Something went wrong!")
        if  affected_child and parent_1 and parent_2:
            for group in variant_data_grouped:
                check_compound(group, affected_child, parent_1, parent_2)

    variant_data = pd.concat(variant_data_grouped)

    return variant_data


def check_filters(variant):
    variant_genes = re.sub("\(.*?\)", "", str(variant["SYMBOL"]))
    genenames = set(variant_genes.split(";"))

    consequences = str(variant["Consequence"])
    found_consequences = consequences.split("&")
    seg_dup = float(variant[duplication_identifier])
    tandem = str(variant[repeat_identifier])
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

    if (tandem != "NA") & (tandem != "") & (tandem != "nan"):
        filter_passed = 0 # tandem repeat
        filter_comment = "tandem repeat"
        return filter_passed, filter_comment

    ## TODO: filter later compound only less than 0.01
    ## TODO: change frequency based on inheritance mode (hom/het)
    if float(maf) <= 0.02:
        if any(term for term in coding_variants if term in found_consequences):
            if not "synonymous_variant" in found_consequences:
                if ("frameshift_variant" not in found_consequences) & (any(term for term in supported_coding_variants if term in found_consequences)):
                    if not np.isnan(variant["AIDIVA_SCORE"]):
                        if len(HPO_query) >= 1:
                            if float(variant["HPO_RELATEDNESS"]) > 0.0:
                                filter_passed = 1
                                filter_comment = "passed all"
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
            if (((gene_variants.loc[index_pair[0], parent_1] == "0/0") & (gene_variants.loc[index_pair[0], parent_2] == "0/1") & (gene_variants.loc[index_pair[0], affected_child] == "0/1")) & ((gene_variants.loc[index_pair[1], parent_1] == "0/1") & (gene_variants.loc[index_pair[1], parent_2] == "0/0") & (gene_variants.loc[index_pair[1], affected_child] == "0/1"))) | (((gene_variants.loc[index_pair[0], parent_1] == "0/1") & (gene_variants.loc[index_pair[0], parent_2] == "0/0") & (gene_variants.loc[index_pair[0], affected_child] == "0/1")) & ((gene_variants.loc[index_pair[1], parent_1] == "0/0") & (gene_variants.loc[index_pair[1], parent_2] == "0/1") & (gene_variants.loc[index_pair[1], affected_child] == "0/1"))):
                gene_variants.loc[index_pair[0], "COMPOUND"] = 1
                gene_variants.loc[index_pair[1], "COMPOUND"] = 1


def check_denovo(variant, family):
    judgement = 0
    check_samples = dict()

    # create data structure for completeness check
    for name in family.keys():
        check_samples[name] = 0

        if "REF." + name in variant.index.tolist() and "ALT." + name in variant.index.tolist():
            zygosity = variant[name]
            if variant["REF." + name] == ".":
                refcoverage = variant["REF." + name].replace(".", "0") # could be numeric or .
            else:
                refcoverage = variant["REF." + name]
            if variant["ALT." + name] == ".":
                altcoverage = variant["ALT." + name].replace(".", "0") # could be numeric or .
            else:
                altcoverage = variant["ALT." + name]
        else:
            #stick with genotype and the others are empty
            zygosity = variant[name]
            refcoverage = "."
            altcoverage = "."

        # check if sample is found in pedigree
        # sample info complete?
        ## TODO: do we need this check???
        if name in check_samples:
            check_samples[name] = 1

        # heterozygous in affected individual - good
        if zygosity == "0/1" and family[name] == "1":
            judgement = 1
            continue

        # hom ref, not affected - good
        elif zygosity == "0/0" and family[name] == "0" :
            if int(altcoverage) <= max(3, (float(altcoverage) + float(refcoverage)) / 10) :
                judgement = 1
            else :
                judgement =0
            continue

        # heterozygous in non-affected - bad
        elif zygosity == "0/1" and family[name] == "0":
            judgement = 0
            break

        # hom ref in affected - bad
        elif zygosity == "0/0" and family[name] == "1":
            judgement = 0
            break

        # homozygous can"t be denovo
        elif zygosity == "1/1":
            judgement = 0
            break

        # now a few more complex steps, if genotype is missing (only non-affected individuals should have missing values)
        elif zygosity == "./." and family[name] == "0":

            # if vcf file was not supplemented by pileup data
            # reject it variants which could not be called in the parents
            if refcoverage == "." or altcoverage == "." or refcoverage == "" or altcoverage == "":
                judgement = 0
                continue

            # which chance has the current read distribution to miss out on an alt read
            # e.g. 10ref, 0alt
            # if coverage > 8, the chance is less than 0.5% to miss out on one alt read
            # http://stattrek.com/online-calculator/poisson.aspx
            refcoverage = float(refcoverage)
            altcoverage = float(altcoverage)
            coverage = refcoverage + altcoverage

            if coverage == 0:
                judgement = 0
                break

            # hom ref, non called genotype
            # poisson for low coverage and percentage for high coverage
            # poisson 10 reads (poisson average rate of success = 5) and alt reads = 0 - should get still accepted
            elif (poisson.cdf( float(altcoverage), float(coverage)/2 ) <= 0.007) and (altcoverage / coverage <= 0.05):
                judgement = 1
                continue

            # not necessary to check for hom alt
            # coverage too low?
            else:
                judgement = 0
                break

        # do not accept missing values for affected individuals
        elif zygosity == "./." and family[name] == "1":

            # except if vcf file was not supplemented by pileup data
            # reject variants which could not be called in the parents
            if refcoverage == "." or altcoverage == ".":
                judgement = 0
                continue

            # do not be that grateful, if there is coverage data
            # if the SNP caller could not call a genotype in an area, where there was coverage
            # we rather trust the SNP caller, than starting to call SNP based on pileup coverage
            else:
                judgement = 0
                break

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

        if "REF." + name in variant.index.tolist() and "ALT." + name in variant.index.tolist():
            zygosity = variant[name]
            if variant["REF." + name] == ".":
                refcoverage = variant["REF." + name].replace(".", "0") # could be numeric or .
            else:
                refcoverage = variant["REF." + name]
            if variant["ALT." + name] == ".":
                altcoverage = variant["ALT." + name].replace(".", "0") # could be numeric or .
            else:
                altcoverage = variant["ALT." + name]
        else:
            #stick with genotype and the others are empty
            zygosity = variant[name]
            refcoverage = "."
            altcoverage = "."

        if name in check_samples:
            check_samples[name] = 1

        # affected family members should have the mutation (hom ref not allowed)
        if zygosity == "0/0" and family[name] == "1":
            judgement = 0
            break

        # affected family members might be het
        elif zygosity == "0/1" and family[name] == "1":
            judgement = 1
            continue

        # affected family members might be hom alt
        # that"s the major difference to de novo...
        elif zygosity == "1/1" and family[name] == "1":
            judgement = 1
            continue

        # non-affected family members must not have the mutation - hom ref is OK
        elif zygosity == "0/0" and family[name] == "0":
            judgement = 1
            continue

        # non-affected family members must not have the mutation - het is bad
        elif zygosity == "0/1" and family[name] == "0":
            judgement = 0
            break

        # non-affected family members must not have the mutation - hom alt is worst
        elif zygosity == "1/1" and family[name] == "0":
            judgement = 0
            break

        # now a few more complex steps, if genotype is missing (only non-affected individuals should have missing values)
        elif zygosity == "./." and family[name] == "0":

            # if vcf file was not supplemented
            # accept variants which could not be called
            if refcoverage == "." or altcoverage == "." or refcoverage == "" or altcoverage == "":
                judgement = 1
                continue

            # do not do any other judgements, i.e.
            # being tolerant for refcoverage >= 8 and altcoverage == 0, because it could be that the carrier is not yet sick
            # also tolerate low coverage
            judgement = 1
            continue

        # accept some missing values for affected individuals
        elif zygosity == "./." and family[name] == "1":

            # if vcf file was not supplemented by pileup data
            # accept variants which could not be called
            if refcoverage == "." or altcoverage == ".":
                judgement = 1
                continue

            ### do not do any other judgements, i.e.
            ### being tolerant for refcoverage >= 8 and altcoverage == 0, because it could be that the carrier is not yet sick
            ### also tolerate low coverage
            ##judgement = 1
            ##continue

            # now a few more complex steps, if genotype is missing (only non-affected individuals should have missing values)

            refcoverage = float(refcoverage)
            altcoverage = float(altcoverage)
            coverage = refcoverage + altcoverage

            if coverage == 0:
                judgement = 0
                break

            # hom ref
            elif poisson.cdf(float(altcoverage), float(coverage) / 2) <= 0.007 and altcoverage / coverage <= 0.05:
                judgement = 0
                break

            # hom alt
            elif poisson.cdf(float(refcoverage), float(coverage / 2)) <= 0.007  and refcoverage / coverage <= 0.05:
                judgement = 1
                continue

            # het
            elif poisson.cdf(float(altcoverage), float(coverage) / 2) >= 0.007 or altcoverage / coverage >= 0.05:
                judgement = 1
                continue

        else:
            pass

    for vals in check_samples.values():
       if vals == 0:
            judgement = 0
            break

    return judgement

def check_recessive(variant, family, family_type):
    judgement = 0
    check_samples = dict()

    # create data structure for completeness check
    for name in family.keys():
        check_samples[name] = 0

        if "REF." + name in variant.index.tolist() and "ALT." + name in variant.index.tolist():
            zygosity = variant[name]
            if variant["REF." + name] == ".":
                refcoverage = variant["REF." + name].replace(".", "0") # could be numeric or .
            else:
                refcoverage = variant["REF." + name]
            if variant["ALT." + name] == ".":
                altcoverage = variant["ALT." + name].replace(".", "0") # could be numeric or .
            else:
                altcoverage = variant["ALT." + name]
        else:
            #stick with genotype and the others are empty
            zygosity = variant[name]
            refcoverage = "."
            altcoverage = "."

        if name in check_samples:
            check_samples[name] = 1

        # affected individuals have to be homozygous
        if zygosity == "1/1" and family[name] == "1":
            judgement = 1
            continue

        # affected individuals should not be hom ref or het
        elif ( zygosity == "0/0" or zygosity == "0/1" ) and family[name] == "1":
            judgement = 0
            break

        # non-affected individuals might be het
        elif zygosity == "0/1" and family[name] == "0":
            judgement = 1
            continue

        # non-affected individuals might be hom ref, if a family is interrogated
        elif zygosity == "0/0" and family[name] == "0" and family_type == "FAMILY":
            judgement = 1
            continue

        # non-affected individuals in a trio are the parents and have to be het
        elif zygosity == "0/0" and family[name] == "0" and family_type == "TRIO":
            judgement = 0
            break

        # non-affected individuals must not be hom alt
        elif zygosity == "1/1" and family[name] == "0":
            judgement = 0
            break

        # now a few more complex steps, if genotype is missing (only non-affected individuals should have missing values)
        elif zygosity == "./." and family[name] == "0":
            # which chance has the current read distribution to miss out on an alt read
            # e.g. 10ref, 0alt
            # if coverage > 8, the chance is less than 0.5% to miss out on one alt read
            # http://stattrek.com/online-calculator/poisson.aspx

            # if vcf file was not supplemented by pileup data
            # accept variants which could not be called
            if refcoverage == "." or altcoverage == "." or refcoverage == "" or altcoverage == "":
                judgement = 1
                continue

            refcoverage = float(refcoverage)
            altcoverage = float(altcoverage)
            coverage = refcoverage + altcoverage

            if coverage == 0:
                judgement = 0
                break

            # hom ref
            elif poisson.cdf(float(altcoverage), float(coverage) / 2) <= 0.007 and altcoverage / coverage <= 0.05:
                judgement = 0
                break

            # hom alt
            elif poisson.cdf(float(refcoverage), float(coverage / 2)) <= 0.007  and refcoverage / coverage <= 0.05:
                judgement = 0
                break

            # het, which is OK
            elif poisson.cdf(float(altcoverage), float(coverage) / 2) >= 0.007 or altcoverage / coverage >= 0.05:
                judgement = 1
                continue

            # coverage too low?
            else:
                # accept missing values in family interrogations
                if family_type == "FAMILY":
                    judgement = 1
                    continue
                # do not accept missing values in trio setups
                elif family_type == "TRIO":
                    judgement = 0
                    break

                # for security reasons
                judgement = 0
                break

        # do not accept missing values for affected individuals
        # they should be called hom alt by the SNP caller
        elif zygosity == "./." and family[name] == "1":
            judgement = 0
            break

    for vals in check_samples.values():
        if vals == 0:
            judgement = 0
            break

    return judgement


def check_xlinked(variant, family):
    judgement = 0
    check_samples = dict()
    inheritance_logic = dict()

    if not ((variant["CHROM"] == "X") | (variant["CHROM"] == "x") | (variant["CHROM"] == "23")):
        return 0

    # create data structure for completeness check
    for name in family.keys():
        check_samples[name] = 0

        if "REF." + name in variant.index.tolist() and "ALT." + name in variant.index.tolist():
            zygosity = variant[name]
            if variant["REF." + name] == ".":
                refcoverage = variant["REF." + name].replace(".", "0") # could be numeric or .
            else:
                refcoverage = variant["REF." + name]
            if variant["ALT." + name] == ".":
                altcoverage = variant["ALT." + name].replace(".", "0") # could be numeric or .
            else:
                altcoverage = variant["ALT." + name]
        else:
            #stick with genotype and the others are empty
            zygosity = variant[name]
            refcoverage = "."
            altcoverage = "."

        if name in check_samples:
            check_samples[name] = 1

        if family[name] == "0":
            inheritance_logic[name] = zygosity

        # affected individuals have to be homozygous
        if zygosity == "1/1" and family[name] == "1":
            judgement = 1
            continue

        # affected individuals should not be hom ref or het
        elif ( zygosity == "0/0" or zygosity == "0/1" ) and family[name] == "1":
            judgement = 0
            break

        # non-affected individuals might be het
        elif zygosity == "0/1" and family[name] == "0":
            judgement = 1
            continue

        # non-affected individuals might be hom ref
        elif zygosity == "0/0" and family[name] == "0":
            judgement = 1
            continue

        # non-affected individuals must not be hom alt
        elif zygosity == "1/1" and family[name] == "0":
            judgement = 0
            break

        # now a few more complex steps, if genotype is missing (only non-affected individuals should have missing values)
        elif zygosity == "./." and family[name] == "0":
            # which chance has the current read distribution to miss out on an alt read
            # e.g. 10ref, 0alt
            # if coverage > 8, the chance is less than 0.5% to miss out on one alt read
            # http://stattrek.com/online-calculator/poisson.aspx

            # if vcf file was not supplemented by pileup data
            # accept variants which could not be called
            if refcoverage == "." or altcoverage == ".":
                judgement = 1
                continue

            refcoverage = float(refcoverage)
            altcoverage = float(altcoverage)
            coverage = refcoverage + altcoverage

            if coverage == 0:
                judgement = 0
                break

            # hom ref
            elif poisson.cdf(float(altcoverage), float(coverage) / 2) <= 0.007 and altcoverage / coverage <= 0.05:
                inheritance_logic[name] = "0/0"
                judgement = 1
                continue

            # hom alt
            if poisson.cdf(float(refcoverage), float(coverage) / 2) <= 0.007 and refcoverage / coverage <= 0.05:
                inheritance_logic[name] = "1/1"
                judgement = 0
                break

            # het, which is OK
            elif poisson.cdf(float(altcoverage), float(coverage) / 2) >= 0.007 or altcoverage / coverage >= 0.05:
                inheritance_logic[name] = "0/1"
                judgement = 1
                continue

            # coverage too low?
            else:
                judgement = 0
                break

        # do not accept missing values for affected individuals
        # they should be called hom alt by the SNP caller
        elif zygosity == "./." and family[name] == "1":
            judgement = 0
            break

    # sanity check
    het_checker = 0
    hom_checker = 0

    for values in inheritance_logic.values():
        if values == "0/1":
            het_checker = 1
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
