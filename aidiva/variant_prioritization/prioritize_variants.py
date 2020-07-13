import argparse
import networkx
import os
import pandas as pd
import pickle
import re

if not __name__=="__main__":
    from . import get_HPO_similarity_score as gs

from itertools import combinations
from operator import itemgetter
from scipy.stats import poisson


variant_consequences = {"transcript_ablation": "non_exonic",
                        "splice_acceptor_variant": "exonic;splicing",
                        "splice_donor_variant": "exonic;splicing",
                        "stop_gained": "exonic",
                        "frameshift_variant": "exonic",
                        "stop_lost": "exonic",
                        "start_lost": "exonic",
                        "transcript_amplification": "non_exonic",
                        "inframe_insertion": "exonic",
                        "inframe_deletion": "exonic",
                        "missense_variant": "exonic",
                        "protein_altering_variant": "exonic",
                        "splice_region_variant": "splicing", # maybe not
                        "incomplete_terminal_codon_variant": "exonic",
                        "start_retained_variant": "exonic", # maybe not
                        "stop_retained_variant": "exonic", # maybe not
                        "synonymous_variant": "exonic",
                        "coding_sequence_variant": "exonic",
                        "mature_miRNA_variant": "non_exonic",
                        "5_prime_UTR_variant": "UTR_exonic",
                        "3_prime_UTR_variant": "UTR_exonic",
                        "non_coding_transcript_exon_variant": "non_exonic",
                        "intron_variant": "intronic",
                        "NMD_transcript_variant": "non_exonic",
                        "non_coding_transcript_variant": "non_exonic",
                        "upstream_gene_variant": "non_exonic",
                        "downstream_gene_variant": "non_exonic",
                        "TFBS_ablation": "non_exonic",
                        "TFBS_amplification": "non_exonic",
                        "TF_binding_site_variant": "non_exonic",
                        "regulatory_region_ablation": "non_exonic",
                        "regulatory_region_amplification": "non_exonic",
                        "feature_elongation": "non_exonic",
                        "regulatory_region_variant": "non_exonic",
                        "feature_truncation": "non_exonic",
                        "intergenic_variant": "non_exonic"}


#def prioritize_variants(in_file, out_file, filtered_file, fam_file, inheritance, family_type, white_list, gene_exclusion):
def prioritize_variants(variant_data, family_file, family_type, hpo_resources_folder, hpo_list, gene_exclusion_list):
    #variant_data = pd.read_csv(in_file, sep="\t", low_memory=False)

    # read the gene exclusion list
    # an empty set will not filter out anything, if gene exclusion list is not provided
    genes2exclude = set()
    genes_known = set()

    ### HPO files load
    gene_2_HPO_f = hpo_resources_folder + "gene2hpo.pkl"
    HPO_graph_file = hpo_resources_folder + "hpo_graph.pkl"
    hpo_dict_file = hpo_resources_folder + "hpo2gene.pkl"
    hpo_list_file = hpo_list
    gene_exclusion_file = gene_exclusion_list

    gene_2_HPO = pickle.load(open(gene_2_HPO_f, "rb"))
    HPO_graph = pickle.load(open(HPO_graph_file, "rb"))

    query_dist = 0

    if gene_exclusion_file:
        gene_exclusion = open(gene_exclusion_file, "r")
        for gene in gene_exclusion:
            gene = gene.rstrip()
            genes2exclude.add(gene)
        gene_exclusion.close()

    #fill HPO query terms
    if hpo_list_file == None:
        hpo_list = "None"
        HPO_query = list()
    else:
        HPO_query = set()
        if os.path.isfile(hpo_list_file):
            with open(hpo_list_file, "r") as w:
                HPO_dict = pickle.load(open(hpo_dict_file,"rb"))
                HPO_query = list()
                for line in w:
                    HPO_term = line.rstrip("\n")
                    try:
                        HPO_query.append(HPO_term)
                    except:
                        print("%s not found in database" % (HPO_term))
            HPO_query= list(set(HPO_query))
            query_dist = gs.precompute_query_distances(HPO_graph, HPO_query, 0)
        else:
            print("The specified HPO list %s is not a valid file" % (hpo_list))
            print("eDiVA will proceed as without any HPO list")

    # read family relationships
    # TODO change to ped file
    family = dict()
    with open(family_file, "r") as fam_file:
        for line in fam_file:
            if line.startswith("sample"):
                continue
            line = line.rstrip("\n")
            splitline = line.split("\t")
            family[splitline[0]] = splitline[1]

    family_type = family_type
    cadd_identifier = "CADD_PHRED"
    duplication_identifier = "segmentDuplication"
    repeat_identifier = "simpleRepeat"
    print(family)

    variant_data[["HPO_RELATEDNESS", "FINAL_AIDIVA_SCORE"]] = variant_data.apply(lambda variant: pd.Series(compute_hpo_relatedness_and_final_score(variant, genes2exclude, HPO_graph, gene_2_HPO, HPO_query, query_dist)), axis=1)
    variant_data["COMPOUND"] = 0
    variant_data[["RECESSIVE", "DOMINANT_DENOVO", "DOMINANT_INHERITED", "XLINKED", "FILTER_PASSED"]] = variant_data.apply(lambda variant: pd.Series(check_inheritance_and_filters(variant, genes2exclude, HPO_query, family, family_type, cadd_identifier, duplication_identifier, repeat_identifier)), axis=1)

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
    variant_data.sort_values(["FINAL_AIDIVA_SCORE"], ascending=[False])
    variant_data.reset_index(inplace=True, drop=True)
    #variant_data.to_csv(out_file, sep="\t", encoding="utf-8", index=False)
    #variant_data[variant_data["FILTER_PASSED"] == 1].to_csv(filtered_out_file, sep="\t", encoding="utf-8", index=False)

    return variant_data


def compute_hpo_relatedness_and_final_score(variant, genes2exclude, HPO_graph, gene_2_HPO, HPO_query, query_dist):
    genecolumn = re.sub("\(.*?\)", "", str(variant["SYMBOL"]))
    genenames = set(genecolumn.split(";"))

    gene_distances = []
    processed_HPO_genes = dict()
    for gene_id in genenames:
        #if not gene_id in genes2exclude: # TODO: decide where to handle the geneexclusion list
        if gene_id in processed_HPO_genes.keys():
            gene_distances.append(processed_HPO_genes[gene_id])
        else:
            #process ex novo
            #get HPOs related to the gene
            gene_HPO_list = gs.extract_HPO_related_to_gene(gene_2_HPO, gene_id)
            (g_dist,query_dist) = gs.list_distance(HPO_graph, HPO_query, gene_HPO_list, query_dist)
            gene_distances.append(g_dist)
            processed_HPO_genes[gene_id] = g_dist
    hpo_relatedness = str(max(gene_distances))

    if float(variant["AIDIVA_SCORE"] != -1.0):
        final_score = str((float(variant["AIDIVA_SCORE"]) + float(hpo_relatedness)) / 2)
    else:
        final_score = float(variant["AIDIVA_SCORE"])

    return [hpo_relatedness, final_score]


def check_inheritance_and_filters(variant, genes2exclude, HPO_list, family, family_type, cadd_identifier, duplication_identifier, repeat_identifier):
    genecolumn = re.sub("\(.*?\)", "", str(variant["SYMBOL"]))
    genenames = set(genecolumn.split(";"))

    consequences = str(variant["Consequence"])
    found_consequences = [variant_consequences[consequence] for consequence in consequences.split("&")]
    seg_dup = float(variant[duplication_identifier])
    tandem = str(variant[repeat_identifier])
    cadd = float(variant[cadd_identifier])

    try:
        #maf = max(float(variant["AA_AF"]), float(variant["AFR_AF"]), float(variant["AMR_AF"]), float(variant["EA_AF"]), float(variant["EAS_AF"]), float(variant["EUR_AF"]), float(variant["SAS_AF"]), float(variant["gnomAD_AFR_AF"]), float(variant["gnomAD_AMR_AF"]), float(variant["gnomAD_ASJ_AF"]), float(variant["gnomAD_EAS_AF"]), float(variant["gnomAD_FIN_AF"]), float(variant["gnomAD_NFE_AF"]), float(variant["gnomAD_OTH_AF"]), float(variant["gnomAD_SAS_AF"]))
        maf = variant["MAX_AF"]
    except Exception as e:
        print("Allele frequency could not be identified, use 0.0 instead")
        maf = 0.0

    dominant_denovo = check_denovo(variant, family)

    ## TODO: do we need a check for affected family members?
    dominant_inherited = check_dominant(variant, family)
    if (variant["CHROM"] == "X") | (variant["CHROM"] == "x") | (variant["CHROM"] == "23"):
        xlinked = check_xlinked(variant, family)
    else:
        xlinked = 0
    recessive = check_recessive(variant, family, family_type)

    # exclude gene, if it is on the exclusion list
    if len(genes2exclude & genenames) > 0:
        for gene in genenames:
            if gene in genes2exclude:
                filter_passed = 0 # gene in exclusion list
                return [recessive, dominant_denovo, dominant_inherited, xlinked, filter_passed]

    if (tandem != "NA") & (tandem != "") & (tandem != "nan"):
        filter_passed = 0 # tandem repeat
        return [recessive, dominant_denovo, dominant_inherited, xlinked, filter_passed]

    ## TODO: filter later compound only less than 0.01
    if float(maf) <= 0.02:
        if (("exonic" in found_consequences) | ("splicing" in found_consequences) | ("exonic;splicing" in found_consequences)):
            if not (("synonymous_variant" in consequences.split("&")) & ("unknown" != consequences) & ("UNKNOWN" != consequences)):
                if (seg_dup == 0):
                    filter_passed = 1
                    if (len(HPO_list) > 1) & ("NONE" not in HPO_list):
                        if float(variant["HPO_RELATEDNESS"]) > 0:
                            filter_passed = 1
                        else:
                            filter_passed = 0 # no relation to reported HPO terms
                # e.g. intronic variants fitting the criteria
                else:
                    filter_passed = 0 # segment duplication
            else:
                filter_passed = 0 # synonymous variant  or unknown effect
        else:
            filter_passed = 0 # not exonic
    else:
        filter_passed = 0 # allele frequency to high

    return [recessive, dominant_denovo, dominant_inherited, xlinked, filter_passed]


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
    parser.add_argument("--out_filename", type=str, dest="out_filename", required=True, help="Name to save the results [required]")
    parser.add_argument("--family", type=str, dest="family", required=True, help="Tab separated list of samples annotated with affection status. [required]")
    parser.add_argument("--family_type", type=str, choices=["TRIO", "FAMILY", "SINGLE"], dest="family_type", required=True, help="Choose if the data you provide is a trio or a larger family [required]")
    parser.add_argument("--gene_exclusion_list", type=str, dest="gene_exclusion_list", required=False, help="List of genes that should be excluded in the prioritization")
    parser.add_argument("--hpo_list", type=str, dest="hpo_list", default=None, required=False, help="List of HPO terms that are observed in the patient. These terms are used to adjust the AIDIVA_SCORE\n")
    parser.add_argument("--hpo_resources", type=str, dest="hpo_resources", default="../../data/", required=True, help="Folder where the HPO resources (HPO_graph,...) are found\n")
    args = parser.parse_args()

    input_data = pd.read_csv(args.in_file, sep="\t", low_memory=False)

    prioritized_variants = prioritize_variants(input_data, args.family, args.family_type, args.hpo_resources, args.hpo_list, args.gene_exclusion_list)
    prioritized_variants.to_csv(args.out_filename, sep="\t", index=False)
