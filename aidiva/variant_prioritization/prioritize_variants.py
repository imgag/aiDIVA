import argparse
import gzip
import json
import logging
import multiprocessing as mp
import networkx as nx
import numpy as np
import os
import pandas as pd
import pysam
import re

from functools import partial
from itertools import combinations

if not __name__=="__main__":
    from . import get_HPO_similarity_score as gs


CODING_VARIANTS = ["stop_gained",
                   "frameshift_variant",
                   "stop_lost",
                   "start_lost",
                   "inframe_insertion",
                   "inframe_deletion",
                   "missense_variant",
                   "protein_altering_variant",
                   "incomplete_terminal_codon_variant",
                   "coding_sequence_variant"]

SYNONYMOUS_VARIANTS = ["synonymous_variant",
                       "start_retained_variant",
                       "stop_retained_variant"]

SPLICE_VARIANTS = ["splice_acceptor_variant",
                   "splice_donor_variant",
                   "splice_donor_5th_base_variant",
                   "splice_region_variant",
                   "splice_donor_region_variant",
                   "splice_polypyrimidine_tract_variant"]

VARIANT_CONSEQUENCES = {"transcript_ablation": 1,
                        "splice_acceptor_variant": 2,
                        "splice_donor_variant": 3,
                        "stop_gained": 4,
                        "frameshift_variant": 5,
                        "stop_lost": 6,
                        "start_lost": 7,
                        "transcript_amplification": 8,
                        "inframe_insertion": 9,
                        "inframe_deletion": 10,
                        "missense_variant": 11,
                        "protein_altering_variant": 12,
                        "splice_region_variant": 13,
                        "splice_donor_5th_base_variant": 14,
                        "splice_donor_region_variant": 15,
                        "splice_polypyrimidine_tract_variant": 16,
                        "incomplete_terminal_codon_variant": 17,
                        "start_retained_variant": 18,
                        "stop_retained_variant": 19,
                        "synonymous_variant": 20,
                        "coding_sequence_variant": 21,
                        "mature_miRNA_variant": 22,
                        "5_prime_UTR_variant": 23,
                        "3_prime_UTR_variant": 24,
                        "non_coding_transcript_exon_variant": 25,
                        "intron_variant": 26,
                        "NMD_transcript_variant": 27,
                        "non_coding_transcript_variant": 28,
                        "upstream_gene_variant": 29,
                        "downstream_gene_variant": 30,
                        "TFBS_ablation": 31,
                        "TFBS_amplification": 32,
                        "TF_binding_site_variant": 33,
                        "regulatory_region_ablation": 34,
                        "regulatory_region_amplification": 35,
                        "feature_elongation": 36,
                        "regulatory_region_variant": 37,
                        "feature_truncation": 38,
                        "intergenic_variant": 39,
                        # use "unknown" consequence as default if new consequence terms are added to the database that are not yet implemented (this prevents the program from exiting with an error)
                        "unknown": 40}


logger = logging.getLogger(__name__)


cadd_identifier = "CADD_PHRED"
duplication_identifier = "segmentDuplication"
repeat_identifier = "simpleRepeat"


def load_from_json(json_file):
    if json_file.endswith(".gz"):
        with gzip.open(json_file, "rt") as json_f:
            loaded_resource = json.load(json_f)

    else:
        with open(json_file, "r") as json_f:
            loaded_resource = json.load(json_f)
    
    return loaded_resource


def parse_ped_file(family_file):
    family_dict = dict()

    if family_file is not None:
        if os.path.isfile(family_file):
            with open(family_file, "r") as family_f:
                for line in family_f:
                    if line == "\n":
                        continue

                    line = line.rstrip()
                    splitline = line.split("\t")

                    if splitline[5] == "2":
                        family_dict[splitline[1]] = 1

                    elif splitline[5] == "1":
                        family_dict[splitline[1]] = 0

                    else:
                        logger.error("There was a problem with the given PED file describing the family relations.")

        else:
            logger.error("The specified family file %s is not a valid file" % (family_file))
            logger.warning("Inheritance assessment will be skipped!")

    return family_dict


def parse_hpo_list(hpo_list_file):
    hpo_query = set()

    if hpo_list_file is not None:
        if os.path.isfile(hpo_list_file):
            with open(hpo_list_file, "r") as hpo_file:
                for line in hpo_file:
                    if line == "\n":
                        continue

                    hpo_term = line.rstrip()
                    hpo_query.add(hpo_term)

            hpo_query = list(hpo_query)
            hpo_query.sort() # makes sure that the hpo terms are ordered (could lead to problems otherwise)

        else:
            hpo_query = hpo_list_file.split(",")
            hpo_query.sort()
            #logger.error("The specified HPO list %s is not a valid file" % (hpo_list_file))

    else:
        logger.warning("HPO score finalization will be skipped!")
    
    return list(hpo_query)


def parse_gene_list(gene_exclusion_file):
    genes2exclude = set()

    if gene_exclusion_file is not None:
        if os.path.isfile(gene_exclusion_file):
            with open(gene_exclusion_file, "r") as exclusion_file:
                for line in exclusion_file:
                    if line.startswith("#"):
                        continue

                    if not line.rstrip():
                        continue

                    gene = line.rstrip().split("\t")[0].upper()
                    genes2exclude.add(gene)

        else:
            logger.error("The specified gene exclusion list %s is not a valid file" % (gene_exclusion_file))
            logger.warning("No genes will be excluded during filtering!")

    return genes2exclude


def get_resource_file(resource_path):
    if resource_path.startswith(".."):
        resource_file = os.path.dirname(__file__) + "/" + resource_path

    else:
        resource_file = resource_path

    return resource_file


def get_feature_completeness(variant, feature_list):
    feature_data = variant[feature_list]
    num_features = len(feature_list)
    num_missing_features = feature_data.isna().sum()
    percentage_missing_features = 0

    if num_missing_features >= 1:
        percentage_missing_features = float(num_missing_features) / float(num_features)

    return percentage_missing_features


def compare_polyphen_and_sift_prediction(variant):
    polyphen_score = variant["PolyPhen"]
    sift_score = variant["SIFT"]

    polyphen_sift_opposed = 0

    if str(polyphen_score) != "nan" and str(polyphen_score) != "" and str(sift_score) != "nan" and str(sift_score) != "":
        polyphen_score = float(polyphen_score)
        sift_score = float(sift_score)

        if polyphen_score < 0.5 and sift_score < 0.5:
            polyphen_sift_opposed = 1

        elif polyphen_score > 0.5 and sift_score > 0.5:
            polyphen_sift_opposed = 1

    return polyphen_sift_opposed


def prioritize_variants(variant_data, internal_parameter_dict, reference, num_cores, build, feature_list, skip_db_check=False, family_file=None, family_type="SINGLE", hpo_list=None, gene_exclusion_list=None):
    # load HPO resources
    hpo_graph_f = get_resource_file(internal_parameter_dict["hpo-graph"])
    hpo_replacement_f = get_resource_file(internal_parameter_dict["hpo2replacement-mapping"])
    gene_2_HPO_f = get_resource_file(internal_parameter_dict["gene2hpo-mapping"])
    hgnc_2_gene_f = get_resource_file(internal_parameter_dict["hgnc2gene-mapping"])
    gene_2_interacting_f = get_resource_file(internal_parameter_dict["gene2interacting-mapping"])
    transcript_length_mapping_f = get_resource_file(internal_parameter_dict["transcript2length-mapping"])

    if not (build == "GRCh37" or build == "GRCh38"):
        logger.error(f"Unrecognized assembly build given: {build}")

    hpo_list_file = hpo_list
    gene_exclusion_file = gene_exclusion_list

    hgnc_2_gene = load_from_json(hgnc_2_gene_f)
    gene_2_interacting = load_from_json(gene_2_interacting_f)
    gene_2_hpo = load_from_json(gene_2_HPO_f)
    hpo_replacement_information = load_from_json(hpo_replacement_f)
    transcript_length_mapping = load_from_json(transcript_length_mapping_f)

    genes2exclude = parse_gene_list(gene_exclusion_file)
    hpo_query = parse_hpo_list(hpo_list_file)
    family = parse_ped_file(family_file)

    hpo_graph = nx.read_gexf(hpo_graph_f)
    information_content_per_node = nx.get_node_attributes(hpo_graph, "IC")
    node_ancestor_mapping = {hpo_term: nx.ancestors(hpo_graph, hpo_term) for hpo_term in hpo_graph}

    variant_data = parallelize_dataframe_processing(variant_data, partial(parallelized_variant_processing, skip_db_check, transcript_length_mapping, family, family_type, genes2exclude, gene_2_hpo, hgnc_2_gene, gene_2_interacting, hpo_graph, hpo_query, information_content_per_node, node_ancestor_mapping, hpo_replacement_information, reference, feature_list), num_cores)
    variant_data = variant_data.sort_values(["FINAL_AIDIVA_SCORE"], ascending=[False])
    variant_data = variant_data.reset_index(drop=True)

    return variant_data


def parallelize_dataframe_processing(variant_data, function, num_cores):
    num_partitions = num_cores * 2

    if len(variant_data) <= num_partitions:
        dataframe_splitted = np.array_split(variant_data, 1)

    else:
        dataframe_splitted = np.array_split(variant_data, num_partitions)

    try:
        pool = mp.Pool(num_cores)
        variant_data = pd.concat(pool.map(function, dataframe_splitted))

    finally:
        pool.close()
        pool.join()

    return variant_data


def parallelized_variant_processing(skip_db_check, transcript_dict, family, family_type, genes2exclude, gene_2_HPO, hgnc_2_gene, gene_2_interacting, HPO_graph, HPO_query, ic_per_nodes, node_ancestor_mapping, hpo_replacement_information, reference, feature_list, variant_data):
    genotype_column = [column for column in variant_data.columns if column.startswith("GT.")]

    if genotype_column:
        variant_data = check_inheritance(variant_data, family_type, family)

    else:
        logger.info(f"Skip inheritance check!")

    variant_data["MISSING_FEATURE_PERCENTAGE"] = variant_data.apply(lambda variant: pd.Series(get_feature_completeness(variant, feature_list)), axis=1)
    variant_data["POLYPHEN_SIFT_OPPOSED"] = variant_data.apply(lambda variant: pd.Series(compare_polyphen_and_sift_prediction(variant)), axis=1)
    
    logger.debug("Investigate Transcript CDS region!")
    variant_data[["CDS_START_PERCENTAGE", "PREDICTED_AIDIVA_SCORE", "AIDIVA_SCORE"]] = variant_data.apply(lambda variant: pd.Series(investigate_transcript_cds_position(variant, transcript_dict)), axis=1)
    
    if not skip_db_check:
        logger.debug("Check databases (ClinVar, HGMD) for known variants!")
        variant_data[["VARIANT_DB_SCORE", "AIDIVA_SCORE"]] = variant_data.apply(lambda variant: pd.Series(check_databases_for_pathogenicity_classification(variant)), axis=1)

    else:
        logger.info(f"Skip variant pathogenicity lookup in existing databases (ClinVar, HGMD)!")

    variant_data[["HPO_RELATEDNESS", "HPO_RELATEDNESS_INTERACTING", "FINAL_AIDIVA_SCORE"]] = variant_data.apply(lambda variant: pd.Series(compute_hpo_relatedness_and_final_score(variant, genes2exclude, gene_2_HPO, hgnc_2_gene, gene_2_interacting, HPO_graph, HPO_query, ic_per_nodes, node_ancestor_mapping, hpo_replacement_information)), axis=1)
    variant_data[["FILTER_PASSED", "FILTER_COMMENT"]] = variant_data.apply(lambda variant: pd.Series(check_filters(variant, genes2exclude, HPO_query, reference)), axis=1)

    return variant_data


def investigate_transcript_cds_position(variant, transcript_dict):
    adjusted_score = np.nan
    start_percentage = np.nan
    predicted_score = variant["AIDIVA_SCORE"]

    if str(variant["CDS_position"]) != "nan" and str(variant["CDS_position"]) != "":
        transcript = str(variant["Feature"])

        if "-" in str(variant["CDS_position"]):
            cds_start = variant["CDS_position"].split("-")[0]

        else:
            cds_start = variant["CDS_position"]

        # filter for the same frameshifts as in the scoring part
        if ("?" not in cds_start) and (("frameshift" in variant["Consequence"]) or (abs(len(variant["REF"]) - len(variant["ALT"])) % 3 != 0)):
            if transcript in transcript_dict.keys():
                cds_length = transcript_dict[transcript]
                if float(cds_length) >= float(cds_start):
                    start_percentage = float(cds_start) / float(cds_length)

                    if 0.95 <= start_percentage < 0.96:
                        adjusted_score = float(predicted_score) - 0.01

                    elif 0.96 <= start_percentage < 0.97:
                        adjusted_score = float(predicted_score) - 0.02

                    elif 0.97 <= start_percentage < 0.98:
                        adjusted_score = float(predicted_score) - 0.03

                    elif 0.98 <= start_percentage < 0.99:
                        adjusted_score = float(predicted_score) - 0.04

                    elif 0.99 <= start_percentage <= 1:
                        adjusted_score = float(predicted_score) - 0.05

                    else:
                        adjusted_score = predicted_score

                else:
                    logger.debug(f"The CDS start position ({cds_start}) was greater than the CDS length ({cds_length})! Skip percentage computation!")
                    start_percentage = np.nan
                    adjusted_score = predicted_score

            else:
                start_percentage = np.nan
                adjusted_score = predicted_score

        else:
            start_percentage = np.nan
            adjusted_score = predicted_score

    else:
        start_percentage = np.nan
        adjusted_score = predicted_score

    return start_percentage, predicted_score, adjusted_score


def check_databases_for_pathogenicity_classification(variant):
    improved_score = np.nan
    database_score = np.nan

    clinvar_classification = str(variant["CLINVAR_DETAILS"]).split("%")[0]
    if "(replaced" in clinvar_classification:
        clinvar_classification = clinvar_classification.split("(replaced")[0]

    elif "/" in clinvar_classification:
        if clinvar_classification.lower() == "benign/likely_benign" or clinvar_classification.lower() == "likely_benign/benign":
            clinvar_classification  = "likely_benign"

        elif clinvar_classification.lower() == "pathogenic/likely_pathogenic" or clinvar_classification.lower() == "likely_pathogenic/pathogenic":
            clinvar_classification = "likely_pathogenic"

        else:
            logger.debug(f"Found unknown ClinVar classification {clinvar_classification}")

    hgmd_classification = str(variant["HGMD_CLASS"])

    if clinvar_classification.lower() == "benign":
        clinvar_score = 0.0

    elif clinvar_classification.lower() == "likely_benign":
        clinvar_score = 0.25

    elif clinvar_classification.lower() == "uncertain_significance":
        clinvar_score = 0.5

    elif clinvar_classification.lower() == "likely_pathogenic":
        clinvar_score = 0.75

    elif clinvar_classification.lower() == "pathogenic":
        clinvar_score = 1.0

    else:
        # ignore missing database scores (if both scores missing use unmodified AIDIVA_SCORE)
        clinvar_score = np.nan

    if hgmd_classification == "DM":
        hgmd_score = 1.0
        #if float(variant["HGMD_RANKSCORE"]) > 0.5:
        #    hgmd_score = 1.0

    elif hgmd_classification == "DM?":
        hgmd_score = 0.75

    #elif hgmd_classification == "DP":
    #    hgmd_score = 0.5

    #elif hgmd_classification == "FP":
    #    hgmd_score = 0.5

    #elif hgmd_classification == "DFP":
    #    hgmd_score = 0.5

    # completely ignore removed entries database entries
    elif hgmd_classification == "R":
        hgmd_score = np.nan

    else:
        # ignore missing database scores (if both scores missing use unmodified AIDIVA_SCORE)
        hgmd_score = np.nan


    # log a warning message if hgmd score and clinvar_score are both not missing but differ
    if (hgmd_score != clinvar_score) and not (np.isnan(hgmd_score) or np.isnan(clinvar_score)):
        logger.warning(f"The scores of the used databases (ClinVar, HGMD) are different: {clinvar_score}, {hgmd_score} (Variant: {variant['#CHROM']}:{variant['POS']}_{variant['REF']}_{variant['ALT']}) \n If the HGMD database was used during the annotation you may want to further investigate this matter, otherwise you can ignore this warning since we just use the predicted value as default value if the entry is missing.")

    if np.isnan(clinvar_score) and not np.isnan(hgmd_score):
        database_score = hgmd_score

    elif np.isnan(hgmd_score) and not np.isnan(clinvar_score):
        database_score = clinvar_score

    elif not np.isnan(hgmd_score) and not np.isnan(clinvar_score):
        database_score = (hgmd_score + clinvar_score) / 2

    else:
        database_score = np.nan

    if not np.isnan(float(variant["AIDIVA_SCORE"])):
        if not np.isnan(database_score):
            # update the predicted score with known information from the databases
            improved_score = 0.8 * float(variant["AIDIVA_SCORE"]) + 0.2 * database_score

        else:
            # use unmodified AIDIVA_SCORE if no database entries were found
            improved_score = float(variant["AIDIVA_SCORE"])

    else:
        improved_score = np.nan

    return [database_score, improved_score]


def compute_hpo_relatedness_and_final_score(variant, genes2exclude, gene_2_HPO, hgnc_2_gene, gene_2_interacting, HPO_graph, HPO_query, ic_per_nodes, node_ancestor_mapping, hpo_replacement_information):
    if HPO_query:
        if np.isnan(variant["AIDIVA_SCORE"]) and ((str(variant["SpliceAI"]) == "nan") or (str(variant["SpliceAI"]) == "")):
            final_score = np.nan
            hpo_relatedness = np.nan
            hpo_relatedness_interacting = np.nan

        else:
            if np.isnan(variant["AIDIVA_SCORE"]):
                if ((str(variant["SpliceAI"]) != "nan") and (str(variant["SpliceAI"]) != "")):
                    pathogenictiy_prediction = float(variant["SpliceAI"])

                else:
                    pathogenictiy_prediction = np.nan

                logger.debug("Using SpliceAI prediction instead of aiDIVA prediction for splicing variant!")

            else:
                pathogenictiy_prediction = float(variant["AIDIVA_SCORE"])

            variant_gene = str(variant["SYMBOL"]).upper()
            hgnc_id = str(variant["HGNC_ID"])

            gene_similarities = []
            gene_similarities_interacting = []

            if variant_gene in gene_2_interacting.keys():
                interacting_genes = [gene_interaction["interacting_gene"] for gene_interaction in gene_2_interacting[variant_gene]]

            else:
                interacting_genes = []

            # we use the hgnc ID to prevent problems if a given gene symbol isn't used anymore
            if (variant_gene not in genes2exclude):
                if variant_gene in gene_2_HPO.keys():
                    gene_HPO_list = gene_2_HPO.get(variant_gene, [])

                else:
                    if (str(hgnc_id) != "nan") and (hgnc_id in hgnc_2_gene.keys()):
                        gene_symbol = hgnc_2_gene[hgnc_id]
                        gene_HPO_list = gene_2_HPO.get(gene_symbol, [])

                    else:
                        logger.debug("The processed gene %s (HGNC:%s) is not found in the current HPO resources" % (variant_gene, hgnc_id))
                        gene_HPO_list = []

                if gene_HPO_list:
                    gene_hpo_similarity = gs.calculate_hpo_set_similarity(HPO_graph, HPO_query, gene_HPO_list, ic_per_nodes, node_ancestor_mapping, hpo_replacement_information)
                    gene_similarities.append(gene_hpo_similarity)

                for interacting_gene in interacting_genes:
                    if (interacting_gene in genes2exclude):
                        continue

                    if interacting_gene in gene_2_HPO.keys():
                        gene_HPO_list = gene_2_HPO.get(interacting_gene, [])

                    else:
                        logger.debug("The processed interacting gene %s is not found in the current HPO resources" % (interacting_gene))
                        gene_HPO_list = []

                    if gene_HPO_list:
                        gene_hpo_similarity = gs.calculate_hpo_set_similarity(HPO_graph, HPO_query, gene_HPO_list, ic_per_nodes, node_ancestor_mapping, hpo_replacement_information)
                        gene_similarities_interacting.append(gene_hpo_similarity)
                    
                hpo_relatedness = max(gene_similarities, default=0.0)
                hpo_relatedness_interacting = max(gene_similarities_interacting, default=0.0)

                # predicted pathogenicity has a higher weight than the HPO relatedness
                # weight optimization suggests the following weights
                final_score = (pathogenictiy_prediction * 0.34 + float(hpo_relatedness) * 0.65 + float(hpo_relatedness_interacting) * 0.01)

            else:
                final_score = np.nan
                hpo_relatedness = np.nan
                hpo_relatedness_interacting = np.nan

    else:
        if np.isnan(variant["AIDIVA_SCORE"]):
            if ((str(variant["SpliceAI"]) == "") or (str(variant["SpliceAI"]) == "nan")):
                pathogenictiy_prediction = np.nan

            else:
                pathogenictiy_prediction = float(variant["SpliceAI"])
                logger.debug("Using SpliceAI prediction instead of aiDIVA prediction for splicing variant!")

        else:
            pathogenictiy_prediction = float(variant["AIDIVA_SCORE"])

        final_score = pathogenictiy_prediction
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

            variant_data["DOMINANT"] = variant_data.apply(lambda variant: check_dominant(variant, family), axis=1)
            variant_data["XLINKED"] = variant_data.apply(lambda variant: check_xlinked(variant, family), axis=1)
            variant_data["RECESSIVE"] = variant_data.apply(lambda variant: check_recessive(variant, family, family_type), axis=1)

            variant_data_grouped = [group for key, group in variant_data.groupby("SYMBOL")]

            affected_child = ""
            parent_1 = ""
            parent_2 = ""

            # TRIO case relies on the assumption that the two parents are healthy and the child is sick
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
                        logger.error("There was a problem processing the loaded PED file")

            if  affected_child and parent_1 and parent_2:
                for group in variant_data_grouped:
                    check_compound(group, affected_child, parent_1, parent_2)

            else:
                logger.error(f"Something went wrong while checking the TRIO for compound heterozygosity!")

            variant_data = pd.concat(variant_data_grouped)

        else:
            logger.error("If family type (TRIO) is used a proper PED file defining the family relations is required!")
            logger.warning("Inheritance assessment will be skipped!")

    elif (family_type == "MULTI") and (family is not None):
        if family is not None:
            variant_data["DOMINANT"] = variant_data.apply(lambda variant: check_dominant(variant, family), axis=1)
            variant_data["XLINKED"] = variant_data.apply(lambda variant: check_xlinked(variant, family), axis=1)
            variant_data["RECESSIVE"] = variant_data.apply(lambda variant: check_recessive(variant, family, family_type), axis=1)

        else:
            logger.error("If family type (MULTI) is used a proper PED file defining the family relations (if known) is required!")
            logger.warning("Inheritance assessment will be skipped!")

    else:
        logger.error("Unsupported family type (%s) was given!" % family_type)
        logger.warning("Inheritance assessment will be skipped!")

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


# Returns longest homopolymer
def longest_homopolymer(sequence):
    if len(sequence) == 0:
        return 0

    runs = ''.join('*' if x == y else ' ' for x, y in zip(sequence, sequence[1:]))
    star_strings = runs.split()
    if len(star_strings) == 0:
        return 1

    return 1 + max(len(stars) for stars in star_strings)


# Returns most frequent base in a string
def frequent_base(sequence):
    # Get the most common base and its percentage of length of the sequence
    nucleotide_list = list(sequence)
    most_common_base = max([nucleotide_list.count(base) for base in set(nucleotide_list)])
    most_common_base_percentage = round(float(most_common_base) / len(sequence), 2)

    return most_common_base_percentage


# Homopolymer filter
def homopolymer_filter(sequence):
    base_type_zero_count = 0  # number of base types (A, C, G, T) not represented in the sequence
    low_complexity_flag = 0  # low complexity = (less than 3 base types represented or homopolymers of length >= 4)
    homopolymer_flag = 0  # homopolymer of length >= 5

    if sequence != '.':

        # Get the longest k-mer
        max_homopolymer_length = longest_homopolymer(sequence)

        # Get the frequency of the homopolymer base
        max_base_freq = frequent_base(sequence)

        # Count base types
        nucleotide_list = list(sequence)
        for base in ['A', 'C', 'G', 'T']:
            if nucleotide_list.count(base) == 0:
                base_type_zero_count += 1

        # Set homopolymer flag
        if (max_homopolymer_length >= 5) or (max_base_freq >= 0.8):
            homopolymer_flag = 1

        # Set low complexity flag
        if (base_type_zero_count > 1) or (max_homopolymer_length >= 4):
            low_complexity_flag = 1

    return homopolymer_flag, low_complexity_flag


def check_filters(variant, genes2exclude, HPO_query, reference):
    variant_genes = re.sub("\(.*?\)", "", str(variant["SYMBOL"]))
    genenames = set(variant_genes.split(";"))

    in_fasta = pysam.FastaFile(reference)

    #consequences = str(variant["Consequence"])
    #found_consequences = consequences.split("&")
    most_severe_consequence = str(variant["MOST_SEVERE_CONSEQUENCE"])

    repeat = str(variant[repeat_identifier])
    filter_comment = ""

    try:
        seg_dup = float(variant[duplication_identifier])

    except Exception as e:
        logger.debug("Use 0.0 for missing segment duplication entries!")
        seg_dup = 0.0
    
    try:
        maf = float(variant["MAX_AF"])

    except Exception as e:
        logger.debug("Allele frequency entry could not be identified, use 0.0 instead!")
        maf = 0.0

    # filter for low confidence regions (these are regions where often false positives were observed)
    if "low_conf_region" in variant["FILTER"]:
        filter_passed = 0 # low confidence region
        filter_comment = "low confidence region"

        return filter_passed, filter_comment

    elif "off-target" in variant["FILTER"]:
        filter_passed = 0 # off target
        filter_comment = "off target"

        return filter_passed, filter_comment

    # exclude gene, if it is in the exclusion list
    if len(genes2exclude & genenames) > 0:
        for gene in genenames:
            if (gene.upper() in genes2exclude):
                filter_passed = 0 # gene in exclusion list
                filter_comment = "gene exclusion"

                return filter_passed, filter_comment

    # let variants with a high FINAL_AIDIVA_SCORE (>=0.7) pass to be more sensitive
    if ((repeat != "NA") and (repeat != "") and (repeat != "nan")) and (float(variant["FINAL_AIDIVA_SCORE"]) < 0.7):
        filter_passed = 0 # tandem repeat
        filter_comment = "tandem repeat"

        return filter_passed, filter_comment

    # let variants with a high FINAL_AIDIVA_SCORE (>=0.7) pass to be more sensitive
    if (seg_dup > 0) and (float(variant["FINAL_AIDIVA_SCORE"]) < 0.7):
        filter_passed = 0 # segmental duplication
        filter_comment = "segmental duplication"

        return filter_passed, filter_comment

    # let variants with a high FINAL_AIDIVA_SCORE (>=0.7) pass to be more sensitive
    if ((len(variant["REF"]) > 1 or len(variant["ALT"]) > 1)) and (float(variant["FINAL_AIDIVA_SCORE"]) < 0.7):
        # make sure to use the correct internal chromsome notation (with chr)
        if "chr" in str(variant["#CHROM"]):
            chrom_id = str(variant["#CHROM"])

        else:
            chrom_id = "chr" + str(variant["#CHROM"])

        # Get sequence context (vicinity) of a variant for homopolymer check (5 bases up- and down-stream)
        # Get fewer bases when variant is at the start or end of the sequence
        if "chr" in str(variant["#CHROM"]):
            chrom_id = str(variant["#CHROM"])

        else:
            chrom_id = "chr" + str(variant["#CHROM"])

        num_bases = 5
        pos_start = max(int(variant["POS"]) - (num_bases + 1), 1)
        pos_end = min(int(variant["POS"]) + num_bases, in_fasta.get_reference_length(chrom_id))

        try:
            sequence_context = in_fasta.fetch(chrom_id, pos_start, pos_end)
            sequence_context = sequence_context.upper()

        except FileNotFoundError:
            sequence_context = '.'

        homopolymer_flag, low_complexity_flag = homopolymer_filter(sequence_context)

        if low_complexity_flag and homopolymer_flag:
            filter_passed = 0
            filter_comment = "homopolymer and low complexity region"

            return filter_passed, filter_comment

        if homopolymer_flag:
            filter_passed = 0
            filter_comment = "homopolymer"

            return filter_passed, filter_comment

        if low_complexity_flag:
            filter_passed = 0
            filter_comment = "low complexity region"

            return filter_passed, filter_comment

    # check if there is something annoted from REPEATMASKER
    # let variants with a high FINAL_AIDIVA_SCORE (>=0.7) pass to be more sensitive
    if (str(variant["REPEATMASKER"]).strip() != "") and (str(variant["REPEATMASKER"]).strip() != ".") and (str(variant["REPEATMASKER"]).strip() != "nan") and (not str(variant['REPEATMASKER']).isspace()) and (float(variant["FINAL_AIDIVA_SCORE"]) < 0.7):
        filter_passed = 0 # masked repeat region
        filter_comment = "masked repeat region"

        return filter_passed, filter_comment

    # MAF threshold could be changed dynamically based on inheritance mode (hom/het)
    if maf <= 0.01 or np.isnan(maf):
        # check if most severe consequence is supported
        if (most_severe_consequence in CODING_VARIANTS) or (most_severe_consequence in SPLICE_VARIANTS):
                if not np.isnan(variant["FINAL_AIDIVA_SCORE"]):
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
                    filter_comment = "missing prediction"

        else:
            filter_passed = 0 # non coding or synonymous

            if most_severe_consequence in SYNONYMOUS_VARIANTS:
                filter_comment = "synonymous variant"

            else:
                filter_comment = "variant type not supported"

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
                logger.debug("Skip variant pair, uncalled genotype in affected sample!")


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
        if name in check_samples:
            check_samples[name] = 1

        else:
            logger.debug(f"It seems that your pedigree is incomplete. The following sample: {name} could not be found in the pedigree!")

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
            logger.debug(f"Skip denovo_check for sample {name}!")
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

        else:
            logger.debug(f"It seems that your pedigree is incomplete. The following sample: {name} could not be found in the pedigree!")

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
            logger.debug(f"Skip denovo_check for sample! Reason one or more sample names where not defined in the pedigree!")
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

        else:
            logger.debug(f"It seems that your pedigree is incomplete. The following sample: {name} could not be found in the pedigree!")

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
        elif zygosity == "0/0" and family[name] == 0 and family_type == "MULTI":
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
            logger.debug(f"Skip denovo_check for sample! Reason one or more sample names where not defined in the pedigree!")
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

    if not ((variant["#CHROM"] == "X") or (variant["#CHROM"] == "x") or (variant["#CHROM"] == "chrX") or (variant["#CHROM"] == "chrx") or (variant["#CHROM"] == "23")):
        return 0

    # create data structure for completeness check
    for name in family.keys():
        check_samples[name] = 0

        zygosity = variant["GT." + name].replace("|", "/")

        if name in check_samples:
            check_samples[name] = 1

        else:
            logger.debug(f"It seems that your pedigree is incomplete. The following sample: {name} could not be found in the pedigree!")

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
            logger.debug(f"Skip denovo_check for sample! Reason one or more sample names where not defined in the pedigree!")
            break

    return judgement


if __name__=="__main__":
    import get_HPO_similarity_score as gs

    parser = argparse.ArgumentParser(description = "Filter variants and finalize the AIDIVA_SCORE based on the given HPO terms (if this information is present)")
    parser.add_argument("--in_file", type=str, dest="in_file", required=True, help="Tab separated input annotated and scored file [required]")
    parser.add_argument("--out_file", type=str, dest="out_filename", required=True, help="Name to save the results [required]")
    parser.add_argument("--family", type=str, dest="family", required=False, help="Tab separated list of samples annotated with affection status.")
    parser.add_argument("--family_type", type=str, choices=["TRIO", "FAMILY", "SINGLE"], dest="family_type", required=False, help="Choose if the data you provide is a trio or a larger family")
    parser.add_argument("--gene_exclusion", type=str, dest="gene_exclusion_list", required=False, help="List of genes that should be excluded in the prioritization")
    parser.add_argument("--hpo_list", type=str, dest="hpo_list", default=None, required=False, help="List of HPO terms that are observed in the patient. These terms are used to adjust the AIDIVA_SCORE\n")
    #parser.add_argument("--hpo_resources", type=str, dest="hpo_resources", default="../../data/", required=True, help="Folder where the HPO resources (HPO_graph,...) are found\n")
    parser.add_argument("--genome_build", type=str, dest="genome_build", default="GRCh38", required=True, help="Version of the genome build to use [GRCh37, GRCH38]\n")
    parser.add_argument("--feature_list", type=str, dest="feature_list", required=True, help="Feature list used in the AIDIVA_SCORE prediction.\n")
    parser.add_argument("--reference", type=str, dest="reference", required=True, help="Path to the refernce genome.\n")
    parser.add_argument("--skip_db_check", action="store_true", required=False, help="Skip db check.\n")
    parser.add_argument("--threads", type=str, dest="threads", default=1, required=False, help="Number of threads to use.\n")
    args = parser.parse_args()

    input_data = pd.read_csv(args.in_file, sep="\t", low_memory=False)

    if args.family:
        family_file = args.family

    else:
        family_file = None

    if args.family_type and args.family_file is not None:
        family_type = args.family_type

    else:
        family_type="SINGLE"

    if args.genome_build:
        genome_build = args.genome_build

    else:
        genome_build = "GRCh38"

    internal_parameter_dict = {"gene2hpo-mapping": "../../data/hpo_resources/gene2hpo.json",
                               "gene2interacting-mapping": "../../data/hpo_resources/gene2interacting.json",
                               "transcript2length-mapping": "../../data/hpo_resources/grch38transcript2length.json",
                               "hgnc2gene-mapping": "../../data/hpo_resources/hgnc2gene.json",
                               "hpo-graph": "../../data/hpo_resources/hpo_graph.gexf",
                               "hpo2replacement-mapping": "../../data/hpo_resources/hpo2replacement.json"}

    if args.threads:
        num_threads = int(args.threads)

    else:
        num_threads = 1

    feature_list = args.feature_list.split(",")

    prioritized_variants = prioritize_variants(input_data, internal_parameter_dict, args.reference, num_threads, genome_build, feature_list, args.skip_db_check, family_file, family_type, args.hpo_list, args.gene_exclusion_list)
    prioritized_variants.to_csv(args.out_filename, sep="\t", index=False)
