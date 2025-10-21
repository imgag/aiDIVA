import argparse
import gzip
import json
import logging
import os
import pandas as pd
import random
import re


logger = logging.getLogger(__name__)


def extract_gene_info(row, sample_id, CONSEQUENCE_MAPPING):
    current_gene = row["SYMBOL"]
    current_ref = row["REF"]
    current_alt = row["ALT"]
    current_consequence = row["MOST_SEVERE_CONSEQUENCE"]
    current_rank = int(row["AIDIVA_RANK"])
    current_score = float(row["FINAL_AIDIVA_SCORE"])
    current_gsvar_variant = row["GSVAR_VARIANT"]

    # fixed for single case samples or if multisample it should now default to the index patient
    current_genotype = row[f"GT_{sample_id}"]

    if current_genotype == "1/1" or current_genotype == "1|1":
        current_genotype = "hom"

    elif current_genotype == "0/1" or current_genotype == "0|1" or current_genotype == "1/0" or current_genotype == "1|0":
        current_genotype = "het"

    else:
        logger.error(f"ERROR: {sample_id} unsupported genotype -> {current_genotype}")
        current_genotype = "nan"

    # workaround to handle protein_altering_variant and coding_sequence_variant
    if current_consequence == "protein_altering_variant" or current_consequence == "coding_sequence_variant":
        if len(current_ref) > 1 and len(current_alt) == 1:
            if abs(len(current_ref) - len(current_alt)) % 3 == 0:
                current_consequence_type = CONSEQUENCE_MAPPING["inframe_deletion"]

            else:
                current_consequence_type = CONSEQUENCE_MAPPING["frameshift_variant"]

        if len(current_ref) == 1 and len(current_alt) > 1:
            if abs(len(current_ref) - len(current_alt)) % 3 == 0:
                current_consequence_type = CONSEQUENCE_MAPPING["inframe_insertion"]

            else:
                current_consequence_type = CONSEQUENCE_MAPPING["frameshift_variant"]

    else:
        if current_consequence in CONSEQUENCE_MAPPING.keys():
            current_consequence_type = CONSEQUENCE_MAPPING[current_consequence]

        else:
            current_consequence_type = "unspecified variant"
            logger.warning(f"Unsupported variant: {current_consequence}, {row['GSVAR_VARIANT']}")

    result_string = f"{current_gene}(Consequence: {current_consequence}, Type: {current_consequence_type}, Rank: {int(current_rank)}, Score: {current_score}, Genotype: {current_genotype}, Variant: {current_gsvar_variant})"

    return result_string


def extract_top_ranking_entries_random_forest_based(sample_id, result_data_random_forest, maximum_rank, CONSTANT_DICTIONARY):
    # get constants
    VARIANT_CONSEQUENCES = CONSTANT_DICTIONARY["VARIANT_CONSEQUENCES"]
    CONSEQUENCE_MAPPING = CONSTANT_DICTIONARY["CONSEQUENCE_MAPPING"]
    top_ranking_results = result_data_random_forest[result_data_random_forest["AIDIVA_RANK"] <= maximum_rank].copy(deep=True)
    top_ranking_results["extracted_gene_info"] = top_ranking_results.apply(lambda row: extract_gene_info(row, sample_id, CONSEQUENCE_MAPPING), axis=1)
    top_ranking_entries = ";".join(top_ranking_results["extracted_gene_info"].to_list())

    logger.debug("Top ranking results random forest:\n", top_ranking_results)
    logger.debug("Top ranking entries random forest:\n", top_ranking_entries)

    return top_ranking_entries


def get_most_severe_consequence(current_consequence, VARIANT_CONSEQUENCES):
    found_consequences = current_consequence.split("&")
    # choose most severe consequence if overlapping consequences are present
    most_severe_consequence = found_consequences[0]

    for consequence in found_consequences:
        if VARIANT_CONSEQUENCES[consequence] < VARIANT_CONSEQUENCES[most_severe_consequence]:
            most_severe_consequence = consequence

        else:
            most_severe_consequence = found_consequences[0]

    return most_severe_consequence


def extract_gene_info_gsvar(row, sample_id, VARIANT_CONSEQUENCES, CONSEQUENCE_MAPPING):
    current_ref = row["ref"]
    current_alt = row["obs"]
    current_rank = row["GSvar_rank"]
    current_score = row["GSvar_score"]
    current_genotype = row[sample_id]
    current_gene_info = row["coding_and_splicing"]

    # extract gene name and consequence from gene_info
    processed_genes = []
    processed_consequences = []
    processed_consequence_types = []

    for gene in current_gene_info.split(","):
        current_gene = gene.split(":")[0]
        current_consequence = gene.split(":")[2]

        if "_variant_variant" in current_consequence:
            current_consequence = current_consequence.replace("_variant_variant", "_variant")

        # only coding variants are considered
        if any(supported_consequence in current_consequence for supported_consequence in CONSEQUENCE_MAPPING.keys()):
            ## TODO: handle mutliple overlapping gene/consequence information for a single variant

            # workaround to handle protein_altering_variant and coding_sequence_variant
            if current_consequence == "protein_altering_variant" or current_consequence == "coding_sequence_variant":
                if len(current_ref) > 1 and len(current_alt) == 1:
                    if abs(len(current_ref) - len(current_alt)) % 3 == 0:
                        current_consequence_type = CONSEQUENCE_MAPPING["inframe_deletion"]

                    else:
                        current_consequence_type = CONSEQUENCE_MAPPING["frameshift_variant"]

                if len(current_ref) == 1 and len(current_alt) > 1:
                    if abs(len(current_ref) - len(current_alt)) % 3 == 0:
                        current_consequence_type = CONSEQUENCE_MAPPING["inframe_insertion"]

                    else:
                        current_consequence_type = CONSEQUENCE_MAPPING["frameshift_variant"]

            else:
                if "&" in current_consequence:
                    current_consequence = get_most_severe_consequence(current_consequence, VARIANT_CONSEQUENCES)

                current_consequence_type = CONSEQUENCE_MAPPING[current_consequence]

            processed_genes.append(current_gene)
            processed_consequences.append(current_consequence)
            processed_consequence_types.append(current_consequence_type)

        else:
            logger.warning(f"Skip gene with unsupported variant consequence! Gene: {current_gene}, Consequence: {current_consequence}")

    if not processed_genes:
        logger.warning(f"Use of gene with unsupported variant type! Gene: {current_gene}, Consequence: {current_consequence}")

        if "&" in current_consequence:
            current_consequence = get_most_severe_consequence(current_consequence)

        processed_genes.append(current_gene)
        processed_consequences.append(current_consequence)
        processed_consequence_types.append("unspecified variant")

    result_gene = processed_genes[0]
    result_consequence = processed_consequences[0]
    resutl_consequence_type = processed_consequence_types[0]

    for i in range(len(processed_genes)):
        if VARIANT_CONSEQUENCES[processed_consequences[i]] < VARIANT_CONSEQUENCES[result_consequence]:
            result_gene = processed_genes[i]
            result_consequence = processed_consequences[i]
            resutl_consequence_type = processed_consequence_types[i]

        else:
            continue

    result_string = f"{result_gene}(Consequence: {result_consequence}, Type: {resutl_consequence_type}, Rank: {int(current_rank)}, Score: {current_score}, Genotype: {current_genotype}, Variant: {row['#chr']}:{row['start']}-{row['end']}_{current_ref}>{current_alt})"
    logger.debug("Result string top ranking entries evidence based:", result_string)

    return result_string


def extract_top_ranking_entries_evidence_based(sample_id, in_data_evidence, maximum_rank, CONSTANT_DICTIONARY):
    # get constants
    VARIANT_CONSEQUENCES = CONSTANT_DICTIONARY["VARIANT_CONSEQUENCES"]
    CONSEQUENCE_MAPPING = CONSTANT_DICTIONARY["CONSEQUENCE_MAPPING"]

    gsvar_header = extract_gsvar_header(in_data_evidence)
    result_data_evidence = pd.read_csv(in_data_evidence, comment="#", names=gsvar_header, sep="\t", low_memory=False, on_bad_lines="warn")
    top_ranking_results_evidence = result_data_evidence[result_data_evidence["GSvar_rank"] <= maximum_rank].copy(deep=True)
    top_ranking_results_evidence["extracted_gene_info"] = top_ranking_results_evidence.apply(lambda row: extract_gene_info_gsvar(row, sample_id, VARIANT_CONSEQUENCES, CONSEQUENCE_MAPPING), axis=1)
    top_ranking_entries_evidence = ";".join(top_ranking_results_evidence["extracted_gene_info"].to_list())

    logger.debug("Top ranking results evidence based:\n", top_ranking_results_evidence)
    logger.debug("Top ranking entries evidence based:\n", top_ranking_entries_evidence)

    return top_ranking_entries_evidence


def extract_gsvar_header(filepath):
    header_line = ""

    with open(filepath, "r") as vcf_file_to_reformat: #, tempfile.NamedTemporaryFile(mode="w+") as tmp:
        # extract header from gsvar file
        for line in vcf_file_to_reformat:
            if line.startswith("##"):
                continue

            elif line.strip().startswith("#chr"):
                header_line = line.strip()
                break # now the variant entries are coming

            else:
                continue

        if header_line == "":
            logger.warning("The GSvar file seems to be corrupted")

    gsvar_header = header_line.split("\t")

    return gsvar_header


# for debugging
def main(infile, sample_id, rank, evidence):
    if evidence:
        gsvar_header = extract_gsvar_header(infile)
        in_data = pd.read_csv(infile, comment="#", names=gsvar_header, sep="\t", low_memory=False, on_bad_lines="warn" )
        top_ranking_genes = extract_top_ranking_entries_evidence_based(sample_id, in_data, rank)
        print("Final top ranking genes evidence:", top_ranking_genes)

    else:
        in_data = pd.read_csv(infile, sep="\t", low_memory=False)
        top_ranking_genes = extract_top_ranking_entries_random_forest_based(sample_id, in_data, rank)
        print("Final top ranking genes random forest:", top_ranking_genes)


## The possibility to directly run the script is mainly meant for testing und debugging
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Create table with prompts for chatGPT")
    parser.add_argument("--in_file", type=str, dest="in_file", required=True, help="Tab separated input file containing sample ids to choose [required]")
    parser.add_argument("--sample_id", type=str, dest="sample_id", required=True, help="Sample id [required]")
    parser.add_argument("--rank", type=str, dest="rank", required=True, help="Maximum rank to include in the variant/gene list [required]")
    parser.add_argument("--evidence", action="store_true", required=False, help="Flag if data is evidence based.")
    args = parser.parse_args()

    main(args.in_file, args.sample_id, int(args.rank), args.evidence)
