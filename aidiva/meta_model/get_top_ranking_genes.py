import argparse
import gzip
import json
import logging
import os
import pandas as pd
import random
import re


logger = logging.getLogger(__name__)


# Find and return shared prefix of two strings
def find_shared_prefix(string_a, string_b):
    prefix_end_position = 0

    length_string_a = len(string_a)
    length_string_b = len(string_b)
    max_length_prefix = min(length_string_a, length_string_b)

    while prefix_end_position < max_length_prefix and string_a[prefix_end_position] == string_b[prefix_end_position]:
        prefix_end_position += 1

    shared_prefix = string_a[:prefix_end_position]

    return shared_prefix


# Find and return shared suffix of two strings
def find_shared_suffix(string_a, string_b):
    suffix_start_position = 0

    length_string_a = len(string_a)
    length_string_b = len(string_b)
    max_length_suffix = min(length_string_a, length_string_b)

    while suffix_start_position < max_length_suffix and string_a[length_string_a-suffix_start_position-1] == string_b[length_string_b-suffix_start_position-1]:
        suffix_start_position += 1

    if suffix_start_position == 0:
        return ""

    shared_suffix = string_a[suffix_start_position:]

    return shared_suffix


def convert_variant_representation(row):
    chrom = row["#CHROM"]
    start_position = row["POS"]
    ref = row["REF"].upper()
    alt = row["ALT"].upper()

    # remove common first base
    if ref != "" and alt != "" and ref[0] == alt[0]:
        ref = ref[1:]
        alt = alt[1:]
        start_position +=1;

    # remove shared suffix
    suffix_length = len(find_shared_suffix(ref, alt));
    if suffix_length > 0:
        ref = ref[:-suffix_length]
        alt = alt[:-suffix_length]

    # remove shared prefix
    prefix_length = len(find_shared_prefix(ref, alt))
    if prefix_length > 0:
        ref = ref[prefix_length:]
        alt = alt[prefix_length:]
        start_position += prefix_length

    # determine start and end
    end_position = start_position
    ref_length = len(ref)
    alt_length = len(alt)

    if alt_length == 1 and ref_length == 1: # SNV
        # nothing to do for SNVs
        pass

    elif ref_length == 0: # insertion
        ref = "-"
        # change insertions from before the coordinate to after the coordinate!!!!
        start_position -= 1
        end_position -= 1

    elif alt_length == 0: # deletion
        end_position = start_position + ref_length - 1
        alt="-"

    elif alt_length >= 1 and ref_length > 1: #complex indel
        end_position = start_position + ref_length - 1

    normalized_variant = f"{chrom}:{start_position}-{end_position}_{ref}>{alt}"

    return normalized_variant


def extract_gene_info(row, sample_id, CONSEQUENCE_MAPPING):
    current_gene = row["SYMBOL"]
    current_ref = row["REF"]
    current_alt = row["ALT"]
    current_consequence = row["MOST_SEVERE_CONSEQUENCE"]
    current_rank = int(row["AIDIVA_RANK"])
    current_score = float(row["FINAL_AIDIVA_SCORE"])

    # convert variant representation to match variant in GSvar file
    if "GSVAR_VARIANT" in list(row.index):
        current_gsvar_variant = row["GSVAR_VARIANT"]

    else:
        current_gsvar_variant = convert_variant_representation(row)

    # fixed for single case samples or if multisample it should now default to the index patient
    if f"GT_{sample_id}" in list(row.index):
        current_genotype = row[f"GT_{sample_id}"]

    else:
        current_genotype = "nan"

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
