import argparse
import logging
import multiprocessing as mp
import numpy as np
import pandas as pd
import tempfile

from functools import partial
from operator import itemgetter


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

USED_INFO_FIELDS = ["INDEL_ID",
                    "CSQ",
                    "FATHMM_XF",
                    "CONDEL",
                    "EIGEN_PHRED",
                    "MutationAssessor",
                    "REVEL",
                    "phyloP_primate",
                    "phyloP_mammal",
                    "phyloP_vertebrate",
                    "phastCons_primate",
                    "phastCons_mammal",
                    "phastCons_vertebrate",
                    "gnomAD_Hom",
                    "gnomAD_AN",
                    "gnomAD_AFR_AF",
                    "gnomAD_AMR_AF",
                    "gnomAD_EAS_AF",
                    "gnomAD_NFE_AF",
                    "gnomAD_SAS_AF",
                    "CAPICE",
                    "CADD",
                    "ADA_SCORE",
                    "RF_SCORE",
                    "oe_lof",
                    "SegDup",
                    "SimpleRepeats",
                    "CLINVAR_DETAILS",
                    #"DETAILS",
                    #"CLASS",
                    "HGMD_CLASS",
                    "HGMD_RANKSCORE",
                    "OMIM",
                    "REPEATMASKER"]

logger = logging.getLogger(__name__)


def reformat_vcf_file_and_read_into_pandas_and_extract_header(filepath):
    header_line = ""
    comment_lines = []

    with open(filepath, "r") as vcf_file_to_reformat, tempfile.NamedTemporaryFile(mode="w+") as tmp:
        # make sure that there are no unwanted linebreaks in the variant entries
        tmp.write(vcf_file_to_reformat.read().replace(r"(\n(?!((((((chr)?[0-9]{1,2}|(chr)?[xXyY]{1}|(chr)?(MT|mt){1})\t)(.+\t){6,}(.+(\n|\Z))))|(#{1,2}.*(\n|\Z))|(\Z))))", ""))
        tmp.seek(0)

        # extract header from vcf file
        for line in tmp:
            if line.strip().startswith("##"):
                comment_lines.append(line.strip())
            if line.strip().startswith("#CHROM"):
                header_line = line.strip()
                comment_lines.append(header_line)
                break # now the variant entries are coming
            else:
                continue

        if header_line == "":
            logger.warning("The VCF file seems to be corrupted")

        vcf_header = header_line.strip().split("\t")
        vcf_as_dataframe = pd.read_csv(tmp.name, names=vcf_header, sep="\t", comment="#", low_memory=False)

    #vcf_as_dataframe = vcf_as_dataframe.rename(columns={"#CHROM": "CHROM"})

    return comment_lines, vcf_as_dataframe


def extract_annotation_header(header):
    annotation_header = [entry.strip().replace("\">", "").split(": ")[1].split("|") for entry in header if entry.startswith("##INFO=<ID=CSQ,")][0]

    return annotation_header


def extract_sample_header(header):
    sample_header = []

    for line in header:
        if line.startswith("##SAMPLE=<"):
            sample_entry = line.strip().replace("##SAMPLE=<", "").replace(">", "")
            sample_id = sample_entry.split(",")[0]
            sample_header.append(sample_id.split("=")[1])

    return sample_header


def extract_columns(cell, process_indel):
    info_fields = str(cell).split(";")

    indel_ID = np.nan
    fathmm_xf = np.nan
    condel = np.nan
    eigen_phred = np.nan
    mutation_assessor = np.nan
    revel = np.nan
    phyloP_primate = np.nan
    phyloP_mammal = np.nan
    phyloP_vertebrate = np.nan
    phastCons_primate = np.nan
    phastCons_mammal = np.nan
    phastCons_vertebrate = np.nan
    gnomAD_hom = np.nan
    gnomAD_an = np.nan
    gnomAD_homAF = np.nan
    gnomAD_afr_af = np.nan
    gnomAD_amr_af = np.nan
    gnomAD_eas_af = np.nan
    gnomAD_nfe_af = np.nan
    gnomAD_sas_af = np.nan
    max_af = np.nan
    capice = np.nan
    cadd = np.nan
    segDup = np.nan
    ada = np.nan
    rf = np.nan
    oe_lof = np.nan
    simpleRepeat = ""
    clinvar_details = ""
    hgmd_class = ""
    hgmd_rankscore = np.nan
    omim_details = ""
    annotation = ""
    repeat_masker = ""

    #details = ""
    #class_orig = ""

    for field in info_fields:
        field_splitted = field.split("=")

        # skip empty INFO field annotations to prevent problems
        if field_splitted[1] != "nan" and field_splitted[1] != "."  and field_splitted[1] != "":

            # skip unsupported INFO fields
            if field_splitted[0] in USED_INFO_FIELDS:

                if field_splitted[0] == "INDEL_ID":
                    indel_ID = field_splitted[1]

                elif field_splitted[0] == "CSQ":
                    annotation = field_splitted[1]

                # TODO: add SIFT and PolyPhen also here (when vep is not used anymore)

                elif field_splitted[0] == "FATHMM_XF":
                    if "&" in field_splitted[1]:
                        fathmm_xf = max([float(value) for value in field_splitted[1].split("&") if (value != "." and value != "nan" and value !="")], default=np.nan)
                    else:
                        fathmm_xf = float(field_splitted[1])

                elif field_splitted[0] == "CONDEL":
                    if "&" in field_splitted[1]:
                        condel = max([float(value) for value in field_splitted[1].split("&") if (value != "." and value != "nan" and value !="")], default=np.nan)
                    else:
                        condel = float(field_splitted[1])

                elif field_splitted[0] == "EIGEN_PHRED":
                    if "&" in field_splitted[1]:
                        eigen_phred = max([float(value) for value in field_splitted[1].split("&") if (value != "." and value != "nan" and value !="")], default=np.nan)
                    else:
                        eigen_phred = float(field_splitted[1])

                elif field_splitted[0] == "MutationAssessor":
                    if "&" in field_splitted[1]:
                        mutation_assessor = max([float(value) for value in field_splitted[1].split("&") if (value != "." and value != "nan" and value !="")], default=np.nan)
                    else:
                        mutation_assessor = float(field_splitted[1])

                elif field_splitted[0] == "REVEL":
                    if "&" in field_splitted[1]:
                        revel = max([float(value) for value in field_splitted[1].split("&") if (value != "." and value != "nan" and value !="")], default=np.nan)
                    else:
                        revel = float(field_splitted[1])

                elif field_splitted[0] == "phyloP_primate":
                    if "&" in field_splitted[1]:
                        phyloP_primate = max([float(value) for value in field_splitted[1].split("&") if (value != "." and value != "nan" and value !="")], default=np.nan)
                    else:
                        phyloP_primate = float(field_splitted[1])

                elif field_splitted[0] == "phyloP_mammal":
                    if "&" in field_splitted[1]:
                        phyloP_mammal = max([float(value) for value in field_splitted[1].split("&") if (value != "." and value != "nan" and value !="")], default=np.nan)
                    else:
                        phyloP_mammal = float(field_splitted[1])

                elif field_splitted[0] == "phyloP_vertebrate":
                    if "&" in field_splitted[1]:
                        phyloP_vertebrate = max([float(value) for value in field_splitted[1].split("&") if (value != "." and value != "nan" and value !="")], default=np.nan)
                    else:
                        phyloP_vertebrate = float(field_splitted[1])

                elif field_splitted[0] == "phastCons_primate":
                    if "&" in field_splitted[1]:
                        phastCons_primate = max([float(value) for value in field_splitted[1].split("&") if (value != "." and value != "nan" and value !="")], default=np.nan)
                    else:
                        phastCons_primate = float(field_splitted[1])

                elif field_splitted[0] == "phastCons_mammal":
                    if "&" in field_splitted[1]:
                        phastCons_mammal = max([float(value) for value in field_splitted[1].split("&") if (value != "." and value != "nan" and value !="")], default=np.nan)
                    else:
                        phastCons_mammal = float(field_splitted[1])

                elif field_splitted[0] == "phastCons_vertebrate":
                    if "&" in field_splitted[1]:
                        phastCons_vertebrate = max([float(value) for value in field_splitted[1].split("&") if (value != "." and value != "nan" and value !="")], default=np.nan)
                    else:
                        phastCons_vertebrate = float(field_splitted[1])

                elif field_splitted[0] == "gnomAD_Hom":
                    if "&" in field_splitted[1]:
                        gnomAD_hom = float(field_splitted[1].split("&")[0])
                    else:
                        gnomAD_hom = float(field_splitted[1])

                elif field_splitted[0] == "gnomAD_AN":
                    if "&" in field_splitted[1]:
                        gnomAD_an = float(field_splitted[1].split("&")[0])
                    else:
                        gnomAD_an = float(field_splitted[1])

                elif field_splitted[0] == "gnomAD_AFR_AF":
                    if "&" in field_splitted[1]:
                        gnomAD_afr_af = max([float(value) for value in field_splitted[1].split("&") if (value != "." and value != "nan" and value !="")], default=np.nan)
                    else:
                        gnomAD_afr_af = float(field_splitted[1])

                elif field_splitted[0] == "gnomAD_AMR_AF":
                    if "&" in field_splitted[1]:
                        gnomAD_amr_af = max([float(value) for value in field_splitted[1].split("&") if (value != "." and value != "nan" and value !="")], default=np.nan)
                    else:
                        gnomAD_amr_af = float(field_splitted[1])

                elif field_splitted[0] == "gnomAD_EAS_AF":
                    if "&" in field_splitted[1]:
                        gnomAD_eas_af = max([float(value) for value in field_splitted[1].split("&") if (value != "." and value != "nan" and value !="")], default=np.nan)
                    else:
                        gnomAD_eas_af = float(field_splitted[1])

                elif field_splitted[0] == "gnomAD_NFE_AF":
                    if "&" in field_splitted[1]:
                        gnomAD_nfe_af = max([float(value) for value in field_splitted[1].split("&") if (value != "." and value != "nan" and value !="")], default=np.nan)
                    else:
                        gnomAD_nfe_af = float(field_splitted[1])

                elif field_splitted[0] == "gnomAD_SAS_AF":
                    if "&" in field_splitted[1]:
                        gnomAD_sas_af = max([float(value) for value in field_splitted[1].split("&") if (value != "." and value != "nan" and value !="")], default=np.nan)
                    else:
                        gnomAD_sas_af = float(field_splitted[1])

                elif field_splitted[0] == "CAPICE":
                    if "&" in field_splitted[1]:
                        capice = max([float(value) for value in field_splitted[1].split("&") if (value != "." and value != "nan" and value !="")], default=np.nan)
                    else:
                        capice = float(field_splitted[1])

                elif field_splitted[0] == "CADD":
                    if "&" in field_splitted[1]:
                        cadd = max([float(value) for value in field_splitted[1].split("&") if (value != "." and value != "nan" and value !="")], default=np.nan)
                    else:
                        cadd = float(field_splitted[1])

                elif field_splitted[0] == "ADA_SCORE":
                    if "&" in field_splitted[1]:
                        ada = max([float(value) for value in field_splitted[1].split("&") if (value != "." and value != "nan" and value !="")], default=np.nan)
                    else:
                        ada = float(field_splitted[1])

                elif field_splitted[0] == "RF_SCORE":
                    if "&" in field_splitted[1]:
                        rf = max([float(value) for value in field_splitted[1].split("&") if (value != "." and value != "nan" and value !="")], default=np.nan)
                    else:
                        rf = float(field_splitted[1])

                elif field_splitted[0] == "oe_lof":
                    oe_lof = min([float(value) for value in field_splitted[1].split("&")  if (value != "." and value != "nan" and value !="")], default=np.nan)

                elif field_splitted[0] == "SegDup":
                    segDup = max([float(value) for value in field_splitted[1].split("&")  if (value != "." and value != "nan" and value !="")], default=np.nan)

                elif field_splitted[0] == "SimpleRepeats":
                    simpleRepeat = field_splitted[1]

                elif field_splitted[0] == "CLINVAR_DETAILS":
                    clinvar_details = field_splitted[1]

                #elif field_splitted[0] == "DETAILS":
                #    details = field_splitted[1]

                #elif field_splitted[0] == "CLASS":
                #    class_orig = field_splitted[1]

                elif field_splitted[0] == "HGMD_CLASS":
                    hgmd_class = field_splitted[1]

                elif field_splitted[0] == "HGMD_RANKSCORE":
                    if "&" in field_splitted[1]:
                        hgmd_rankscore = max([float(value) for value in field_splitted[1].split("&") if (value != "." and value != "nan" and value !="")], default=np.nan)
                    else:
                        hgmd_rankscore = float(field_splitted[1])

                elif field_splitted[0] == "OMIM":
                    omim_details = field_splitted[1]

                elif field_splitted[0] == "REPEATMASKER":
                    repeat_masker = field_splitted[1]

                else:
                    logger.error(f"Could not recognize INFO field {field}")

            else:
                logger.debug(f"Skip unused INFO field value {field}!")

        else:
            logger.debug(f"Skip empty (NaNs are handled as empty) INFO field {field}")

    if (gnomAD_hom > 0.0) and (gnomAD_an > 0.0):
        gnomAD_homAF = gnomAD_hom / gnomAD_an

    max_af = max([gnomAD_afr_af, gnomAD_amr_af, gnomAD_eas_af, gnomAD_nfe_af, gnomAD_sas_af], default=np.nan)

    if process_indel:
        extracted_columns = [indel_ID, annotation, fathmm_xf, condel, eigen_phred, mutation_assessor, revel, phyloP_primate, phyloP_mammal, phyloP_vertebrate, phastCons_primate, phastCons_mammal, phastCons_vertebrate, gnomAD_homAF, gnomAD_afr_af, gnomAD_amr_af, gnomAD_eas_af, gnomAD_nfe_af, gnomAD_sas_af, max_af, capice, cadd, oe_lof, segDup, simpleRepeat, ada, rf, repeat_masker, clinvar_details, hgmd_class, hgmd_rankscore, omim_details] #, details, class_orig] #, abb_score]

    else:
        extracted_columns = [annotation, fathmm_xf, condel, eigen_phred, mutation_assessor, revel, phyloP_primate, phyloP_mammal, phyloP_vertebrate, phastCons_primate, phastCons_mammal, phastCons_vertebrate, gnomAD_homAF, gnomAD_afr_af, gnomAD_amr_af, gnomAD_eas_af, gnomAD_nfe_af, gnomAD_sas_af, max_af, capice, cadd, oe_lof, segDup, simpleRepeat, ada, rf, repeat_masker, clinvar_details, hgmd_class, hgmd_rankscore, omim_details] #, details, class_orig, abb_score]

    return extracted_columns


def extract_vep_annotation(cell, annotation_header):
    annotation_fields = str(cell["CSQ"]).split(",")
    new_cols = []

    # choose the most severe annotation variant
    # if new consequence terms were added to the database that are not yet handled from AIdiva use default consequence "unknown" with lowest severity value
    consequences = [min([VARIANT_CONSEQUENCES.get(consequence) if consequence in VARIANT_CONSEQUENCES.keys() else VARIANT_CONSEQUENCES.get("unknown") for consequence in field.split("|")[annotation_header.index("Consequence")].split("&")]) for field in annotation_fields]

    target_index = min(enumerate(consequences), key=itemgetter(1))[0]
    new_cols = annotation_fields[target_index].strip().split("|")

    return new_cols


def extract_sample_information(row, sample, sample_header=None):
    if str(row["FORMAT"]) == "MULTI":
        if sample_header is not None:
            sample_dict = {}
            sample_information = []
            for sample_entry in row[sample].split(","):
                splitted_entry = sample_entry.split("=")
                sample_dict[splitted_entry[0]] = splitted_entry[1]

            for sample_id in sample_header:
                splitted_sample_information = sample_dict[sample_id].split("|")
                if splitted_sample_information[0] == "wt":
                    sample_information.append("0/0")
                elif splitted_sample_information[0] == "hom":
                    sample_information.append("1/1")
                elif splitted_sample_information[0] == "het":
                    sample_information.append("0/1")
                else:
                    logger.warning("Genotype not recognized! (%s)" % (splitted_sample_information[0]))

                sample_information.append(splitted_sample_information[1])
                sample_information.append(splitted_sample_information[2])
                sample_information.append(sample_id + "=" + sample_dict[sample_id])
        else:
            logger.error("Format is MULTI but no sample_header is given!")
    else:
        format_entries = str(row["FORMAT"]).strip().split(":")
        sample_fields = str(row[sample + ".full"]).strip().split(":")

        if len(format_entries) != len(sample_fields):
            num_missing_entries = abs(len(format_entries) - len(sample_fields))
            for i in range(num_missing_entries):
                sample_fields.append(".")

        if "GT" in format_entries:
            sample_gt_information = sample_fields[format_entries.index("GT")]
        else:
            sample_gt_information = "./."

        if "DP" in format_entries:
            sample_dp_information = sample_fields[format_entries.index("DP")]
        else:
            sample_dp_information = "."

        if "AD" in format_entries:
            sample_ref_information = sample_fields[format_entries.index("AD")].split(",")[0]
            sample_alt_information = sample_fields[format_entries.index("AD")].split(",")[1]
        else:
            sample_ref_information = "."
            sample_alt_information = "."

        if "GQ" in format_entries:
            sample_gq_information = sample_fields[format_entries.index("GQ")]
        else:
            sample_gq_information = "."

        if (sample_ref_information != ".") and (sample_alt_information != "."):
            divisor = (int(sample_ref_information) + int(sample_alt_information))
            if divisor == 0:
                sample_af_information = 0
            else:
                sample_af_information = (int(sample_alt_information) / divisor)

        else:
            sample_af_information = "."

        sample_information = [sample_gt_information, sample_dp_information, sample_ref_information, sample_alt_information, sample_af_information, sample_gq_information]

    return sample_information


# Find and return common prefix of two strings
def find_common_prefix(string_a, string_b):
    prefix_end_position = 0

    length_string_a = len(string_a)
    length_string_b = len(string_b)
    max_length_prefix = min(length_string_a, length_string_b)

    while prefix_end_position < max_length_prefix and string_a[prefix_end_position] == string_b[prefix_end_position]:
        prefix_end_position += 1

    shared_prefix = string_a[:prefix_end_position]

    return shared_prefix


# Find and return common suffix of two strings
def find_common_suffix(string_a, string_b):
    suffix_start_position = 0

    length_string_a = len(string_a)
    length_string_b = len(string_b)
    max_length_suffix = min(length_string_a, length_string_b)

    while suffix_start_position < max_length_suffix and string_a[length_string_a-suffix_start_position-1] == string_b[length_string_b-suffix_start_position-1]:
        suffix_start_position += 1

    if suffix_start_position == 0:
        return ""

    common_suffix = string_a[suffix_start_position:]

    return common_suffix


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

    # remove common suffix
    suffix_length = len(find_common_suffix(ref, alt));
    if suffix_length > 0:
        ref = ref[:-suffix_length]
        alt = alt[:-suffix_length]

    # remove common prefix
    prefix_length = len(find_common_prefix(ref, alt))
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


def add_INFO_fields_to_dataframe(process_indel, expanded_indel, vcf_as_dataframe):
    indel_annotation_columns = ["INDEL_ID", "CSQ", "FATHMM_XF", "CONDEL", "EIGEN_PHRED", "MutationAssessor", "REVEL", "phyloP_primate", "phyloP_mammal", "phyloP_vertebrate", "phastCons_primate", "phastCons_mammal", "phastCons_vertebrate", "homAF", "gnomAD_AFR_AF", "gnomAD_AMR_AF", "gnomAD_EAS_AF", "gnomAD_NFE_AF", "gnomAD_SAS_AF", "MAX_AF", "CAPICE", "CADD_PHRED", "oe_lof", "segmentDuplication", "simpleRepeat", "ada_score", "rf_score", "REPEATMASKER", "CLINVAR_DETAILS", "HGMD_CLASS", "HGMD_RANKSCORE", "OMIM"]
    snp_annotation_columns = ["CSQ", "FATHMM_XF", "CONDEL", "EIGEN_PHRED", "MutationAssessor", "REVEL", "phyloP_primate", "phyloP_mammal", "phyloP_vertebrate", "phastCons_primate", "phastCons_mammal", "phastCons_vertebrate", "homAF", "gnomAD_AFR_AF", "gnomAD_AMR_AF", "gnomAD_EAS_AF", "gnomAD_NFE_AF", "gnomAD_SAS_AF", "MAX_AF", "CAPICE", "CADD_PHRED", "oe_lof", "segmentDuplication", "simpleRepeat", "ada_score", "rf_score", "REPEATMASKER", "CLINVAR_DETAILS", "HGMD_CLASS", "HGMD_RANKSCORE", "OMIM"]

    if process_indel:
        if not expanded_indel:
            vcf_as_dataframe["GSVAR_VARIANT"] = vcf_as_dataframe.apply(lambda row: convert_variant_representation(row), axis=1)
            vcf_as_dataframe[indel_annotation_columns] = vcf_as_dataframe["INFO"].apply(lambda x: pd.Series(extract_columns(x, process_indel)))

        else:
            vcf_as_dataframe[indel_annotation_columns] = vcf_as_dataframe["INFO"].apply(lambda x: pd.Series(extract_columns(x, process_indel)))

    else:
        vcf_as_dataframe["GSVAR_VARIANT"] = vcf_as_dataframe.apply(lambda row: convert_variant_representation(row), axis=1)
        vcf_as_dataframe[snp_annotation_columns] = vcf_as_dataframe["INFO"].apply(lambda x: pd.Series(extract_columns(x, process_indel)))

    vcf_as_dataframe = vcf_as_dataframe.drop(columns=["INFO"])

    return vcf_as_dataframe


def add_VEP_annotation_to_dataframe(annotation_header, vcf_as_dataframe):
    vcf_as_dataframe[annotation_header] = vcf_as_dataframe.apply(lambda x: pd.Series(extract_vep_annotation(x, annotation_header)), axis=1)
    vcf_as_dataframe = vcf_as_dataframe.drop(columns=["CSQ"])

    return vcf_as_dataframe


def add_sample_information_to_dataframe(sample_ids, sample_header, vcf_as_dataframe):
    for sample in sample_ids:
        if (sample == "trio") or (sample == "multi"):
            sample_header_multi = []
            for sample_id in sample_header:
                sample_header_multi.append("GT." + sample_id)
                sample_header_multi.append("DP." + sample_id)
                sample_header_multi.append("ALT." + sample_id)
                sample_header_multi.append(sample_id + ".full")
            vcf_as_dataframe[sample_header_multi] = vcf_as_dataframe.apply(lambda x: pd.Series(extract_sample_information(x, sample, sample_header)), axis=1)

        else:
            vcf_as_dataframe.rename(columns={sample: sample + ".full"}, inplace=True)
            sample_header_default = ["GT." + sample, "DP." + sample, "REF." + sample, "ALT." + sample, "AF." + sample, "GQ." + sample]
            vcf_as_dataframe[sample_header_default] = vcf_as_dataframe.apply(lambda x: pd.Series(extract_sample_information(x, sample)), axis=1)
            vcf_as_dataframe = vcf_as_dataframe.drop(columns=[sample + ".full"])

    vcf_as_dataframe = vcf_as_dataframe.drop(columns=["FORMAT"])

    return vcf_as_dataframe


def convert_vcf_to_pandas_dataframe(input_file, process_indel, expanded_indel, num_cores):
    header, vcf_as_dataframe = reformat_vcf_file_and_read_into_pandas_and_extract_header(input_file)

    logger.debug("Convert VCF file")

    sample_ids = []
    # FORMAT column has index 8 (counted from 0) and sample columns follow afterwards (sample names are unique)
    # Check if FORMAT column exists
    if len(vcf_as_dataframe.columns) > 8:
        for i in range(9, len(vcf_as_dataframe.columns)):
            sample_ids.append(vcf_as_dataframe.columns[i])

        sample_header = extract_sample_header(header)

    annotation_header = extract_annotation_header(header)

    if not vcf_as_dataframe.empty:
        vcf_as_dataframe = parallelize_dataframe_processing(vcf_as_dataframe, partial(add_INFO_fields_to_dataframe, process_indel, expanded_indel), num_cores)
        vcf_as_dataframe = parallelize_dataframe_processing(vcf_as_dataframe, partial(add_VEP_annotation_to_dataframe, annotation_header), num_cores)

        if len(vcf_as_dataframe.columns) > 8:
            if "FORMAT" in vcf_as_dataframe.columns:
                vcf_as_dataframe = parallelize_dataframe_processing(vcf_as_dataframe, partial(add_sample_information_to_dataframe, sample_ids, sample_header), num_cores)

            else:
                # This warning is always triggered when the expanded indel vcf file is processed. If it is only triggered once in this case it can be ignored.
                logger.warning("It seems that your VCF file does contain sample information but no FORMAT description!")

        else:
            logger.error("MISSING SAMPLE INFORMATION!")

        # replace empty strings or only spaces with NaN
        vcf_as_dataframe = vcf_as_dataframe.replace(r"^\s*$", np.nan, regex=True)

    else:
        logger.error("The given VCF file is empty!")

    return vcf_as_dataframe


def parallelize_dataframe_processing(vcf_as_dataframe, function, num_cores):
    num_partitions = num_cores * 2

    if len(vcf_as_dataframe) <= num_partitions:
        dataframe_splitted = np.array_split(vcf_as_dataframe, 1)

    else:
        dataframe_splitted = np.array_split(vcf_as_dataframe, num_partitions)

    try:
        pool = mp.Pool(num_cores)
        vcf_as_dataframe = pd.concat(pool.map(function, dataframe_splitted))

    finally:
        pool.close()
        pool.join()

    return vcf_as_dataframe


def write_vcf_to_csv(vcf_as_dataframe, out_file):
    vcf_as_dataframe.to_csv(out_file, sep="\t", encoding="utf-8", index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--in_data", type=str, dest="in_data", metavar="input.vcf", required=True, help="VCF file to convert file\n")
    parser.add_argument("--out_data", type=str, dest="out_data", metavar="output.csv", required=True, help="CSV file containing the converted VCF file\n")
    parser.add_argument("--indel", action="store_true", required=False, help="Flag to indicate whether the file to convert consists of indel variants or not.\n")
    parser.add_argument("--expanded", action="store_true", required=False, help="Flag to indicate whether the file to convert consists of expanded indel variants.\n")
    parser.add_argument("--threads", type=int, dest="threads", metavar="1", required=False, help="Number of threads to use.")
    args = parser.parse_args()

    if args.threads is not None:
        num_cores = int(args.threads)

    else:
        num_cores = 1

    vcf_as_dataframe = convert_vcf_to_pandas_dataframe(args.in_data, args.indel, args.expanded, num_cores)
    write_vcf_to_csv(vcf_as_dataframe, args.out_data)
