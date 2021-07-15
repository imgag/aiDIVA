import pandas as pd
import numpy as np
import multiprocessing as mp
import tempfile
import argparse
from functools import partial
from operator import itemgetter


variant_consequences = {"transcript_ablation": 1,
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
                        "incomplete_terminal_codon_variant": 14,
                        "start_retained_variant": 15,
                        "stop_retained_variant": 16,
                        "synonymous_variant": 17,
                        "coding_sequence_variant": 18,
                        "mature_miRNA_variant": 19,
                        "5_prime_UTR_variant": 20,
                        "3_prime_UTR_variant": 21,
                        "non_coding_transcript_exon_variant": 22,
                        "intron_variant": 23,
                        "NMD_transcript_variant": 24,
                        "non_coding_transcript_variant": 25,
                        "upstream_gene_variant": 26,
                        "downstream_gene_variant": 27,
                        "TFBS_ablation": 28,
                        "TFBS_amplification": 29,
                        "TF_binding_site_variant": 30,
                        "regulatory_region_ablation": 31,
                        "regulatory_region_amplification": 32,
                        "feature_elongation": 33,
                        "regulatory_region_variant": 34,
                        "feature_truncation": 35,
                        "intergenic_variant": 36}


def reformat_vcf_file_and_read_into_pandas_and_extract_header(filepath):
    header_line = ""
    comment_lines = []

    vcf_file_to_reformat = open(filepath, "r")

    # make sure that there are no unwanted linebreaks in the variant entries
    tmp = tempfile.NamedTemporaryFile(mode="w+")
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
        print("ERROR: The VCF file seems to be corrupted")

    vcf_header = header_line.strip().split("\t")
    vcf_as_dataframe = pd.read_csv(tmp.name, names=vcf_header, sep="\t", comment="#", low_memory=False)

    vcf_file_to_reformat.close()
    tmp.close()

    vcf_as_dataframe = vcf_as_dataframe.rename(columns={"#CHROM": "CHROM"})
    #vcf_as_dataframe = vcf_as_dataframe.drop(columns=["ID", "QUAL", "FILTER"])

    return comment_lines, vcf_as_dataframe


def extract_annotation_header(header):
    annotation_header = [entry.strip().replace("\">", "").split(": ")[1].split("|") for entry in header if entry.startswith("##INFO=<ID=CSQ,")][0]

    return annotation_header


def extract_columns(cell, process_indel):
    info_fields = str(cell).split(";")

    indel_ID = np.nan
    fathmm_xf = np.nan
    condel = np.nan
    eigen_phred = np.nan
    mutation_assessor = np.nan
    gnomAD_hom = np.nan
    gnomAD_an = np.nan
    gnomAD_homAF = np.nan
    capice = np.nan
    annotation = ""

    for field in info_fields:
        if process_indel:
            if field.startswith("INDEL_ID"):
                indel_ID = field.split("INDEL_ID=")[1]
        if field.startswith("CSQ="):
            annotation = field.split("CSQ=")[1]
        if field.startswith("FATHMM_XF="):
            if field.split("FATHMM_XF=")[1] != "nan":
                fathmm_xf = field.split("FATHMM_XF=")[1]
        if field.startswith("CONDEL="):
            if field.split("CONDEL=")[1] != "nan":
                condel = field.split("CONDEL=")[1]
        if field.startswith("EIGEN_PHRED="):
            if field.split("EIGEN_PHRED=")[1] != "nan":
                eigen_phred = field.split("EIGEN_PHRED=")[1]
        if field.startswith("MutationAssessor="):
            if field.split("MutationAssessor=")[1] != "nan":
                mutation_assessor = field.split("MutationAssessor=")[1]
        if field.startswith("gnomAD_Hom"):
            if field.split("gnomAD_Hom=")[1] != "nan":
                gnomAD_hom = float(field.split("gnomAD_Hom=")[1])
        if field.startswith("gnomAD_AN"):
            if field.split("gnomAD_AN=")[1] != "nan":
                gnomAD_an = float(field.split("gnomAD_AN=")[1])
        if field.startswith("CAPICE"):
            if field.split("CAPICE=")[1] != "nan":
                capice = float(field.split("CAPICE=")[1])

    if (gnomAD_hom > 0.0) and (gnomAD_an > 0.0):
        gnomAD_homAF = gnomAD_hom / gnomAD_an

    if process_indel:
        extracted_columns = [indel_ID, annotation, fathmm_xf, condel, eigen_phred, mutation_assessor, gnomAD_homAF, capice]
    else:
       extracted_columns = [annotation, fathmm_xf, condel, eigen_phred, mutation_assessor, gnomAD_homAF, capice]

    return extracted_columns


def extract_vep_annotation(cell, annotation_header):
    annotation_fields = str(cell["CSQ"]).split(",")
    new_cols = []
    consequences = []

    # take the most severe annotation variant
    for field in annotation_fields:
        consequences.append(min([variant_consequences.get(x) for x in field.split("|")[annotation_header.index("Consequence")].split("&")]))

    target_index = min(enumerate(consequences), key=itemgetter(1))[0]
    new_cols = annotation_fields[target_index].strip().split("|")

    return new_cols


def extract_sample_information(row, sample):
    sample_header = str(row["FORMAT"]).strip().split(":")
    sample_fields = str(row[sample + ".full"]).strip().split(":")

    if len(sample_header) != len(sample_fields):
        num_missing_entries = abs(len(sample_header) - len(sample_fields))
        for i in range(num_missing_entries):
            sample_fields.append(".")

    if len(sample_header) != len(sample_fields):
        num_missing_entries = abs(len(sample_header) - len(sample_fields))
        for i in range(num_missing_entries):
            sample_fields.append(".")

    if "GT" in sample_header:
        sample_gt_information = sample_fields[sample_header.index("GT")]
    else:
        sample_gt_information = "./."

    if "DP" in sample_header:
        sample_dp_information = sample_fields[sample_header.index("DP")]
    else:
        sample_dp_information = "."

    if "AD" in sample_header:
        sample_ref_information = sample_fields[sample_header.index("AD")].split(",")[0]
        sample_alt_information = sample_fields[sample_header.index("AD")].split(",")[1]
    else:
        sample_ref_information = "."
        sample_alt_information = "."

    if "GQ" in sample_header:
        sample_gq_information = sample_fields[sample_header.index("GQ")]
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


def add_INFO_fields_to_dataframe(process_indel, vcf_as_dataframe):
    if process_indel:
        vcf_as_dataframe[["INDEL_ID", "CSQ", "FATHMM_XF", "CONDEL", "EIGEN_PHRED", "MutationAssessor", "homAF", "CAPICE"]] = vcf_as_dataframe["INFO"].apply(lambda x: pd.Series(extract_columns(x, process_indel)))
    else:
        vcf_as_dataframe[["CSQ", "FATHMM_XF", "CONDEL", "EIGEN_PHRED", "MutationAssessor", "homAF", "CAPICE"]] = vcf_as_dataframe["INFO"].apply(lambda x: pd.Series(extract_columns(x, process_indel)))

    vcf_as_dataframe = vcf_as_dataframe.drop(columns=["INFO"])

    return vcf_as_dataframe


def add_VEP_annotation_to_dataframe(annotation_header, vcf_as_dataframe):
    vcf_as_dataframe[annotation_header] = vcf_as_dataframe.apply(lambda x: pd.Series(extract_vep_annotation(x, annotation_header)), axis=1)
    vcf_as_dataframe = vcf_as_dataframe.drop(columns=["CSQ"])

    return vcf_as_dataframe


def add_sample_information_to_dataframe(sample_ids, vcf_as_dataframe):
    for sample in sample_ids:
        vcf_as_dataframe.rename(columns={sample: sample + ".full"}, inplace=True)
        sample_header = ["GT." + sample, "DP." + sample, "REF." + sample, "ALT." + sample, "AF." + sample, "GQ." + sample]
        vcf_as_dataframe[sample_header] = vcf_as_dataframe.apply(lambda x: pd.Series(extract_sample_information(x, sample)), axis=1)
        vcf_as_dataframe = vcf_as_dataframe.drop(columns=[sample + ".full"])

    vcf_as_dataframe = vcf_as_dataframe.drop(columns=["FORMAT"])

    return vcf_as_dataframe


def convert_vcf_to_pandas_dataframe(input_file, process_indel, num_cores):
    header, vcf_as_dataframe = reformat_vcf_file_and_read_into_pandas_and_extract_header(input_file)

    sample_ids = []
    # FORMAT column has index 8 (counted from 0) and sample columns follow afterwards (sample names are unique)
    ## TODO: add condition to check if FORMAT column exists
    for i in range(9, len(vcf_as_dataframe.columns)):
        sample_ids.append(vcf_as_dataframe.columns[i])

    annotation_header = extract_annotation_header(header)

    if not vcf_as_dataframe.empty:
        vcf_as_dataframe = parallelize_dataframe_processing(vcf_as_dataframe, partial(add_INFO_fields_to_dataframe, process_indel), num_cores)
        vcf_as_dataframe = parallelize_dataframe_processing(vcf_as_dataframe, partial(add_VEP_annotation_to_dataframe, annotation_header), num_cores)

        if len(vcf_as_dataframe.columns) > 8:
            if "FORMAT" in vcf_as_dataframe.columns:
                vcf_as_dataframe = parallelize_dataframe_processing(vcf_as_dataframe, partial(add_sample_information_to_dataframe, sample_ids), num_cores)
            else:
                # This warning is always triggered when the expanded indel vcf file is processed. If it is only triggered once in this case it can be ignored.
                print("WARNING: It seems that your VCF file does contain sample information but not FORMAT description!")
        else:
            print("MISSING SAMPLE INFORMATION!")

        # replace empty strings or only spaces with NaN
        vcf_as_dataframe = vcf_as_dataframe.replace(r"^\s*$", np.nan, regex=True)
    else:
        print("WARNING: The given VCF file is empty!")

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
    parser.add_argument("--threads", type=int, dest="threads", metavar="1", nargs="?", const=1, required=True, help="Number of threads to use.")
    args = parser.parse_args()

    vcf_as_dataframe = convert_vcf_to_pandas_dataframe(args.in_data, args.indel, args.threads)
    write_vcf_to_csv(vcf_as_dataframe, args.out_data)
