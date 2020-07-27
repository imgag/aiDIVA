import pandas as pd
import numpy as np
import tempfile
import argparse
from itertools import takewhile
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
    vcf_file_to_reformat = open(filepath, "r")

    header_line = ""
    comment_lines = []

    with open(filepath, "r") as temp_vcf:
        for line in temp_vcf:
            if line.strip().startswith("##"):
                comment_lines.append(line.strip())
            if line.strip().startswith("#CHROM"):
                header_line = line.strip()
            else:
                continue

        if header_line == "":
            print("ERROR: The VCF file seems to be corrupted")

    tmp = tempfile.NamedTemporaryFile(mode="w")
    tmp.write(vcf_file_to_reformat.read().replace(r"(\n(?!((((([0-9]{1,2}|[xXyY]{1}|(MT|mt){1})\t)(.+\t){6,}(.+(\n|\Z))))|(#{1,2}.*(\n|\Z))|(\Z))))", ""))

    vcf_header = header_line.strip().split("\t")

    vcf_as_dataframe = pd.read_csv(tmp.name, names=vcf_header, sep="\t", comment="#", low_memory=False)

    vcf_file_to_reformat.close()
    tmp.close()

    vcf_as_dataframe = vcf_as_dataframe.rename(columns={"#CHROM": "CHROM"})
    vcf_as_dataframe = vcf_as_dataframe.drop(columns=["ID", "QUAL", "FILTER"])

    return comment_lines, vcf_as_dataframe
