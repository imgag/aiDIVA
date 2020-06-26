import pandas as pd
import numpy as np
import tempfile
import argparse
from itertools import takewhile
from operator import itemgetter



def split_vcf_file_in_indel_and_snps_set(filepath, filepath_snp, filepath_indel):
    vcf_file_to_reformat = open(filepath, "r")
    outfile_snps = open(filepath_snp, "w")
    outfile_indel = open(filepath_indel, "w")

    # make sure that there are no unwanted linebreaks in the variant entries
    tmp = tempfile.NamedTemporaryFile(mode="w+")
    tmp.write(vcf_file_to_reformat.read().replace(r"(\n(?!((((([0-9]{1,2}|[xXyY]{1}|(MT|mt){1})\t)(.+\t){6,}(.+(\n|\Z))))|(#{1,2}.*(\n|\Z))|(\Z))))", ""))
    tmp.seek(0)

    # extract header from vcf file
    indel_ID = 0
    for line in tmp:
        splitted_line = line.split("\t")
        if line.strip().startswith("##"):
            outfile_snps.write(line)
            outfile_indel.write(line)
            continue
        if line.strip().startswith("#CHROM"):
            outfile_snps.write(line)
            outfile_indel.write(line)
            continue

        # skip empty lines if the VCF file is not correctly formatted (eg. if there are multiple blank lines in the end of the file)
        if line == "\n":
            continue

        # remove variants with multiple alternative alleles reported (to make tests for the master thesis easier)
        # TODO decide how to handle them in general
        if "," in splitted_line[4]:
            print("Variant was removed!")
            print("REASON: Too many alternative alleles reported!")
            continue
        else:
            ref_length = len(splitted_line[3])
            alt_length = max([len(alt) for alt in splitted_line[4].split(",")])

            if (ref_length == 1) & (alt_length == 1):
                outfile_snps.write(line)
            elif (ref_length > 1) | (alt_length > 1):
                indel_ID += 1
                #splitted_line[7] = splitted_line[7].replace("\n", "") + ";indel_ID=indel_" + str(indel_ID) + "\n"
                if splitted_line[7].endswith("\n"):
                    splitted_line[7] = splitted_line[7].replace("\n", "") + ";indel_ID=indel_" + str(indel_ID) + "\n"
                else:
                    splitted_line[7] = splitted_line[7].replace("\n", "") + ";indel_ID=indel_" + str(indel_ID)
                outfile_indel.write("\t".join(splitted_line))
            else:
                print("Something was not rigtht!")

    vcf_file_to_reformat.close()
    outfile_snps.close()
    outfile_indel.close()
    tmp.close()


if __name__=="__main__":
    parser = argparse.ArgumentParser("Script to generate the HPO resources needed in the prioritization step of AIdiva")
    parser.add_argument("--in_file", type=str, dest="in_file", metavar="input.csv", required=True, help="Input file\n")
    parser.add_argument("--snp_file", type=str, dest="snp_file", metavar="snps.csv", required=True, help="File to save the SNP variants\n")
    parser.add_argument("--indel_file", type=str, dest="indel_file", metavar="indels.csv", required=True, help="File to save the InDel variants\n")
    args = parser.parse_args()

    split_vcf_file_in_indel_and_snps_set(args.in_file, args.snp_file, args.indel_file)
