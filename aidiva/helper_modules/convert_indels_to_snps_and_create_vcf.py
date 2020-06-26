import pandas as pd
import numpy as np
import argparse
import random
from Bio import SeqIO


def write_header(out_file):
    out_file.write("##fileformat=VCFv4.1\n")
    out_file.write("##FILTER=<ID=PASS,Description=\"All filters passed\">\n")
    out_file.write("##FILTER=<ID=LowQual,Description=\"Low quality\">\n")
    out_file.write("##ALT=<ID=NON_REF,Description=\"Represents any possible alternative allele at this location\">\n")
    out_file.write("##INFO=<ID=indel_ID,Number=.,Type=String,Description=\"Unique ID to identify the to which this SNP belongs.\">\n")
    out_file.write("##contig=<ID=1,length=249250621,assembly=hg19>\n")
    out_file.write("##contig=<ID=2,length=243199373,assembly=hg19>\n")
    out_file.write("##contig=<ID=3,length=198022430,assembly=hg19>\n")
    out_file.write("##contig=<ID=4,length=191154276,assembly=hg19>\n")
    out_file.write("##contig=<ID=5,length=180915260,assembly=hg19>\n")
    out_file.write("##contig=<ID=6,length=171115067,assembly=hg19>\n")
    out_file.write("##contig=<ID=7,length=159138663,assembly=hg19>\n")
    out_file.write("##contig=<ID=8,length=146364022,assembly=hg19>\n")
    out_file.write("##contig=<ID=9,length=141213431,assembly=hg19>\n")
    out_file.write("##contig=<ID=10,length=135534747,assembly=hg19>\n")
    out_file.write("##contig=<ID=11,length=135006516,assembly=hg19>\n")
    out_file.write("##contig=<ID=12,length=133851895,assembly=hg19>\n")
    out_file.write("##contig=<ID=13,length=115169878,assembly=hg19>\n")
    out_file.write("##contig=<ID=14,length=107349540,assembly=hg19>\n")
    out_file.write("##contig=<ID=15,length=102531392,assembly=hg19>\n")
    out_file.write("##contig=<ID=16,length=90354753,assembly=hg19>\n")
    out_file.write("##contig=<ID=17,length=81195210,assembly=hg19>\n")
    out_file.write("##contig=<ID=18,length=78077248,assembly=hg19>\n")
    out_file.write("##contig=<ID=19,length=59128983,assembly=hg19>\n")
    out_file.write("##contig=<ID=20,length=63025520,assembly=hg19>\n")
    out_file.write("##contig=<ID=21,length=48129895,assembly=hg19>\n")
    out_file.write("##contig=<ID=22,length=51304566,assembly=hg19>\n")
    out_file.write("##contig=<ID=X,length=155270560,assembly=hg19>\n")
    out_file.write("##contig=<ID=Y,length=59373566,assembly=hg19>\n")
    out_file.write("##contig=<ID=MT,length=16571,assembly=hg19>\n")
    out_file.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")


# TODO include FORMAT and SAMPLE information
def write_data_information_to_file(input_data, outfile, ref_folder, header):
    #print(input_data.columns)
    data_grouped = [group for key, group in input_data.groupby("CHROM")]

    random.seed(14038)

    for line in header:
        if line.startswith("#CHROM"):
            outfile.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        else:
            outfile.write(line)


    for group in data_grouped:
        #print(group)
        ref_seq = str(SeqIO.read(ref_folder + "Homo_sapiens.GRCh37.dna.chromosome." + str(group["CHROM"].iloc[0]) + ".fa", "fasta").seq)
        for row in group.itertuples():
            window_start = int(row.POS) - 3
            window_end = int(row.POS) + len(row.REF) + 2
            extended_ref_seq = ref_seq[window_start:window_end]

            for i in range(abs(window_end-window_start)):
                alt_variant = ""
                if (extended_ref_seq[i] == "A") | (extended_ref_seq[i] == "T"):
                    alt_variant = random.choice(["G", "C"])
                elif (extended_ref_seq[i] == "G") | (extended_ref_seq[i] == "C"):
                    alt_variant = random.choice(["A", "T"])
                else:
                    print("ERROR: Something went wrong!")

                outfile.write(str(row.CHROM).strip() + "\t" +
                            str(window_start + i + 1).strip() + "\t" +
                            "." + "\t" +
                            str(extended_ref_seq[i]).strip() + "\t" +
                            str(alt_variant).strip() + "\t" +
                            "." + "\t" +
                            "." + "\t" +
                            str(row.INFO).strip() + "\n")

    data_combined = pd.concat(data_grouped)


def import_csv_data(in_data):
    header_line = ""
    comment_lines = []

    input_vcf = open(in_data, "r")

    # extract header from vcf file
    for line in input_vcf:
        if line.strip().startswith("##"):
            comment_lines.append(line)
        if line.strip().startswith("#CHROM"):
            header_line = line.strip()
            comment_lines.append(line)
            break # now the variant entries are coming
        else:
            continue

    if header_line == "":
        print("ERROR: The VCF file seems to be corrupted")

    # reset file pointer to begin reading at the beginning
    input_vcf.close()

    data = pd.read_csv(in_data, names=header_line.split("\t"), sep="\t", comment="#", low_memory=False)
    data.fillna(".", inplace=True)

    data = data.rename(columns={"#CHROM": "CHROM"})

    #if "indel_ID" not in data.columns:
    #    data["indel_ID"] = data.index + 1
    #    data["indel_ID"] = data.apply(lambda row: "indel_" + str(row["indel_ID"]), axis=1)

    return data, comment_lines


def convert_csv_to_vcf(in_data, out_data, ref_folder):
    input_data, header = import_csv_data(in_data)
    outfile = open(out_data, "w", newline="")
    #write_header(outfile)
    write_data_information_to_file(input_data, outfile, ref_folder, header)
    outfile.close()


if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--in_data", type=str, dest="in_data", metavar="input.csv", required=True, help="CSV file to convert to VCF\n")
    parser.add_argument("--out_data", type=str, dest="out_data", metavar="output.vcf", required=True, help="output VCF file\n")
    parser.add_argument("--hg19_path", type=str, dest="hg19_path", metavar="/path/to/hg19/Homo_sapiens.GRCh37.dna.chromosome.[ID].fa", required=True, help="Path were the reference hg19 is found.\n")
    args = parser.parse_args()

    hg19_folder = args.hg19_path
    convert_csv_to_vcf(args.in_data, args.out_data, args.hg19_path)
