import pandas as pd
import numpy as np
import argparse
import random
from Bio import SeqIO


def write_data_information_to_file(input_data, outfile, ref_sequence, header):
    data_grouped = [group for key, group in input_data.groupby("CHROM")]

    random.seed(14038)

    for line in header:
        if line.startswith("#CHROM"):
            outfile.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        else:
            outfile.write(line)
    ref_seq_records = SeqIO.index(ref_sequence, "fasta")

    for group in data_grouped:
        if "chr" in str(group["CHROM"].iloc[0]):
            chrom_id = str(group["CHROM"].iloc[0])
        else:
            chrom_id = "chr" + str(group["CHROM"].iloc[0])

        ref_seq = str(ref_seq_records[chrom_id].seq)
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

    #data_combined = pd.concat(data_grouped)


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
