import argparse
import logging
import pandas as pd
import pysam
import random


logger = logging.getLogger(__name__)


def write_data_information_to_file(input_data, outfile, ref_sequence, header):
    data_grouped = [group for key, group in input_data.groupby("CHROM")]
    in_fasta = pysam.FastaFile(ref_sequence)
    random.seed(14038)

    for line in header:
        if line.startswith("#CHROM"):
            outfile.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        else:
            outfile.write(line)

    for group in data_grouped:
        if "chr" in str(group["CHROM"].iloc[0]):
            chrom_id = str(group["CHROM"].iloc[0])
        else:
            chrom_id = "chr" + str(group["CHROM"].iloc[0])

        for row in group.itertuples():
            # make sure that the window starts at the beginning of the reference sequence
            window_start = max(int(row.POS) - 3, 1)

            # make sure that the window won't exceed the reference sequence length
            window_end = min(int(row.POS) + len(row.REF) + 2, in_fasta.get_reference_length(chrom_id))
            extended_ref_seq = in_fasta.fetch(chrom_id, window_start, window_end)

            for i in range(abs(window_end-window_start)):
                alt_variant = ""
                if (extended_ref_seq[i] == "A") or (extended_ref_seq[i] == "T"):
                    alt_variant = random.choice(["G", "C"])
                elif (extended_ref_seq[i] == "G") or (extended_ref_seq[i] == "C"):
                    alt_variant = random.choice(["A", "T"])
                elif (extended_ref_seq[i] == "N"):
                    logger.debug("Reference base was skipped because it was 'N'!")
                    continue
                else:
                    logger.error("The given reference sequence seems to be corrupted!")

                outfile.write(str(row.CHROM).strip() + "\t" + str(window_start + i + 1).strip() + "\t" + "." + "\t" + str(extended_ref_seq[i]).strip() + "\t" + str(alt_variant).strip() + "\t" + "." + "\t" + "." + "\t" + str(row.INFO).strip() + "\n")


def import_vcf_data(in_data):
    header_line = ""
    comment_lines = []

    with open(in_data, "r") as input_vcf:
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
            logger.error("The VCF seems to be corrupted, missing header line!")
        
        # reset file pointer to begin reading at the beginning (is done by closing the file before reading again)

    data = pd.read_csv(in_data, names=header_line.split("\t"), sep="\t", comment="#", low_memory=False)
    data.fillna(".", inplace=True)
    data = data.rename(columns={"#CHROM": "CHROM"})

    return data, comment_lines


def convert_indel_vcf_to_expanded_indel_vcf(in_data, out_data, ref_folder):
    input_data, header = import_vcf_data(in_data)
    with open(out_data, "w", newline="") as outfile:
        write_data_information_to_file(input_data, outfile, ref_folder, header)
    

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--in_data", type=str, dest="in_data", metavar="input.csv", required=True, help="InDel VCF file to expand\n")
    parser.add_argument("--out_data", type=str, dest="out_data", metavar="output.vcf", required=True, help="Output VCF file\n")
    parser.add_argument("--ref_path", type=str, dest="ref_path", metavar="/path/to/hg19/Homo_sapiens.GRCh37.dna.chromosome.[ID].fa", required=True, help="Path were the reference genome is found.\n")
    args = parser.parse_args()

    ref_folder = args.ref_path
    convert_indel_vcf_to_expanded_indel_vcf(args.in_data, args.out_data, args.ref_path)
