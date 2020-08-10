import pandas as pd
import numpy as np
import argparse


def write_header(out_file):
    out_file.write("##fileformat=VCFv4.1\n")
    out_file.write("##INFO=<ID=AIDIVA_SCORE,Number=.,Type=String,Description=\"Pathogenicity score of variant predicted with a random forest model.\">\n")
    out_file.write("##INFO=<ID=AIDIVA_SCORE_FINAL,Number=.,Type=String,Description=\"AIDIVA_SCORE score adjusted with the HPO relatedness.\">\n")
    out_file.write("##INFO=<ID=AIDIVA_HPO_RELATEDNESS,Number=.,Type=String,Description=\"Score that shows the HPO relatedness.\">\n")
    out_file.write("##INFO=<ID=AIDIVA_DOMINANT,Number=.,Type=String,Description=\"Flag indicating dominant inheritance.\">\n")
    out_file.write("##INFO=<ID=AIDIVA_DENOVO,Number=.,Type=String,Description=\"Flag indicating a denovo variant.\">\n")
    out_file.write("##INFO=<ID=AIDIVA_RECESSIVE,Number=.,Type=String,Description=\"Flag indicating recessive inheritance.\">\n")
    out_file.write("##INFO=<ID=AIDIVA_XLINKED,Number=.,Type=String,Description=\"Flag indicating xlinked inheritance.\">\n")
    out_file.write("##INFO=<ID=AIDIVA_COMPOUND,Number=.,Type=String,Description=\"Flag indicating a possible candidate for compound inheritance.\">\n")
    out_file.write("##INFO=<ID=AIDIVA_FILTER,Number=.,Type=String,Description=\"Flag indicating if all filters are passed.\">\n")
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
    out_file.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

def write_result_vcf(input_data, vcf_file):
    input_data.sort_values(["CHROM", "POS"], ascending=[True, True], inplace=True)
    input_data.reset_index(inplace=True, drop=True)

    with open(vcf_file, "w") as out:
        write_header(out)

        inheritance_information = ("DOMINANT_DENOVO" in input_data.columns) & ("DOMINANT_INHERITED" in input_data.columns) & ("RECESSIVE" in input_data.columns) & ("XLINKED" in input_data.columns) & ("COMPOUND" in input_data.columns)
        for row in input_data.itertuples():
            if inheritance_information:
                info_entry = "AIDIVA_SCORE=" + str(row.AIDIVA_SCORE) + ";AIDIVA_SCORE_FINAL=" + str(row.FINAL_AIDIVA_SCORE) + ";AIDIVA_DOMINANT=" + str(row.DOMINANT_INHERITED) + ";AIDIVA_DENOVO=" + str(row.DOMINANT_DENOVO) + ";AIDIVA_RECESSIVE=" + str(row.RECESSIVE) + ";AIDIVA_XLINKED=" + str(row.XLINKED) + ";AIDIVA_COMPOUND=" + str(row.COMPOUND) + ";AIDIVA_HPO_RELATEDNESS=" + str(row.HPO_RELATEDNESS) + ";AIDIVA_FILTER=" + str(row.FILTER_PASSED)
            else:
                info_entry = "AIDIVA_SCORE=" + str(row.AIDIVA_SCORE) + ";AIDIVA_SCORE_FINAL=" + str(row.FINAL_AIDIVA_SCORE) + ";AIDIVA_HPO_RELATEDNESS=" + str(row.HPO_RELATEDNESS) + ";AIDIVA_FILTER=" + str(row.FILTER_PASSED)
            out.write(str(row.CHROM).strip() + "\t" + str(row.POS) + "\t" + "." + "\t" + str(row.REF) + "\t" + str(row.ALT) + "\t" + "." + "\t" + "." + "\t" + info_entry + "\n")


if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--in", type=str, dest="in_data", metavar="in.csv", required=True, help="CSV file with the results from the AIdiva run that should be converted to a VCF file\n")
    parser.add_argument("--out", type=str, dest="out_data", metavar="out.vcf", required=True, help="VCF file witht the results of the AIdiva run\n")
    args = parser.parse_args()

    in_data = pd.read_csv(args.in_data, sep="\t", low_memory=False)
    write_result_vcf(in_data, args.out_data)
