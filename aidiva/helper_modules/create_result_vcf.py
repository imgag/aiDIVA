import pandas as pd
import numpy as np
import argparse


def write_header(out_file, single):
    out_file.write("##fileformat=VCFv4.1\n")
    if not single:
        out_file.write("##INFO=<ID=AIDIVA,Number=4,Type=String,Description=\"AIdiva scores: AIdiva-score,AIdiva-final-score,AIdiva-hpo-relatedness,AIdiva-filter. (AIdiva-score is the pathogenicity prediction from the random forest; AIdiva-final-score is the predction finalized with the given HPO terms; AIdiva-hpo-relatedness indicates how strong the currrent variant is associated with the given HPO terms; AIdiva-filter 0 or 1 wether all internal filters were passed or not)\">\n")
        out_file.write("##INFO=<ID=AIDIVA_INHERITANCE,Number=4,Type=String,Description=\"AIdiva inheritance flags: dominant,denovo,recessive,xlinked,compound. (Each value can be 0 or 1)\">\n")
    else:
        out_file.write("##INFO=<ID=AIDIVA,Number=4,Type=String,Description=\"AIdiva scores: AIdiva-score,AIdiva-final-score,AIdiva-hpo-relatedness,AIdiva-filter. (AIdiva-score is the pathogenicity prediction from the random forest; AIdiva-final-score is the predction finalized with the given HPO terms; AIdiva-hpo-relatedness indicates how strong the currrent variant is associated with the given HPO terms; AIdiva-filter 0 or 1 wether all internal filters were passed or not)\">\n")

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

def write_result_vcf(input_data, vcf_file, single):
    input_data.sort_values(["CHROM", "POS"], ascending=[True, True], inplace=True)
    input_data.reset_index(inplace=True, drop=True)

    with open(vcf_file, "w") as out:
        write_header(out, single)

        for row in input_data.itertuples():
            if not single:
                info_entry = "AIDIVA=" + str(row.AIDIVA_SCORE) + "," + str(row.FINAL_AIDIVA_SCORE) + "," + str(row.HPO_RELATEDNESS) + "," + str(row.FILTER_PASSED) + ";AIDIVA_INHERITANCE=" + str(row.DOMINANT_INHERITED) + "," + str(row.DOMINANT_DENOVO) + "," + str(row.RECESSIVE) + "," + str(row.XLINKED) + "," + str(row.COMPOUND)
            else:
                info_entry = "AIDIVA=" + str(row.AIDIVA_SCORE) + "," + str(row.FINAL_AIDIVA_SCORE) + "," + str(row.HPO_RELATEDNESS) + "," + str(row.FILTER_PASSED)

            out.write(str(row.CHROM).strip() + "\t" + str(row.POS) + "\t" + "." + "\t" + str(row.REF) + "\t" + str(row.ALT) + "\t" + "." + "\t" + "." + "\t" + info_entry + "\n")


if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--in", type=str, dest="in_data", metavar="in.csv", required=True, help="CSV file with the results from the AIdiva run that should be converted to a VCF file\n")
    parser.add_argument("--out", type=str, dest="out_data", metavar="out.vcf", required=True, help="VCF file witht the results of the AIdiva run\n")
    parser.add_argument("--not_single", action="store_false", required=False, help="Flag to indicate whether the file is single sample or not.\n")
    args = parser.parse_args()

    in_data = pd.read_csv(args.in_data, sep="\t", low_memory=False)
    write_result_vcf(in_data, args.out_data, args.not_single)
