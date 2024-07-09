import pandas as pd
import numpy as np
import argparse
import logging


logger = logging.getLogger(__name__)


def write_header(out_file, assembly_build, single):
    out_file.write("##fileformat=VCFv4.2\n")
    if not single:
        out_file.write("##INFO=<ID=AIDIVA,Number=5,Type=String,Description=\"aiDIVA scores: aiDIVA-score,aiDIVA-final-score,aiDIVA-hpo-relatedness,aiDIVA-hpo-relatedness-interacting,aiDIVA-filter. (aiDIVA-score is the pathogenicity prediction from the random forest; aiDIVA-final-score is the predction finalized with the given HPO terms; aiDIVA-hpo-relatedness indicates how strong the currrent variant is associated with the given HPO terms; aiDIVA-filter 0 or 1 wether all internal filters were passed or not)\">\n")
        out_file.write("##INFO=<ID=AIDIVA_INHERITANCE,Number=4,Type=String,Description=\"aiDIVA inheritance flags: dominant,denovo,recessive,xlinked,compound. (Each value can be 0 or 1)\">\n")
        out_file.write("##INFO=<ID=AIDIVA_INHERITANCE_COMMENT,Number=1,Type=String,Description=\"aiDIVA inheritance flags (dominant,denovo,recessive,xlinked,compound) in written form\">\n")

    else:
        out_file.write("##INFO=<ID=AIDIVA,Number=5,Type=String,Description=\"aiDIVA scores: aiDIVA-score,aiDIVA-final-score,aiDIVA-hpo-relatedness,aiDIVA-hpo-relatedness-interacting,aiDIVA-filter. (aiDIVA-score is the pathogenicity prediction from the random forest; aiDIVA-final-score is the predction finalized with the given HPO terms; aiDIVA-hpo-relatedness indicates how strong the currrent variant is associated with the given HPO terms; aiDIVA-filter 0 or 1 wether all internal filters were passed or not)\">\n")
        out_file.write("##INFO=<ID=AIDIVA_INHERITANCE,Number=2,Type=String,Description=\"aiDIVA inheritance flags: recessive,compound. (Each value can be 0 or 1)\">\n")
        out_file.write("##INFO=<ID=AIDIVA_INHERITANCE_COMMENT,Number=1,Type=String,Description=\"aiDIVA inheritance flags (recessive,compound) in written form\">\n")

    if assembly_build == "GRCh38":
        out_file.write("##reference=GRCh38.fa\n")
        out_file.write("##contig=<ID=chr1,length=248956422>\n")
        out_file.write("##contig=<ID=chr2,length=242193529>\n")
        out_file.write("##contig=<ID=chr3,length=198295559>\n")
        out_file.write("##contig=<ID=chr4,length=190214555>\n")
        out_file.write("##contig=<ID=chr5,length=181538259>\n")
        out_file.write("##contig=<ID=chr6,length=170805979>\n")
        out_file.write("##contig=<ID=chr7,length=159345973>\n")
        out_file.write("##contig=<ID=chr8,length=145138636>\n")
        out_file.write("##contig=<ID=chr9,length=138394717>\n")
        out_file.write("##contig=<ID=chr10,length=133797422>\n")
        out_file.write("##contig=<ID=chr11,length=135086622>\n")
        out_file.write("##contig=<ID=chr12,length=133275309>\n")
        out_file.write("##contig=<ID=chr13,length=114364328>\n")
        out_file.write("##contig=<ID=chr14,length=107043718>\n")
        out_file.write("##contig=<ID=chr15,length=101991189>\n")
        out_file.write("##contig=<ID=chr16,length=90338345>\n")
        out_file.write("##contig=<ID=chr17,length=83257441>\n")
        out_file.write("##contig=<ID=chr18,length=80373285>\n")
        out_file.write("##contig=<ID=chr19,length=58617616>\n")
        out_file.write("##contig=<ID=chr20,length=64444167>\n")
        out_file.write("##contig=<ID=chr21,length=46709983>\n")
        out_file.write("##contig=<ID=chr22,length=50818468>\n")
        out_file.write("##contig=<ID=chrX,length=156040895>\n")
        out_file.write("##contig=<ID=chrY,length=57227415>\n")

    elif assembly_build == "GRCh37":
        out_file.write("##reference=GRCh37.fa\n")
        out_file.write("##contig=<ID=chr1,length=249250621>\n")
        out_file.write("##contig=<ID=chr2,length=243199373>\n")
        out_file.write("##contig=<ID=chr3,length=198022430>\n")
        out_file.write("##contig=<ID=chr4,length=191154276>\n")
        out_file.write("##contig=<ID=chr5,length=180915260>\n")
        out_file.write("##contig=<ID=chr6,length=171115067>\n")
        out_file.write("##contig=<ID=chr7,length=159138663>\n")
        out_file.write("##contig=<ID=chr8,length=146364022>\n")
        out_file.write("##contig=<ID=chr9,length=141213431>\n")
        out_file.write("##contig=<ID=chr10,length=135534747>\n")
        out_file.write("##contig=<ID=chr11,length=135006516>\n")
        out_file.write("##contig=<ID=chr12,length=133851895>\n")
        out_file.write("##contig=<ID=chr13,length=115169878>\n")
        out_file.write("##contig=<ID=chr14,length=107349540>\n")
        out_file.write("##contig=<ID=chr15,length=102531392>\n")
        out_file.write("##contig=<ID=chr16,length=90354753>\n")
        out_file.write("##contig=<ID=chr17,length=81195210>\n")
        out_file.write("##contig=<ID=chr18,length=78077248>\n")
        out_file.write("##contig=<ID=chr19,length=59128983>\n")
        out_file.write("##contig=<ID=chr20,length=63025520>\n")
        out_file.write("##contig=<ID=chr21,length=48129895>\n")
        out_file.write("##contig=<ID=chr22,length=51304566>\n")
        out_file.write("##contig=<ID=chrX,length=155270560>\n")
        out_file.write("##contig=<ID=chrY,length=59373566>\n")

    else:
        logger.error(f"Unsupported assembly build: {assembly_build}")

    out_file.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

def write_result_vcf(input_data, vcf_file, assembly_build, single):
    with open(vcf_file, "w") as out_file:
        if input_data is not None:
            input_data = input_data.rename(columns={"#CHROM": "CHROM"})
            input_data = input_data.sort_values(["CHROM", "POS"], ascending=[True, True])
            input_data = input_data.reset_index(drop=True)
            colnames = input_data.columns

            write_header(out_file, assembly_build, single)

            for row in input_data.itertuples():
                if ("AIDIVA_SCORE" in colnames):
                    if (str(row.AIDIVA_SCORE) == "nan"):
                        aidiva_score = "."

                    else:
                        aidiva_score = str(row.AIDIVA_SCORE)

                else:
                    aidiva_score = "."

                if ("FINAL_AIDIVA_SCORE" in colnames):
                    if (str(row.FINAL_AIDIVA_SCORE) == "nan"):
                        final_aidiva_score = "."

                    else:
                        final_aidiva_score = str(row.FINAL_AIDIVA_SCORE)

                else:
                    final_aidiva_score = "."

                if ("HPO_RELATEDNESS" in colnames):
                    if (str(row.HPO_RELATEDNESS) == "nan"):
                        hpo_relatedness = "."

                    else:
                        hpo_relatedness = str(row.HPO_RELATEDNESS)

                else:
                    hpo_relatedness = "."

                if ("HPO_RELATEDNESS_INTERACTING" in colnames):
                    if (str(row.HPO_RELATEDNESS_INTERACTING) == "nan"):
                        hpo_relatedness_interacting = "."

                    else:
                        hpo_relatedness_interacting = str(row.HPO_RELATEDNESS_INTERACTING)

                else:
                    hpo_relatedness_interacting = "."

                if ("FILTER_PASSED" in colnames):
                    if (str(row.FILTER_PASSED) == "nan"):
                        filter_passed = "."

                    else:
                        filter_passed = str(row.FILTER_PASSED)

                else:
                    filter_passed = "."

                if ("DOMINANT" in colnames):
                    if (str(row.DOMINANT) == "nan"):
                        dominant = "."

                    else:
                        dominant = str(row.DOMINANT)

                else:
                    dominant = "."

                if ("DOMINANT_DENOVO" in colnames):
                    if (str(row.DOMINANT_DENOVO) == "nan"):
                        dominant_denovo = "."

                    else:
                        dominant_denovo = str(row.DOMINANT_DENOVO)

                else:
                    dominant_denovo = "."

                if ("RECESSIVE" in colnames):
                    if (str(row.RECESSIVE) == "nan"):
                        recessive = "."

                    else:
                        recessive = str(row.RECESSIVE)

                else:
                    recessive = "."

                if ("XLINKED" in colnames):
                    if (str(row.XLINKED) == "nan"):
                        xlinked = "."

                    else:
                        xlinked = str(row.XLINKED)

                else:
                    xlinked = "."

                if ("COMPOUND" in colnames):
                    if (str(row.COMPOUND) == "nan"):
                        compound = "."

                    else:
                        compound = str(row.COMPOUND)

                else:
                    compound = "."

                if ("INHERITANCE" in colnames):
                    if (str(row.INHERITANCE) == "nan") or (str(row.INHERITANCE) == ""):
                        inheritance_comment = "."

                    else:
                        inheritance_comment = str(row.INHERITANCE)

                else:
                    inheritance_comment = "."

                if not single:
                    info_entry = "AIDIVA=" + aidiva_score + "," + final_aidiva_score + "," + hpo_relatedness + "," + hpo_relatedness_interacting + "," + filter_passed + ";AIDIVA_INHERITANCE=" + dominant + "," + dominant_denovo + "," + recessive + "," + xlinked + "," + compound + ";AIDIVA_INHERITANCE_COMMENT=" + inheritance_comment

                else:
                    info_entry = "AIDIVA=" + aidiva_score + "," + final_aidiva_score + "," + hpo_relatedness + "," + hpo_relatedness_interacting + "," + filter_passed + ";AIDIVA_INHERITANCE=" + recessive + "," + compound + ";AIDIVA_INHERITANCE_COMMENT=" + inheritance_comment

                out_file.write(str(row.CHROM).strip() + "\t" + str(row.POS) + "\t" + "." + "\t" + str(row.REF) + "\t" + str(row.ALT) + "\t" + "." + "\t" + "." + "\t" + info_entry + "\n")

        else:
            write_header(out_file, assembly_build, single)
            logger.warning(f"Wrote an empty result file!")

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--in", type=str, dest="in_data", metavar="in.csv", required=True, help="CSV file with the results from the aiDIVA run that should be converted to a VCF file\n")
    parser.add_argument("--out", type=str, dest="out_data", metavar="out.vcf", required=True, help="VCF file witht the results of the aiDIVA run\n")
    parser.add_argument("--not_single", action="store_false", required=False, help="Flag to indicate whether the file is single sample or not.\n")
    args = parser.parse_args()

    in_data = pd.read_csv(args.in_data, sep="\t", low_memory=False)
    write_result_vcf(in_data, args.out_data, args.not_single)
