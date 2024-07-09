## TODO: Delete file it is not needed anymore
import csv
import sys


chroms = str(sys.argv[1]).split(",")
grch37_outfile = sys.argv[2]
grch38_outfile = sys.argv[3]

with open(grch37_outfile, "w") as writer_grch37, open(grch38_outfile, "w") as writer_grch38:
    writer_grch37.write("##fileformat=VCFv4.1\n")
    writer_grch37.write("##INFO=<ID=ADA_SCORE,Number=1,Type=String,Description=\"Ada score from dbscSNV (hg19 coordinates)\">\n")
    writer_grch37.write("##INFO=<ID=RF_SCORE,Number=1,Type=String,Description=\"Rf score from dbscSNV (hg19 coordinates)\">\n")
    writer_grch37.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

    writer_grch38.write("##fileformat=VCFv4.1\n")
    writer_grch38.write("##INFO=<ID=ADA_SCORE,Number=1,Type=String,Description=\"Ada score from dbscSNV (hg38 coordinates)\">\n")
    writer_grch38.write("##INFO=<ID=RF_SCORE,Number=1,Type=String,Description=\"Rf score from dbscSNV (hg38 coordinates)\">\n")
    writer_grch38.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

    for chrom in chroms:
        reader = csv.DictReader(open(chrom, "rt"), dialect="excel-tab")
        for row in reader:

            # skip missing scores to keep the annotation file small
            if str(row["ada_score"]) != "" and str(row["rf_score"]) != "":
                writer_grch37.write(f"{row['chr']}\t{row['pos']}\t.\t{row['ref']}\t{row['alt']}\t.\t.\tADA_SCORE={str(row['ada_score'])};RF_SCORE={str(row['rf_score'])}\n")
                if row["hg38_pos"] != ".":
                    writer_grch38.write(f"{row['hg38_chr']}\t{row['hg38_pos']}\t.\t{row['ref']}\t{row['alt']}\t.\t.\tADA_SCORE={str(row['ada_score'])};RF_SCORE={str(row['rf_score'])}\n")

                else:
                    print("WARNING: Skip GRCh38 missing dbscSNV score!")

            else:
                print("WARNING: Skip missing scores!")