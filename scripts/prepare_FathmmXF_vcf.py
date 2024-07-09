import csv
import gzip
import sys

infile = sys.argv[1]
outfile = sys.argv[2]

reader = csv.DictReader(gzip.open(infile, "rt"), dialect="excel-tab", fieldnames=["#CHROM", "POS", "REF", "ALT", "FATHMM_XF"])

with open(outfile, "w") as writer:
    writer.write("##fileformat=VCFv4.1\n")
    writer.write("##INFO=<ID=FATHMM_XF,Number=1,Type=String,Description=\"FATHMM-XF score\">\n")
    writer.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

    for row in reader:
        # skip missing scores to keep the annotation file small
        if str(row["FATHMM_XF"]) != "":
            writer.write(f"{row['#CHROM']}\t{row['POS']}\t.\t{row['REF']}\t{row['ALT']}\t.\t.\FATHMM_XF={str(row['FATHMM_XF'])}\n")

        else:
            print("WARNING: Skip missing FathmmXF score!")