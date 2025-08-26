######################################################################
### Small script to prepare the Condel data for easier annotation! ###
######################################################################

import csv
import gzip
import sys


infile = sys.argv[1]
outfile = sys.argv[2]

reader = csv.DictReader(gzip.open(infile, "rt"), dialect="excel-tab")

with open(outfile, "w") as writer:
    writer.write("##fileformat=VCFv4.1\n")
    writer.write("##INFO=<ID=CONDEL,Number=1,Type=String,Description=\"Condel score from fannsDB\">\n")
    writer.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

    for row in reader:
        if str(row["CONDEL"]) != "":
            writer.write(f"{row['CHR']}\t{row['START']}\t.\t{row['REF']}\t{row['ALT']}\t.\t.\tEIGEN_PHRED={str(row['CONDEL'])}\n")

        else:
            print("WARNING: Skip missing Condel score!")
