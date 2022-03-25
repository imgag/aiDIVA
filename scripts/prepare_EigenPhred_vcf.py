import csv
import gzip
import sys


infile = sys.argv[1]
outfile = sys.argv[2]

reader = csv.DictReader(gzip.open(infile, "rt"), dialect="excel-tab")

with open(outfile, "w") as writer:
    writer.write("##fileformat=VCFv4.1\n")
    writer.write("##INFO=<ID=EIGEN_PHRED,Number=1,Type=String,Description=\"Eigen-phred score\">\n")
    writer.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

    for row in reader:
        if str(row["Eigen-phred"]) != "":
            writer.write(f"{row['chr']}\t{row['position']}\t.\t{row['ref']}\t{row['alt']}\t.\t.\tEIGEN_PHRED={str(row['Eigen-phred'])}\n")
        else:
            print("WARNING: Skip missing Eigen-phred score!")
