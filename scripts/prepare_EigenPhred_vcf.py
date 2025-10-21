#####################################################################
### Small script to prepare the Eigen data for easier annotation! ###
#####################################################################

import csv
import gzip
import sys


infile = sys.argv[1]
outfile = sys.argv[2]

file_header = ['chr', 'pos', "ref", "alt", "unknown1", "unknown2", "unknown3", "unknown4", "unknown5", "unknown6", "unknown7", "unknown8", "unknown9", "unknown10", "unknown11", "unknown12", 
"unknown13", "unknown14", "Eigen-phred", "unknown16", "unknown17"]
reader = csv.DictReader(gzip.open(infile, "rt"), dialect="excel-tab", fieldnames=file_header)

with open(outfile, "w") as writer:
    writer.write("##fileformat=VCFv4.1\n")
    writer.write("##INFO=<ID=EIGEN_PHRED,Number=1,Type=String,Description=\"Eigen-phred score\">\n")
    writer.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

    for row in reader:
        # The header line in the table from Zenodo is not the first line, but instead the last line. To prevent a faulty VCF entry we skip this line.
        if str(row['chr']) == "chr":
            continue

        if str(row["Eigen-phred"]) != "":
            writer.write(f"{row['chr']}\t{row['pos']}\t.\t{row['ref']}\t{row['alt']}\t.\t.\tEIGEN_PHRED={str(row['Eigen-phred'])}\n")

        else:
            print("WARNING: Skip missing Eigen-phred score!")
