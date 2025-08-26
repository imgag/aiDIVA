####################################################################
### Small script to prepare the CADD data for easier annotation! ###
####################################################################

import gzip
import sys


infile = sys.argv[1]
outfile = sys.argv[2]

with gzip.open(infile, "rt") as infile, open(outfile, "w") as writer:
    writer.write("##fileformat=VCFv4.1\n")
    writer.write("##INFO=<ID=CADD_PHRED,Number=1,Type=String,Description=\"CADD phred score\">\n")
    writer.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

    for line in infile:
        if line.startswith("#"):
            continue

        split = line.replace("\n", "").split("\t")

        if str(split[4]) != "":
            writer.write(f"{split[0]}\t{split[1]}\t.\t{split[2]}\t{split[3]}\t.\t.\tCADD_PHRED={str(split[5])}\n")

        else:
            print("WARNING: Skip missing CADD phred score!")
