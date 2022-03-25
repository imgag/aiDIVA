import csv
import gzip
import sys
import zipfile


infile = sys.argv[1]
outfile_grch37 = sys.argv[2]
outfile_grch38 = sys.argv[3]

reader = csv.DictReader(open(infile, "r"), dialect="excel")

with open(outfile_grch37, "w") as grch37_writer, open(outfile_grch38, "w") as grch38_writer:
    grch37_writer.write("##fileformat=VCFv4.1\n")
    grch37_writer.write("##INFO=<ID=REVEL,Number=1,Type=String,Description=\"REVEL score\">\n")
    grch37_writer.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

    grch38_writer.write("##fileformat=VCFv4.1\n")
    grch38_writer.write("##INFO=<ID=REVEL,Number=1,Type=String,Description=\"REVEL score\">\n")
    grch38_writer.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

    for row in reader:
        if str(row["REVEL"]) != "":
            grch37_writer.write(f"{row['chr']}\t{row['hg19_pos']}\t.\t{row['ref']}\t{row['alt']}\t.\t.\tREVEL={str(row['REVEL'])}\n")
            if row["grch38_pos"] != ".":
                grch38_writer.write(f"{row['chr']}\t{row['grch38_pos']}\t.\t{row['ref']}\t{row['alt']}\t.\t.\tREVEL={str(row['REVEL'])}\n")
            else:
                print("WARNING: Skip GRCh38 missing REVEL score!")
        else:
            print("WARNING: Skip missing REVEL score!")
