import csv
import sys

chroms = str(sys.argv[1]).split(",")
outfile = sys.argv[2]

with open(outfile, "w") as writer:
    writer.write("##fileformat=VCFv4.1\n")
    writer.write("##INFO=<ID=FATHMM_XF,Number=1,Type=String,Description=\"FATHMM-XF score\">\n")
    writer.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

    for chrom in chroms:
        reader = csv.DictReader(open(chrom, "rt"), dialect="excel")
        for row in reader:
            mutation = str(row["Mutation"]).split(",")

            # skip missing scores to keep the annotation file small
            if str(row["FI score"]) != "":
                writer.write(f"{mutation[1]}\t{mutation[2]}\t.\t{mutation[3]}\t{mutation[4]}\t.\t.\tMutationAssessor={str(row['FI score'])}\n")
            else:
                print("WARNING: Skip missing MutationAssessor score!")