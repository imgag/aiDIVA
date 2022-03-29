import gzip
from re import S
from statistics import quantiles
import sys


infile = sys.argv[1]
outfile = sys.argv[2]

with gzip.open(infile, "rt") as gnomad, open(outfile, "w") as prepared_gnomad:

    for line in gnomad:
        if line.startswith("#"):
            if line.startswith("##fileformat="):
                prepared_gnomad.write(line)
            elif line.startswith("##hailversion="):
                prepared_gnomad.write(line)
            elif line.startswith("##INFO=<ID=AN,"):
                prepared_gnomad.write(line)
            elif line.startswith("##INFO=AF,"):
                prepared_gnomad.write(line)
            elif line.startswith("##INFO=<ID=nhomalt,"):
                prepared_gnomad.write("##INFO=<ID=Hom,Number=A,Type=Integer,Description=\"Count of homozygous individuals in samples\">\n")
            elif line.startswith("##contig="):
                prepared_gnomad.write(line)
            elif line.startswith("#CHROM"):
                prepared_gnomad.write(line)
            else:
                continue
        else:
            splitted_line = line.replace("\n", "").split("\t")
            info_fields = splitted_line[7].split(";")

            an_entry = ""
            af_entry = ""
            hom_entry = ""

            for entry in info_fields:
                if entry.startswith("AN="):
                    an_entry = entry.split("=")[1]
                elif entry.startswith("AF="):
                    af_entry = entry.split("=")[1]
                elif entry.startswith("nhomalt="):
                    hom_entry = entry.split("=")[1]
                else:
                    continue
            
            chrom = splitted_line[0]
            pos = splitted_line[1]
            id = "."
            ref = splitted_line[3]
            alt = splitted_line[4]
            qual = "."
            filt = "."
            info = f"AN={an_entry};AF={af_entry};Hom={hom_entry}"

            if str(an_entry) != "" and str(af_entry) != "" and str(hom_entry) != "":
                prepared_gnomad.write(f"{chrom}\t{pos}\t{id}\t{ref}\t{alt}\t{qual}\t{filt}\t{info}\n")
            else:
                print("WARNING: Skip missing gnomAD entry!")
