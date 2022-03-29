import gzip
import sys


infile = sys.argv[1]
outfile = sys.argv[2]

with gzip.open(infile, "rt") as infile, open(outfile, "w") as writer:
    writer.write("##fileformat=VCFv4.1\n")
    writer.write("##INFO=<ID=DETAILS,Number=.,Type=String,Description=\"ClinVar disease/significance annotation\">\n")
    writer.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    
    for line in infile:
        if line.startswith("#"):
            continue

        split = line.replace("\n", "").split("\t")

        if split[4] == ".":
            continue

        info_fields = split[7].split(";")
        
        disease_name = ""
        disease_name_incl = ""
        sig_info = ""
        sig_info_conf = ""
        sig_info_incl = ""

        for field in info_fields:
            if field.startswith("CLNDN="):
                disease_name = field.split("=")[1]
            elif field.startswith("CLNDNINCL="):
                disease_name_incl = field.split("=")[1]
            elif field.startswith("CLNSIG="):
                sig_info = field.split("=")[1]
            elif field.startswith("CLNSIGCONF="):
                sig_info_conf = field.split("=")[1]
            elif field.startswith("CLNSIGINCL="):
                sig_info_incl = field.split("=")[1]
            else:
                continue
        
        # TODO: resolve conflicting interpretations of pathogenicity if possible and include them in the resulting file

        if str(sig_info) != "" and sig_info != "Conflicting_interpretations_of_pathogenicity":
            writer.write(f"{split[0]}\t{split[1]}\t{split[2]}\t{split[3]}\t{split[4]}\t{split[5]}\t{split[6]}\tDETAILS={str(sig_info)}%20{str(disease_name)}%3D{str(disease_name_incl)}\n")
        else:
            print("WARNING: Skip conflicting or empty significance entry!")