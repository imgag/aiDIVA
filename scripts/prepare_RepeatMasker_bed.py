import gzip
import re
import sys


infile = sys.argv[1]
outfile = sys.argv[2]

with gzip.open(infile, "rt") as infile, open(outfile, "w") as writer:
    
    for line in infile:
        current_line = line.lstrip().replace("\n", "")
        current_line = re.sub(r"\s+", " ", current_line)

        if current_line == "":
            continue

        if current_line.startswith("SW"):
            continue

        if current_line.startswith("score"):
            continue

        split = current_line.split(" ")
        
        if str(split[9]) != "" and str(split[10]) != "":
            writer.write(f"{split[4]}\t{str(int(split[5])-1)}\t{split[6]}\t{split[9]} ({split[10]})\n")
        else:
            print("WARNING: Skip missing repeat!")