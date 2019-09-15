import sys
import pandas


data = pandas.read_csv(sys.argv[1], sep="\t")
data.fillna('.', inplace=True)

header = """##fileformat=VCFv4.1
##contig=<ID=1,length=249250621,assembly=hg19>
##contig=<ID=2,length=243199373,assembly=hg19>
##contig=<ID=3,length=198022430,assembly=hg19>
##contig=<ID=4,length=191154276,assembly=hg19>
##contig=<ID=5,length=180915260,assembly=hg19>
##contig=<ID=6,length=171115067,assembly=hg19>
##contig=<ID=7,length=159138663,assembly=hg19>
##contig=<ID=8,length=146364022,assembly=hg19>
##contig=<ID=9,length=141213431,assembly=hg19>
##contig=<ID=10,length=135534747,assembly=hg19>
##contig=<ID=11,length=135006516,assembly=hg19>
##contig=<ID=12,length=133851895,assembly=hg19>
##contig=<ID=13,length=115169878,assembly=hg19>
##contig=<ID=14,length=107349540,assembly=hg19>
##contig=<ID=15,length=102531392,assembly=hg19>
##contig=<ID=16,length=90354753,assembly=hg19>
##contig=<ID=17,length=81195210,assembly=hg19>
##contig=<ID=18,length=78077248,assembly=hg19>
##contig=<ID=19,length=59128983,assembly=hg19>
##contig=<ID=20,length=63025520,assembly=hg19>
##contig=<ID=21,length=48129895,assembly=hg19>
##contig=<ID=22,length=51304566,assembly=hg19>
##contig=<ID=X,length=155270560,assembly=hg19>
##contig=<ID=Y,length=59373566,assembly=hg19>
##contig=<ID=MT,length=16571,assembly=hg19>
##INFO=<ID=RANK,Number=.,Type=String,Description="Rank that indicates whether or not the variant is pathogenic">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO\n"""

with open(sys.argv[2], "w", newline="") as outfile:
    outfile.write(header)
    count = 0
    for index, row in data.iterrows():
        print(count)
        outfile.write(str(row['Chr']) + "\t" + str(row['Position']) + "\t" + str(row['dbsnpIdentifier']) + "\t" + str(row['Reference']) + "\t" + str(row['Alteration']) + "\t" + "." + "\t" + "." + "\t" + "RANK=" + str(row['rank']) + "\n")
        count += 1

print('Finished')