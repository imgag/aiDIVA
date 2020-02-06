import pandas as pd
import argparse


def write_header(out_file):
    out_file.write("##fileformat=VCFv4.1\n")
    out_file.write("##FILTER=<ID=PASS,Description=\"All filters passed\">\n")
    out_file.write("##FILTER=<ID=LowQual,Description=\"Low quality\">\n")
    out_file.write("##ALT=<ID=NON_REF,Description=\"Represents any possible alternative allele at this location\">\n")
    out_file.write("##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency, for each ALT allele, in the same order as listed\">\n")
    out_file.write("##INFO=<ID=RANK,Number=.,Type=String,Description=\"Rank that indicates whether or not the variant is pathogenic (1 -> pathogenic, 0 -> benign)\">\n")
    out_file.write("##contig=<ID=1,length=249250621,assembly=hg19>\n")
    out_file.write("##contig=<ID=2,length=243199373,assembly=hg19>\n")
    out_file.write("##contig=<ID=3,length=198022430,assembly=hg19>\n")
    out_file.write("##contig=<ID=4,length=191154276,assembly=hg19>\n")
    out_file.write("##contig=<ID=5,length=180915260,assembly=hg19>\n")
    out_file.write("##contig=<ID=6,length=171115067,assembly=hg19>\n")
    out_file.write("##contig=<ID=7,length=159138663,assembly=hg19>\n")
    out_file.write("##contig=<ID=8,length=146364022,assembly=hg19>\n")
    out_file.write("##contig=<ID=9,length=141213431,assembly=hg19>\n")
    out_file.write("##contig=<ID=10,length=135534747,assembly=hg19>\n")
    out_file.write("##contig=<ID=11,length=135006516,assembly=hg19>\n")
    out_file.write("##contig=<ID=12,length=133851895,assembly=hg19>\n")
    out_file.write("##contig=<ID=13,length=115169878,assembly=hg19>\n")
    out_file.write("##contig=<ID=14,length=107349540,assembly=hg19>\n")
    out_file.write("##contig=<ID=15,length=102531392,assembly=hg19>\n")
    out_file.write("##contig=<ID=16,length=90354753,assembly=hg19>\n")
    out_file.write("##contig=<ID=17,length=81195210,assembly=hg19>\n")
    out_file.write("##contig=<ID=18,length=78077248,assembly=hg19>\n")
    out_file.write("##contig=<ID=19,length=59128983,assembly=hg19>\n")
    out_file.write("##contig=<ID=20,length=63025520,assembly=hg19>\n")
    out_file.write("##contig=<ID=21,length=48129895,assembly=hg19>\n")
    out_file.write("##contig=<ID=22,length=51304566,assembly=hg19>\n")
    out_file.write("##contig=<ID=X,length=155270560,assembly=hg19>\n")
    out_file.write("##contig=<ID=Y,length=59373566,assembly=hg19>\n")
    out_file.write("##contig=<ID=MT,length=16571,assembly=hg19>\n")
    out_file.write("#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO\n")


def write_data_information_to_file(outfile, input_data):
    for index, row in input_data.iterrows():
        outfile.write(str(row['Chr']) + '\t' + 
                    str(row['Position']) + '\t' + 
                    str(row['dbsnpIdentifier']) + '\t' + 
                    str(row['Reference']) + '\t' + 
                    str(row['Alteration']) + '\t' + 
                    '.' + '\t' + 
                    '.' + '\t' + 
                    #'AF=' + str(row['AlleleFrequency'])[:6] + ';' + 
                    'RANK=' + str(row['rank'])[:6] + '\n')


def import_csv_data(in_data):
    data = pd.read_csv(in_data, sep='\t')
    data.fillna('.', inplace=True)
    
    return data


def convert_csv_to_vcf(out_data, in_data):
    input_data = import_csv_data(in_data)
    outfile = open(out_data, 'w', newline='')
    write_header(outfile)
    write_data_information_to_file(outfile, input_data)
    outfile.close()


if __name__=='__main__':  
    parser = argparse.ArgumentParser()
    parser.add_argument('--in_data', type=str, dest='in_data', metavar='input.csv', required=True, help='CSV file to convert to VCF\n')
    parser.add_argument('--out_data', type=str, dest='out_data', metavar='output.vcf', required=True, help='output VCF file\n')
    args = parser.parse_args()
    
    convert_csv_to_vcf(args.out_data, args.in_data)
