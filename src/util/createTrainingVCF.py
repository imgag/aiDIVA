import sys
import pandas


data = pandas.read_csv(sys.argv[1], sep=',')
data.fillna('.', inplace=True)

header = """##fileformat=VCFv4.1
##FILTER=<ID=PASS,Description="All filters passed">
##FILTER=<ID=LowQual,Description="Low quality">
##ALT=<ID=NON_REF,Description="Represents any possible alternative allele at this location">
##FORMAT=<ID=AB,Number=1,Type=Float,Description="Allele balance for each het genotype">
##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description="Minimum DP observed within the GVCF block">
##FORMAT=<ID=MLPSAC,Number=A,Type=Integer,Description="Maximum likelihood expectation (MLE) for the alternate allele count, in the same order as listed, for each individual sample">
##FORMAT=<ID=MLPSAF,Number=A,Type=Float,Description="Maximum likelihood expectation (MLE) for the alternate allele fraction, in the same order as listed, for each individual sample">
##FORMAT=<ID=MQ0,Number=1,Type=Integer,Description="Number of Mapping Quality Zero Reads per sample">
##FORMAT=<ID=PGT,Number=1,Type=String,Description="Physical phasing haplotype information, describing how the alternate alleles are phased in relation to one another">
##FORMAT=<ID=PID,Number=1,Type=String,Description="Physical phasing ID information, where each unique ID within a given sample (but not across samples) connects records within a phasing group">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
##FORMAT=<ID=SB,Number=4,Type=Integer,Description="Per-sample component statistics which comprise the Fisher's Exact Test to detect strand bias.">
##GVCFBlock=minGQ=0(inclusive),maxGQ=1(exclusive)
##INFO=<ID=ABHet,Number=1,Type=Float,Description="Allele Balance for heterozygous calls (ref/(ref+alt))">
##INFO=<ID=ABHom,Number=1,Type=Float,Description="Allele Balance for homozygous calls (A/(A+O)) where A is the allele (ref or alt) and O is anything other">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
##INFO=<ID=BaseCounts,Number=4,Type=Integer,Description="Counts of each base">
##INFO=<ID=BaseQRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt Vs. Ref base qualities">
##INFO=<ID=CCC,Number=1,Type=Integer,Description="Number of called chromosomes">
##INFO=<ID=ClippingRankSum,Number=1,Type=Float,Description="Z-score From Wilcoxon rank sum test of Alt vs. Ref number of hard clipped bases">
##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP Membership">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
##INFO=<ID=DS,Number=0,Type=Flag,Description="Were any of the samples downsampled?">
##INFO=<ID=END,Number=1,Type=Integer,Description="Stop position of the interval">
##INFO=<ID=FS,Number=1,Type=Float,Description="Phred-scaled p-value using Fisher's exact test to detect strand bias">
##INFO=<ID=GC,Number=1,Type=Integer,Description="GC content around the variant (see docs for window size details)">
##INFO=<ID=GQ_MEAN,Number=1,Type=Float,Description="Mean of all GQ values">
##INFO=<ID=GQ_STDDEV,Number=1,Type=Float,Description="Standard deviation of all GQ values">
##INFO=<ID=HRun,Number=1,Type=Integer,Description="Largest Contiguous Homopolymer Run of Variant Allele In Either Direction">
##INFO=<ID=HW,Number=1,Type=Float,Description="Phred-scaled p-value for Hardy-Weinberg violation">
##INFO=<ID=HWP,Number=1,Type=Float,Description="P value from test of Hardy Weinberg Equilibrium">
##INFO=<ID=HaplotypeScore,Number=1,Type=Float,Description="Consistency of the site with at most two segregating haplotypes">
##INFO=<ID=InbreedingCoeff,Number=1,Type=Float,Description="Inbreeding coefficient as estimated from the genotype likelihoods per-sample when compared against the Hardy-Weinberg expectation">
##INFO=<ID=LikelihoodRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt Vs. Ref haplotype likelihoods">
##INFO=<ID=LowMQ,Number=3,Type=Float,Description="3-tuple: <fraction of reads with MQ=0>,<fraction of reads with MQ<=10>,<total number of reads>">
##INFO=<ID=MLEAC,Number=A,Type=Integer,Description="Maximum likelihood expectation (MLE) for the allele counts (not necessarily the same as the AC), for each ALT allele, in the same order as listed">
##INFO=<ID=MLEAF,Number=A,Type=Float,Description="Maximum likelihood expectation (MLE) for the allele frequency (not necessarily the same as the AF), for each ALT allele, in the same order as listed">
##INFO=<ID=MQ,Number=1,Type=Float,Description="RMS Mapping Quality">
##INFO=<ID=MQ0,Number=1,Type=Integer,Description="Total Mapping Quality Zero Reads">
##INFO=<ID=MQRankSum,Number=1,Type=Float,Description="Z-score From Wilcoxon rank sum test of Alt vs. Ref read mapping qualities">
##INFO=<ID=MVLR,Number=1,Type=Float,Description="Mendelian violation likelihood ratio: L[MV] - L[No MV]">
##INFO=<ID=NCC,Number=1,Type=Integer,Description="Number of no-called samples">
##INFO=<ID=OND,Number=1,Type=Float,Description="Overall non-diploid ratio (alleles/(alleles+non-alleles))">
##INFO=<ID=PercentNBase,Number=1,Type=Float,Description="Percentage of N bases in the pileup">
##INFO=<ID=QD,Number=1,Type=Float,Description="Variant Confidence/Quality by Depth">
##INFO=<ID=ReadPosRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias">
##INFO=<ID=SNPEFF_AMINO_ACID_CHANGE,Number=1,Type=String,Description="Old/New amino acid for the highest-impact effect resulting from the current variant (in HGVS style)">
##INFO=<ID=SNPEFF_CODON_CHANGE,Number=1,Type=String,Description="Old/New codon for the highest-impact effect resulting from the current variant">
##INFO=<ID=SNPEFF_EFFECT,Number=1,Type=String,Description="The highest-impact effect resulting from the current variant (or one of the highest-impact effects, if there is a tie)">
##INFO=<ID=SNPEFF_EXON_ID,Number=1,Type=String,Description="Exon ID for the highest-impact effect resulting from the current variant">
##INFO=<ID=SNPEFF_FUNCTIONAL_CLASS,Number=1,Type=String,Description="Functional class of the highest-impact effect resulting from the current variant: [NONE, SILENT, MISSENSE, NONSENSE]">
##INFO=<ID=SNPEFF_GENE_BIOTYPE,Number=1,Type=String,Description="Gene biotype for the highest-impact effect resulting from the current variant">
##INFO=<ID=SNPEFF_GENE_NAME,Number=1,Type=String,Description="Gene name for the highest-impact effect resulting from the current variant">
##INFO=<ID=SNPEFF_IMPACT,Number=1,Type=String,Description="Impact of the highest-impact effect resulting from the current variant [MODIFIER, LOW, MODERATE, HIGH]">
##INFO=<ID=SNPEFF_TRANSCRIPT_ID,Number=1,Type=String,Description="Transcript ID for the highest-impact effect resulting from the current variant">
##INFO=<ID=SOR,Number=1,Type=Float,Description="Symmetric Odds Ratio of 2x2 contingency table to detect strand bias">
##INFO=<ID=Samples,Number=.,Type=String,Description="List of polymorphic samples">
##INFO=<ID=TDT,Number=A,Type=Float,Description="Test statistic from Wittkowski transmission disequilibrium test.">
##INFO=<ID=VariantType,Number=1,Type=String,Description="Variant type description">
##INFO=<ID=hiConfDeNovo,Number=1,Type=String,Description="High confidence possible de novo mutation (GQ >= 20 for all trio members)=[comma-delimited list of child samples]">
##INFO=<ID=loConfDeNovo,Number=1,Type=String,Description="Low confidence possible de novo mutation (GQ >= 10 for child, GQ > 0 for parents)=[comma-delimited list of child samples]">
##INFO=<ID=RANK,Number=.,Type=String,Description="Rank that indicates whether or not the variant is pathogenic (1 -> pathogenic, 0 -> benign)">
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
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO  FORMAT  NA12878 NA12891 NA12892\n"""

with open(sys.argv[2], 'w', newline='') as outfile:
    outfile.write(header)
    count = 0
    for index, row in data.iterrows():
        print(count)
        outfile.write(str(row['#Chr']) + '\t' + 
                      str(row['#Position']) + '\t' + 
                      str(row['#dbsnpIdentifier']) + '\t' + 
                      str(row['#Reference']) + '\t' + 
                      str(row['#Alteration']) + '\t' + 
                      '.' + '\t' + 
                      '.' + '\t' + 
                      'AF=' + str(row['#AlleleFrequency'])[:6] + ';' + 'RANK=' + str(row['Rank'])[:6] + '\t' +
                      'GT:AB:AD:DP:GQ' + '\t' + 
                      str(row['NA12878'])[:6] + ':' + '.' + ':' + str(row['REF']) + ',' + str(row['ALT']) + ':' + str(row['DP']) + ':' + str(row['GQ']) + '\t' +
                      str(row['NA12891'])[:6] + ':' + '.' + ':' + str(row['REF.1']) + ',' + str(row['ALT.1']) + ':' + str(row['DP.1']) + ':' + str(row['GQ.1']) + '\t' +
                      str(row['NA12892'])[:6] + ':' + '.' + ':' + str(row['REF.2']) + ',' + str(row['ALT.2']) + ':' + str(row['DP.2']) + ':' + str(row['GQ.2']) + '\t' +
                      '\n')
        count += 1

print('Finished')
