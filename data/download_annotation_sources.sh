#!/bin/bash

folder=`pwd`
annotation_sources=$folder/annotation_resources/
mkdir -p $annotation_sources

# Install REVEL for VEP - https://sites.google.com/site/revelgenomics/downloads
cd $annotation_sources
mkdir -p REVEL
cd REVEL
wget -c https://rothsj06.u.hpc.mssm.edu/revel/revel_all_chromosomes.csv.zip
unzip -p revel_all_chromosomes.csv.zip | tr ',' '\t' | sed '1s/.*/#&/' | bgzip > revel_all_chromosomes.tsv.gz
tabix -f -s 1 -b 2 -e 2 revel_all_chromosomes.tsv.gz

# Install fathmm-XF for AIdiva (custom VEP annotation) - https://fathmm.biocompute.org.uk/fathmm-xf
cd $annotation_sources
mkdir -p fathmm-XF
cd fathmm-XF
wget -c http://fathmm.biocompute.org.uk/fathmm-xf/fathmm_xf_coding.vcf.gz
python3 $folder/create_fathmm-XF_vcf.py fathmm_xf_coding.vcf.gz hg19_fathmm_xf_coding.vcf
bgzip hg19_fathmm_xf_coding.vcf
tabix -p vcf hg19_fathmm_xf_coding.vcf.gz
rm fathmm_xf_coding.vcf.gz

# Install ABB score for AIdiva (custom VEP annotation) - https://github.com/Francesc-Muyas/ABB
cd $annotation_sources
mkdir -p ABB
cd ABB
wget -c https://public_docs.crg.es/sossowski/publication_data/ABB/ABB_SCORE.txt
python3 $folder/create_ABB-SCORE_bed.py ABB_SCORE.txt hg19_ABB-SCORE.bed
bgzip hg19_ABB-SCORE.bed
tabix -p bed hg19_ABB-SCORE.bed.gz
rm ABB_SCORE.txt

# Install Eigen phred for AIdiva (custom VEP annotation)
cd $annotation_sources
mkdir -p Eigen
cd Eigen
wget -c http://web.corral.tacc.utexas.edu/WGSAdownload/resources/Eigen/Eigen_hg19_combined.tab.chr1.gz
wget -c http://web.corral.tacc.utexas.edu/WGSAdownload/resources/Eigen/Eigen_hg19_combined.tab.chr2.gz
wget -c http://web.corral.tacc.utexas.edu/WGSAdownload/resources/Eigen/Eigen_hg19_combined.tab.chr3.gz
wget -c http://web.corral.tacc.utexas.edu/WGSAdownload/resources/Eigen/Eigen_hg19_combined.tab.chr4.gz
wget -c http://web.corral.tacc.utexas.edu/WGSAdownload/resources/Eigen/Eigen_hg19_combined.tab.chr5.gz
wget -c http://web.corral.tacc.utexas.edu/WGSAdownload/resources/Eigen/Eigen_hg19_combined.tab.chr6.gz
wget -c http://web.corral.tacc.utexas.edu/WGSAdownload/resources/Eigen/Eigen_hg19_combined.tab.chr7.gz
wget -c http://web.corral.tacc.utexas.edu/WGSAdownload/resources/Eigen/Eigen_hg19_combined.tab.chr8.gz
wget -c http://web.corral.tacc.utexas.edu/WGSAdownload/resources/Eigen/Eigen_hg19_combined.tab.chr9.gz
wget -c http://web.corral.tacc.utexas.edu/WGSAdownload/resources/Eigen/Eigen_hg19_combined.tab.chr10.gz
wget -c http://web.corral.tacc.utexas.edu/WGSAdownload/resources/Eigen/Eigen_hg19_combined.tab.chr11.gz
wget -c http://web.corral.tacc.utexas.edu/WGSAdownload/resources/Eigen/Eigen_hg19_combined.tab.chr12.gz
wget -c http://web.corral.tacc.utexas.edu/WGSAdownload/resources/Eigen/Eigen_hg19_combined.tab.chr13.gz
wget -c http://web.corral.tacc.utexas.edu/WGSAdownload/resources/Eigen/Eigen_hg19_combined.tab.chr14.gz
wget -c http://web.corral.tacc.utexas.edu/WGSAdownload/resources/Eigen/Eigen_hg19_combined.tab.chr15.gz
wget -c http://web.corral.tacc.utexas.edu/WGSAdownload/resources/Eigen/Eigen_hg19_combined.tab.chr16.gz
wget -c http://web.corral.tacc.utexas.edu/WGSAdownload/resources/Eigen/Eigen_hg19_combined.tab.chr17.gz
wget -c http://web.corral.tacc.utexas.edu/WGSAdownload/resources/Eigen/Eigen_hg19_combined.tab.chr18.gz
wget -c http://web.corral.tacc.utexas.edu/WGSAdownload/resources/Eigen/Eigen_hg19_combined.tab.chr19.gz
wget -c http://web.corral.tacc.utexas.edu/WGSAdownload/resources/Eigen/Eigen_hg19_combined.tab.chr20.gz
wget -c http://web.corral.tacc.utexas.edu/WGSAdownload/resources/Eigen/Eigen_hg19_combined.tab.chr21.gz
wget -c http://web.corral.tacc.utexas.edu/WGSAdownload/resources/Eigen/Eigen_hg19_combined.tab.chr22.gz
python3 $folder/create_eigen_vcf.py $(ls -m *.tab.chr*.gz | tr -d '[:space:]') hg19_Eigen-phred_coding_chrom1-22.vcf
#python3 $src/Tools/create_eigen_vcf-files.py $(ls -m *.tab.chr*.gz | tr -d '[:space:]') hg19_Eigen-phred_coding_chrom1-22.vcf hg19_Eigen-phred_noncoding_chrom1-22.vcf
bgzip hg19_Eigen-phred_coding_chrom1-22.vcf
#bgzip hg19_Eigen-phred_noncoding_chrom1-22.vcf
tabix -p vcf hg19_Eigen-phred_coding_chrom1-22.vcf.gz
#tabix -p vcf hg19_Eigen-phred_noncoding_chrom1-22.vcf.gz
rm *.tab.chr*.gz

# Install Condel for AIdiva (custom VEP annotation)
cd $annotation_sources
mkdir -p Condel
cd Condel
wget -c https://bbglab.irbbarcelona.org/fannsdb/downloads/fannsdb.tsv.gz
python3 $folder/create_Condel_vcf.py fannsdb.tsv.gz hg19_precomputed_Condel.vcf
bgzip hg19_precomputed_Condel.vcf
tabix -p vcf hg19_precomputed_Condel.vcf.gz
rm fannsdb.tsv.gz

# Install MutationAssessor for AIdiva (custom VEP annotation)
cd $annotation_sources
mkdir -p MutationAssessor
cd MutationAssessor
wget -c http://mutationassessor.org/r3/MA_scores_rel3_hg19_full.tar.gz
tar -xzf MA_scores_rel3_hg19_full.tar.gz
python3 $folder/create_MutationAssessor_vcf.py $(ls -m MA_scores_rel3_hg19_full/*.csv | tr -d '[:space:]') hg19_precomputed_MutationAssessor_unsort.vcf
sort -k1,1 -k2,2n hg19_precomputed_MutationAssessor_unsort.vcf > hg19_precomputed_MutationAssessor.vcf
bgzip hg19_precomputed_MutationAssessor.vcf
tabix -p vcf hg19_precomputed_MutationAssessor.vcf.gz
rm MA_scores_rel3_hg19_full.tar.gz
rm hg19_precomputed_MutationAssessor_unsort.vcf
rm -r MA_scores_rel3_hg19_full/

# Install segment duplication and simple repeats for AIdiva (custom VEP annotation)
cd $annotation_sources
mkdir -p UCSC
cd UCSC
wget -O hg19_genomicSuperDups.txt.gz ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/genomicSuperDups.txt.gz
gunzip hg19_genomicSuperDups.txt.gz
cut -f2,3,4,27 hg19_genomicSuperDups.txt > hg19_genomicSuperDups.bed
grep -v "#" hg19_genomicSuperDups.bed | sort -k1,1 -k2,2n -k3,3n -t$'\t' | bgzip -c > hg19_genomicSuperDups.bed.gz
tabix -p bed hg19_genomicSuperDups.bed.gz
wget -O hg19_simpleRepeat.txt.gz ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/simpleRepeat.txt.gz
gunzip hg19_simpleRepeat.txt.gz
cut -f2,3,4,11 hg19_simpleRepeat.txt > hg19_simpleRepeat.bed
grep -v '#' hg19_simpleRepeat.bed | sort -k1,1 -k2,2n -k3,3n -t$'\t' | bgzip -c > hg19_simpleRepeat.bed.gz
tabix -p bed hg19_simpleRepeat.bed.gz

# Install phastCons46way vertebrate for AIdiva (custom VEP annotation)
cd $annotation_sources
mkdir -p phastCons
cd phastCons
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr1.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr10.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr11.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr11_gl000202_random.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr12.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr13.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr14.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr15.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr16.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr17.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr17_ctg5_hap1.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr17_gl000203_random.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr17_gl000204_random.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr17_gl000205_random.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr17_gl000206_random.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr18.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr18_gl000207_random.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr19.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr19_gl000208_random.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr19_gl000209_random.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr1_gl000191_random.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr1_gl000192_random.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr2.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr20.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr21.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr21_gl000210_random.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr22.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr3.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr4.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr4_ctg9_hap1.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr4_gl000193_random.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr4_gl000194_random.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr5.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr6.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr6_apd_hap1.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr6_cox_hap2.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr6_dbb_hap3.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr6_mann_hap4.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr6_mcf_hap5.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr6_qbl_hap6.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr6_ssto_hap7.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr7.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr7_gl000195_random.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr8.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr8_gl000196_random.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr8_gl000197_random.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr9.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr9_gl000198_random.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr9_gl000199_random.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr9_gl000200_random.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr9_gl000201_random.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chrM.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chrUn_gl000211.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chrUn_gl000212.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chrUn_gl000213.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chrUn_gl000214.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chrUn_gl000215.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chrUn_gl000216.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chrUn_gl000217.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chrUn_gl000218.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chrUn_gl000219.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chrUn_gl000220.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chrUn_gl000221.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chrUn_gl000222.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chrUn_gl000223.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chrUn_gl000224.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chrUn_gl000225.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chrUn_gl000226.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chrUn_gl000227.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chrUn_gl000228.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chrUn_gl000229.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chrUn_gl000230.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chrUn_gl000231.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chrUn_gl000232.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chrUn_gl000233.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chrUn_gl000234.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chrUn_gl000235.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chrUn_gl000236.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chrUn_gl000237.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chrUn_gl000238.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chrUn_gl000239.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chrUn_gl000240.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chrUn_gl000241.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chrUn_gl000242.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chrUn_gl000243.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chrUn_gl000244.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chrUn_gl000245.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chrUn_gl000246.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chrUn_gl000247.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chrUn_gl000248.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chrUn_gl000249.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chrX.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chrY.phastCons46way.wigFix.gz
wget -c http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes
wget -c http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/wigToBigWig
wget -c http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigWigCat
chmod +x wigToBigWig
chmod +x bigWigCat
gunzip *.wigFix.gz
for f in *.wigFix; do ./wigToBigWig -fixedSummaries -keepAllChromosomes $f hg19.chrom.sizes $(basename $f ".wigFix").bw; done;
rm *.wigFix
./bigWigCat hg19_phastCons46way_vertebrate.bw *.phastCons46way.bw
rm *.phastCons46way.bw

# Install phastCons46way primate for AIdiva (custom VEP annotation)
cd $annotation_sources
mkdir -p phastCons
cd phastCons
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr1.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr10.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr11.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr11_gl000202_random.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr12.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr13.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr14.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr15.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr16.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr17.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr17_ctg5_hap1.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr17_gl000203_random.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr17_gl000204_random.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr17_gl000205_random.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr17_gl000206_random.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr18.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr18_gl000207_random.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr19.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr19_gl000208_random.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr19_gl000209_random.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr1_gl000191_random.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr1_gl000192_random.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr2.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr20.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr21.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr21_gl000210_random.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr22.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr3.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr4.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr4_ctg9_hap1.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr4_gl000193_random.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr4_gl000194_random.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr5.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr6.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr6_apd_hap1.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr6_cox_hap2.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr6_dbb_hap3.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr6_mann_hap4.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr6_mcf_hap5.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr6_qbl_hap6.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr6_ssto_hap7.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr7.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr7_gl000195_random.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr8.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr8_gl000196_random.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr8_gl000197_random.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr9.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr9_gl000198_random.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr9_gl000199_random.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr9_gl000200_random.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr9_gl000201_random.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chrM.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chrUn_gl000211.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chrUn_gl000212.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chrUn_gl000213.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chrUn_gl000214.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chrUn_gl000215.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chrUn_gl000216.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chrUn_gl000217.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chrUn_gl000218.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chrUn_gl000219.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chrUn_gl000220.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chrUn_gl000221.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chrUn_gl000222.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chrUn_gl000223.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chrUn_gl000224.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chrUn_gl000225.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chrUn_gl000226.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chrUn_gl000227.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chrUn_gl000228.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chrUn_gl000229.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chrUn_gl000230.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chrUn_gl000231.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chrUn_gl000232.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chrUn_gl000233.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chrUn_gl000234.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chrUn_gl000235.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chrUn_gl000236.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chrUn_gl000237.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chrUn_gl000238.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chrUn_gl000239.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chrUn_gl000240.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chrUn_gl000241.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chrUn_gl000242.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chrUn_gl000243.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chrUn_gl000244.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chrUn_gl000245.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chrUn_gl000246.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chrUn_gl000247.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chrUn_gl000248.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chrUn_gl000249.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chrX.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chrY.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes
wget -c http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/wigToBigWig
wget -c http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigWigCat
chmod +x wigToBigWig
chmod +x bigWigCat
gunzip *.wigFix.gz
for f in *.wigFix; do ./wigToBigWig -fixedSummaries -keepAllChromosomes $f hg19.chrom.sizes $(basename $f ".wigFix").bw; done;
rm *.wigFix
./bigWigCat hg19_phastCons46way_primate.bw *.phastCons46way.primates.bw
rm *.phastCons46way.primates.bw

# Install phastCons46way mammal for AIdiva (custom VEP annotation)
cd $annotation_sources
mkdir -p phastCons
cd phastCons
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr1.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr10.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr11.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr11_gl000202_random.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr12.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr13.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr14.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr15.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr16.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr17.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr17_ctg5_hap1.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr17_gl000203_random.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr17_gl000204_random.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr17_gl000205_random.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr17_gl000206_random.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr18.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr18_gl000207_random.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr19.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr19_gl000208_random.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr19_gl000209_random.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr1_gl000191_random.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr1_gl000192_random.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr2.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr20.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr21.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr21_gl000210_random.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr22.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr3.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr4.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr4_ctg9_hap1.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr4_gl000193_random.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr4_gl000194_random.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr5.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr6.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr6_apd_hap1.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr6_cox_hap2.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr6_dbb_hap3.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr6_mann_hap4.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr6_mcf_hap5.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr6_qbl_hap6.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr6_ssto_hap7.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr7.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr7_gl000195_random.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr8.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr8_gl000196_random.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr8_gl000197_random.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr9.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr9_gl000198_random.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr9_gl000199_random.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr9_gl000200_random.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr9_gl000201_random.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chrM.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chrUn_gl000211.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chrUn_gl000212.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chrUn_gl000213.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chrUn_gl000214.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chrUn_gl000215.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chrUn_gl000216.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chrUn_gl000217.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chrUn_gl000218.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chrUn_gl000219.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chrUn_gl000220.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chrUn_gl000221.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chrUn_gl000222.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chrUn_gl000223.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chrUn_gl000224.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chrUn_gl000225.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chrUn_gl000226.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chrUn_gl000227.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chrUn_gl000228.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chrUn_gl000229.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chrUn_gl000230.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chrUn_gl000231.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chrUn_gl000232.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chrUn_gl000233.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chrUn_gl000234.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chrUn_gl000235.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chrUn_gl000236.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chrUn_gl000237.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chrUn_gl000238.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chrUn_gl000239.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chrUn_gl000240.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chrUn_gl000241.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chrUn_gl000242.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chrUn_gl000243.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chrUn_gl000244.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chrUn_gl000245.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chrUn_gl000246.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chrUn_gl000247.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chrUn_gl000248.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chrUn_gl000249.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chrX.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chrY.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes
wget -c http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/wigToBigWig
wget -c http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigWigCat
chmod +x wigToBigWig
chmod +x bigWigCat
gunzip *.wigFix.gz
for f in *.wigFix; do ./wigToBigWig -fixedSummaries -keepAllChromosomes $f hg19.chrom.sizes $(basename $f ".wigFix").bw; done;
rm *.wigFix
./bigWigCat hg19_phastCons46way_mammal.bw *.phastCons46way.placental.bw
rm *.phastCons46way.placental.bw

# Install phyloP46way primate for AIdiva (custom VEP annotation)
cd $annotation_sources
mkdir -p phyloP
cd phyloP
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr1.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr10.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr11.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr11_gl000202_random.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr12.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr13.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr14.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr15.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr16.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr17.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr17_ctg5_hap1.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr17_gl000203_random.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr17_gl000204_random.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr17_gl000205_random.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr17_gl000206_random.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr18.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr18_gl000207_random.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr19.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr19_gl000208_random.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr19_gl000209_random.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr1_gl000191_random.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr1_gl000192_random.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr2.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr20.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr21.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr21_gl000210_random.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr22.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr3.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr4.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr4_ctg9_hap1.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr4_gl000193_random.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr4_gl000194_random.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr5.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr6.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr6_apd_hap1.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr6_cox_hap2.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr6_dbb_hap3.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr6_mann_hap4.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr6_mcf_hap5.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr6_qbl_hap6.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr6_ssto_hap7.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr7.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr7_gl000195_random.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr8.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr8_gl000196_random.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr8_gl000197_random.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr9.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr9_gl000198_random.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr9_gl000199_random.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr9_gl000200_random.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr9_gl000201_random.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chrM.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chrUn_gl000211.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chrUn_gl000212.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chrUn_gl000213.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chrUn_gl000214.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chrUn_gl000215.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chrUn_gl000216.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chrUn_gl000217.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chrUn_gl000218.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chrUn_gl000219.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chrUn_gl000220.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chrUn_gl000221.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chrUn_gl000222.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chrUn_gl000223.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chrUn_gl000224.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chrUn_gl000225.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chrUn_gl000227.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chrUn_gl000228.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chrUn_gl000229.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chrUn_gl000230.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chrUn_gl000231.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chrUn_gl000232.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chrUn_gl000233.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chrUn_gl000234.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chrUn_gl000235.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chrUn_gl000236.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chrUn_gl000237.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chrUn_gl000238.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chrUn_gl000239.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chrUn_gl000240.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chrUn_gl000241.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chrUn_gl000242.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chrUn_gl000243.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chrUn_gl000244.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chrUn_gl000245.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chrUn_gl000246.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chrUn_gl000247.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chrUn_gl000248.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chrUn_gl000249.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chrX.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chrY.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes
wget -c http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/wigToBigWig
wget -c http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigWigCat
chmod +x wigToBigWig
chmod +x bigWigCat
gunzip *.wigFix.gz
for f in *.wigFix; do ./wigToBigWig -fixedSummaries -keepAllChromosomes $f hg19.chrom.sizes $(basename $f ".wigFix").bw; done;
rm *.wigFix
./bigWigCat hg19_phyloP46way_primate.bw *.phyloP46way.primate.bw
rm *.phyloP46way.primate.bw

# Install phyloP46way mammal for AIdiva (custom VEP annotation)
cd $annotation_sources
mkdir -p phyloP
cd phyloP
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr1.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr10.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr11.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr11_gl000202_random.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr12.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr13.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr14.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr15.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr16.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr17.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr17_ctg5_hap1.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr17_gl000203_random.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr17_gl000204_random.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr17_gl000205_random.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr17_gl000206_random.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr18.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr18_gl000207_random.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr19.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr19_gl000208_random.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr19_gl000209_random.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr1_gl000191_random.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr1_gl000192_random.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr2.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr20.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr21.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr21_gl000210_random.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr22.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr3.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr4.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr4_ctg9_hap1.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr4_gl000193_random.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr4_gl000194_random.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr5.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr6.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr6_apd_hap1.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr6_cox_hap2.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr6_dbb_hap3.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr6_mann_hap4.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr6_mcf_hap5.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr6_qbl_hap6.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr6_ssto_hap7.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr7.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr7_gl000195_random.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr8.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr8_gl000196_random.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr8_gl000197_random.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr9.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr9_gl000198_random.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr9_gl000199_random.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr9_gl000200_random.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr9_gl000201_random.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chrM.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chrUn_gl000211.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chrUn_gl000212.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chrUn_gl000213.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chrUn_gl000214.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chrUn_gl000215.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chrUn_gl000216.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chrUn_gl000217.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chrUn_gl000218.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chrUn_gl000219.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chrUn_gl000220.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chrUn_gl000221.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chrUn_gl000222.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chrUn_gl000223.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chrUn_gl000224.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chrUn_gl000225.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chrUn_gl000227.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chrUn_gl000228.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chrUn_gl000229.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chrUn_gl000230.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chrUn_gl000231.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chrUn_gl000232.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chrUn_gl000233.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chrUn_gl000234.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chrUn_gl000235.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chrUn_gl000236.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chrUn_gl000237.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chrUn_gl000238.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chrUn_gl000239.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chrUn_gl000240.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chrUn_gl000241.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chrUn_gl000242.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chrUn_gl000243.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chrUn_gl000244.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chrUn_gl000245.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chrUn_gl000246.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chrUn_gl000247.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chrUn_gl000248.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chrUn_gl000249.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chrX.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chrY.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes
wget -c http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/wigToBigWig
wget -c http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigWigCat
chmod +x wigToBigWig
chmod +x bigWigCat
gunzip *.wigFix.gz
for f in *.wigFix; do ./wigToBigWig -fixedSummaries -keepAllChromosomes $f hg19.chrom.sizes $(basename $f ".wigFix").bw; done;
rm *.wigFix
./bigWigCat hg19_phyloP46way_mammal.bw *.phyloP46way.placental.bw
rm *.phyloP46way.placental.bw

# Install phyloP46way vertebrate for AIdiva (custom VEP annotation)
cd $annotation_sources
mkdir -p phyloP
cd phyloP
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr1.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr10.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr11.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr11_gl000202_random.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr12.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr13.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr14.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr15.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr16.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr17.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr17_ctg5_hap1.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr17_gl000203_random.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr17_gl000204_random.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr17_gl000205_random.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr17_gl000206_random.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr18.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr18_gl000207_random.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr19.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr19_gl000208_random.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr19_gl000209_random.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr1_gl000191_random.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr1_gl000192_random.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr2.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr20.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr21.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr21_gl000210_random.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr22.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr3.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr4.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr4_ctg9_hap1.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr4_gl000193_random.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr4_gl000194_random.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr5.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr6.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr6_apd_hap1.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr6_cox_hap2.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr6_dbb_hap3.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr6_mann_hap4.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr6_mcf_hap5.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr6_qbl_hap6.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr6_ssto_hap7.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr7.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr7_gl000195_random.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr8.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr8_gl000196_random.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr8_gl000197_random.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr9.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr9_gl000198_random.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr9_gl000199_random.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr9_gl000200_random.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr9_gl000201_random.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chrM.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chrUn_gl000211.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chrUn_gl000212.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chrUn_gl000213.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chrUn_gl000214.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chrUn_gl000215.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chrUn_gl000216.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chrUn_gl000217.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chrUn_gl000218.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chrUn_gl000219.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chrUn_gl000220.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chrUn_gl000221.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chrUn_gl000222.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chrUn_gl000223.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chrUn_gl000224.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chrUn_gl000225.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chrUn_gl000227.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chrUn_gl000228.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chrUn_gl000229.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chrUn_gl000230.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chrUn_gl000231.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chrUn_gl000232.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chrUn_gl000233.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chrUn_gl000234.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chrUn_gl000235.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chrUn_gl000236.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chrUn_gl000237.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chrUn_gl000238.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chrUn_gl000239.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chrUn_gl000240.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chrUn_gl000241.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chrUn_gl000242.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chrUn_gl000243.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chrUn_gl000244.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chrUn_gl000245.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chrUn_gl000246.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chrUn_gl000247.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chrUn_gl000248.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chrUn_gl000249.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chrX.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chrY.phyloP46way.wigFix.gz
wget -c http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes
wget -c http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/wigToBigWig
wget -c http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigWigCat
chmod +x wigToBigWig
chmod +x bigWigCat
gunzip *.wigFix.gz
for f in *.wigFix; do ./wigToBigWig -fixedSummaries -keepAllChromosomes $f hg19.chrom.sizes $(basename $f ".wigFix").bw; done;
rm *.wigFix
./bigWigCat hg19_phyloP46way_vertebrate.bw *.phyloP46way.bw
rm *.phyloP46way.bw
