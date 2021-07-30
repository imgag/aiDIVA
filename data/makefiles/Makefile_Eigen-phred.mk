help:
	@cat Makefile_Eigen-phred.mk

download:
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

convert:
	python3 create_eigen_vcf-files.py $$(ls -m *.tab.chr*.gz | tr -d '[:space:]') hg19_Eigen-phred_coding_chrom1-22.vcf hg19_Eigen-phred_noncoding_chrom1-22.vcf
	bgzip hg19_Eigen-phred_coding_chrom1-22.vcf
	bgzip hg19_Eigen-phred_noncoding_chrom1-22.vcf
	tabix -p vcf hg19_Eigen-phred_coding_chrom1-22.vcf.gz
	tabix -p vcf hg19_Eigen-phred_noncoding_chrom1-22.vcf.gz
	rm *.tab.chr*.gz
