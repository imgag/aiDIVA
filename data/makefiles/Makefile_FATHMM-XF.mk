help:
	@cat Makefile_FATHMM-XF.mk

download:
	wget -c http://fathmm.biocompute.org.uk/fathmm-xf/fathmm_xf_coding.vcf.gz

convert:
	python3 create_fathmm-XF_vcf.py fathmm_xf_coding.vcf.gz hg19_fathmm_xf_coding.vcf
	bgzip hg19_fathmm_xf_coding.vcf
	tabix -p vcf hg19_fathmm_xf_coding.vcf.gz
	rm fathmm_xf_coding.vcf.gz
