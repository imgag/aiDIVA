help:
	@cat Makefile_Condel.mk

download:
	wget -c https://bbglab.irbbarcelona.org/fannsdb/downloads/fannsdb.tsv.gz

convert:
	python3 create_Condel_vcf-file.py fannsdb.tsv.gz hg19_precomputed_Condel.vcf
	bgzip hg19_precomputed_Condel.vcf
	tabix -p vcf hg19_precomputed_Condel.vcf.gz
	rm fannsdb.tsv.gz
