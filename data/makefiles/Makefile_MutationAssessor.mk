help:
	@cat Makefile_MutationAssessor.mk

download:
	wget -c http://mutationassessor.org/r3/MA_scores_rel3_hg19_full.tar.gz

convert:
	tar -xzf MA_scores_rel3_hg19_full.tar.gz
	python3 create_MutationAssessor_vcf-file.py $$(ls -m MA_scores_rel3_hg19_full/*.csv | tr -d '[:space:]') hg19_precomputed_MutationAssessor.vcf
	bgzip hg19_precomputed_MutationAssessor.vcf
	tabix -p vcf hg19_precomputed_MutationAssessor.vcf.gz
	rm MA_scores_rel3_hg19_full.tar.gz
	rm -r MA_scores_rel3_hg19_full/
