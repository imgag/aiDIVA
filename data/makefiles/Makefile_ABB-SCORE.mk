help:
	@cat Makefile_ABB-SCORE.mk

download:
	wget -c https://public_docs.crg.es/sossowski/publication_data/ABB/ABB_SCORE.txt

convert:
	python3 create_ABB-SCORE_bed.py ABB_SCORE.txt hg19_ABB-SCORE.bed
	bgzip hg19_ABB-SCORE.bed
	tabix -p bed hg19_ABB-SCORE.bed.gz
	rm ABB_SCORE.txt
