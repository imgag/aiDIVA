help:
	@cat Makefile_simpleRepeat.mk

download_hg19:
	wget -O hg19.simpleRepeat.txt.gz ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/simpleRepeat.txt.gz

download_hg38:
	wget -O hg38.simpleRepeat.txt.gz ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/simpleRepeat.txt.gz

convert_hg19:
	gunzip hg19.simpleRepeat.txt.gz
	cut -f2,3,4,11 hg19.simpleRepeat.txt > hg19.simpleRepeat.bed
	grep -v '#' hg19.simpleRepeat.bed | sort -k1,1 -k2,2n -k3,3n -t '	' | bgzip -c > hg19.simpleRepeat.bed.gz
	tabix -p bed hg19.simpleRepeat.bed.gz

convert_hg38:
	gunzip hg38.simpleRepeat.txt.gz
	cut -f2,3,4,11 hg38.simpleRepeat.txt > hg38.simpleRepeat.bed
	grep -v '#' hg38.simpleRepeat.bed | sort -k1,1 -k2,2n -k3,3n -t '	' | bgzip -c > hg38.simpleRepeat.bed.gz
	tabix -p bed hg38.simpleRepeat.bed.gz
