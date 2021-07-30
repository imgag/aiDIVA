help:
	@cat Makefile_segmentDuplication.mk

download_hg19:
	wget -O hg19.genomicSuperDups.txt.gz ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/genomicSuperDups.txt.gz

download:
	wget -O hg38.genomicSuperDups.txt.gz ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/genomicSuperDups.txt.gz

convert_hg19:
	gunzip hg19.genomicSuperDups.txt.gz
	cut -f2,3,4,27 hg19.genomicSuperDups.txt > hg19.genomicSuperDups.bed
	grep -v '#' hg19.genomicSuperDups.bed | sort -k1,1 -k2,2n -k3,3n -t '	' | bgzip -c > hg19.genomicSuperDups.bed.gz
	tabix -p bed hg19.genomicSuperDups.bed.gz

convert_hg38:
	gunzip hg38.genomicSuperDups.txt.gz
	cut -f2,3,4,27 hg38.genomicSuperDups.txt > hg38.genomicSuperDups.bed
	grep -v '#' hg38.genomicSuperDups.bed | sort -k1,1 -k2,2n -k3,3n -t '	' | bgzip -c > hg38.genomicSuperDups.bed.gz
	tabix -p bed hg38.genomicSuperDups.bed.gz
