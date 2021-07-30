help:
	@cat Makefile_gnomAD_OE.mk

download:
	wget -c https://storage.googleapis.com/gnomad-public/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz

convert:
	zcat gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz | sed -n '1d;p' | awk -v OFS="\t" -F "\t" '{ print $$75,$$76,$$77,$$24}' > gnomad_OE.bed
	cat gnomad_OE.bed | sort -k1,1 -k2,2n -k3,3n -t '	' | sed '/NA/s/\bNA//g' | bgzip -c > gnomAD_OE_sorted.bed.gz
	tabix -p bed gnomAD_OE_sorted.bed.gz