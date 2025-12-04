# Preparation of Annotation Resources
This document provides links to the annotation sources we used. If it is necessary to further prepare these files it is explained and shown in the respective section with a coode snippet.

The preparation of the annotation resources needs a lot of disk space, you should make sure that you have at least 1000GB of free space available. Furthermore it is adviced to prepare one resource after another and to remove all files that are not needed after preparation to save disk space.

## Necessary Tools
Make sure that the following tools are installed on your system:
- samtools [https://www.htslib.org/doc/samtools.html](https://www.htslib.org/doc/samtools.html)
- bcftools [https://www.htslib.org/doc/bcftools.html](https://www.htslib.org/doc/bcftools.html)
- tabix [https://www.htslib.org/doc/tabix.html](https://www.htslib.org/doc/tabix.html)
- bgzip [https://www.htslib.org/doc/bgzip.html](https://www.htslib.org/doc/bgzip.html)

<br>

For some of the annotation sources there are no GRCh38/hg38 files available in these cases we lifted the GRCh37/hg19 sources manually to the new assembly using the tool [CrossMap](https://github.com/liguowang/CrossMap). For this manual liftover we also provide code snippets.


## Reference Genome Assembly
GRCh38:
<br>
[http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz](http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz)

```
gunzip hg38.fa.gz

samtools dict hg38.fa > hg38.dict
samtools faidx hg38.fa
```

## Annotation Databases
INFO: Please be adviced that some of the following annotation sources are only free for non-commercial and academic use!

If for some reason you choose to exclude one or more of the annotation sources shown here you have to make sure to modify the source code accordingly. If an annotation source is used as a feature for the random forest (RF) and you exclude it in your analysis you have to train a new model before you can use aiDIVA-RF.


### CADD
GRCh38:
<br>
https://kircherlab.bihealth.org/download/CADD/v1.6/GRCh38/whole_genome_SNVs.tsv.gz

```
wget -c -O grch38_whole_genome_SNVs.tsv.gz https://kircherlab.bihealth.org/download/CADD/v1.6/GRCh38/whole_genome_SNVs.tsv.gz

python3 prepare_CADD_vcf.py grch38_whole_genome_SNVs.tsv.gz grch38_CADD_snvs_v16.vcf

bgzip grch38_CADD_snvs_v16.vcf
tabix -p vcf grch38_CADD_snvs_v16.vcf.gz

rm grch38_whole_genome_SNVs.tsv.gz
```


### Capice
GRCh38:
<br>
https://zenodo.org/record/3928295/files/capice_v1.0_build37_snvs.tsv.gz

Needs manual liftover!!!
https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz

```
wget -c https://zenodo.org/record/3928295/files/capice_v1.0_build37_snvs.tsv.gz

python3 prepare_Capice_vcf.py capice_v1.0_build37_snvs.tsv.gz grch37_capice_v1_snvs.vcf

CrossMap vcf --chromid l hg19ToHg38.over.chain.gz grch37_capice_v1_snvs.vcf hg38.fa grch38_capice_v1_snvs_unsort.vcf

# Make sure that the created VCF is sorted
cat grch38_capice_v1_snvs_unsort.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > grch38_capice_v1_snvs.vcf

bgzip grch38_capice_v1_snvs.vcf
tabix -p vcf grch38_capice_v1_snvs.vcf.gz

rm grch37_capice_v1_snvs.vcf.gz
rm capice_v1.0_build37_snvs.tsv.gz
```


### ClinVar
GRCh38:
<br>
https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz

```
wget -c https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz

python3 prepare_ClinVar_vcf.py clinvar.vcf.gz grch38_clinvar.vcf

bgzip grch38_clinvar.vcf
tabix -p vcf grch38_clinvar.vcf.gz

rm clinvar.vcf.gz
```


### Condel
!!! Account needed !!!
<br>
<br>
Please first download the FannsDB database from here: [https://bbglab.irbbarcelona.org/fannsdb/](https://bbglab.irbbarcelona.org/fannsdb/) 
<br>
<br>
Afterwards the `prepare_Condel_vcf.py` script found in the _scripts_ folder can be used to create the VCF annotation file. For the annotation the VCF file needs to be _bgzipped_ and _indexed_

Needs manual liftover!!!
https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz

The following code snippet shows the necessary commands

GRCh38:

```
python3 prepare_Condel_vcf.py fannsdb.tsv.gz grch37_fannsDB_Condel.vcf

bgzip grch37_fannsDB_Condel.vcf
tabix -p vcf grch37_fannsDB_Condel.vcf.gz

CrossMap vcf --chromid l hg19ToHg38.over.chain.gz grch37_fannsDB_Condel.vcf hg38.fa grch38_fannsDB_Condel_unsort.vcf

# Make sure that the created VCF is sorted
cat grch38_fannsDB_Condel_unsort.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > grch38_fannsDB_Condel.vcf

bgzip grch38_fannsDB_Condel.vcf
tabix -p vcf grch38_fannsDB_Condel.vcf.gz

rm grch37_fannsDB_Condel.vcf.gz
rm fannsdb.tsv.gz
```


### Eigen phred
GRCh38:
<br>
https://zenodo.org/record/3256671/files/Eigen_hg38_coding_annot_04092016.tab.gz

```
wget -c https://zenodo.org/record/3256671/files/Eigen_hg38_coding_annot_04092016.tab.gz

python3 prepare_EigenPhred_vcf.py Eigen_hg38_coding_annot_04092016.tab.gz grch38_eigen_phred_coding.vcf

bgzip grch38_eigen_phred_coding.vcf
tabix grch38_eigen_phred_coding.vcf.gz

rm Eigen_hg38_coding_annot_04092016.tab.gz
```


### FATHMM XF
GRCh38:
<br>
https://fathmm.biocompute.org.uk/fathmm-xf/fathmm_xf_coding_hg38.vcf.gz

```
wget -c https://fathmm.biocompute.org.uk/fathmm-xf/fathmm_xf_coding_hg38.vcf.gz

python3 prepare_FathmmXF_vcf.py fathmm_xf_coding_hg38.vcf.gz grch38_fathmm_xf_coding.vcf

bgzip grch38_fathmm_xf_coding.vcf
tabix -p vcf grch38_fathmm_xf_coding.vcf.gz

rm fathmm_xf_coding_hg38.vcf.gz
```


### GnomAD (oe_lof and homAF)
GRCh38:
\
https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz

Needs manual liftover!!!
https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz

```
wget -c https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz

zcat gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz | sed -n '1d;p' | awk -v OFS="\t" -F "\t" '{ print $75,$76,$77,$24}' > gnomad_OE.bed
cat gnomad_OE.bed | sort -k1,1 -k2,2n -k3,3n -t '	' | sed '/NA/s/\bNA//g' > gnomAD_OE_sorted.bed

CrossMap bed hg19ToHg38.over.chain.gz gnomAD_OE_sorted.bed gnomAD_OE_grch38.bed

cat gnomad_OE_grch38.bed | sort -k1,1 -k2,2n -k3,3n -t '	' | sed '/NA/s/\bNA//g' | awk 'NF==4' > gnomAD_OE_grch38_sorted.bed

rm gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz
rm gnomAD_OE_sorted.bed
```

### \[optional\] HGMD (needs a license)
The possibility to use HGMD in aiDIVA is optional due to the fact that you need a license for it. If you choose to include HGMD in the aiDIVA analysis you can use the public/professional version of HGMD.


### Low Confidence Regions
GRCh38:
<br>
https://github.com/imgag/megSAP/raw/master/data/misc/low_conf_regions.bed

```
wget -c -O grch38_low_conf_region.bed https://github.com/imgag/megSAP/raw/master/data/misc/low_conf_regions.bed
```


### MutationAssessor
GRCh38:
<br>
--- NOTE ---
Unfortunately the official MutationAssessor (r3) resource that we utilized for training the aiDIVA-RF model is not available anymore.
We trained a new model with the same settings as before excluding the MutationAssessor feature. This new model can be found togehter with the previous models and in the release tab.
If you use this model without MutationAssessor you need to also remove the feature from the feature list in the configuration file.

As an alternative source for MutationAssessor you can use the dbNFSP database (the instructions shown here are for the academic version of the database). We show here the preparation for the latest academic legacy version (v4.9a) found [here](https://sites.google.com/site/jpopgen/dbNSFP). 
You can also use the newer versions from [here](https://www.dbnsfp.org/) the newer versions (>= 5.x) does not need any further preparations, but you have to download the large BGZF file that contains all chromosomes in one file.
To access the newer versions it is mandatory to create an account!

```
wget -c -O dbNSFP4.9a.zip https://usf.box.com/shared/static/0tq7q3b8ucaxxkmfyvnb0ss7g58ptgcl

unzip dbNSFP4.9a.zip
zcat dbNSFP4.9a_variant.chr1.gz | head -n1 > h

zgrep -h -v ^#chr dbNSFP4.9a_variant.chr* | sort -k1,1 -k2,2n - | cat h - | bgzip -c > dbNSFP4.9a_grch38.gz
tabix -s 1 -b 2 -e 2 dbNSFP4.9a_grch38.gz
```

<!-- 
GRCh38:
<br>
https://web.archive.org/web/20160802103335/http://mutationassessor.org/r3/MA_scores_rel3_hg19_full.tar.gz

Needs manual liftover!!!
https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz
```
wget -c https://web.archive.org/web/20160802103335/http://mutationassessor.org/r3/MA_scores_rel3_hg19_full.tar.gz
tar -xzf MA_scores_rel3_hg19_full.tar.gz

python3 prepare_MutationAssessor_vcf.py $(ls -m MA_scores_rel3_hg19_full/*.csv | tr -d '[:space:]') grch37_precomputed_MutationAssessor_unsort.vcf


# Make sure that the created VCF is sorted
cat grch37_precomputed_MutationAssessor_unsort.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > grch37_precomputed_MutationAssessor.vcf

bgzip grch37_precomputed_MutationAssessor.vcf
tabix -p vcf grch37_precomputed_MutationAssessor.vcf.gz

# Make sure that the created VCF is sorted
cat grch38_precomputed_MutationAssessor_unsort.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > grch37_precomputed_MutationAssessor.vcf

bgzip grch38_precomputed_MutationAssessor.vcf
tabix -p vcf grch38_precomputed_MutationAssessor.vcf.gz

rm MA_scores_rel3_hg19_full.tar.gz
rm -r MA_scores_rel3_hg19_full/
rm grch37_precomputed_MutationAssessor_unsort.vcf
rm grch37_precomputed_MutationAssessor.vcf.gz
```
-->


### PhastCons
For PhastCons is no further preparation necessary. Just make sure that it is correctly specified in the YAML configuration file.

GRCH38:
<br>
phastCons_primate: https://hgdownload.soe.ucsc.edu/goldenPath/hg38/phastCons17way/hg38.phastCons17way.bw
<br>
phastCons_mammal: https://hgdownload.soe.ucsc.edu/goldenPath/hg38/phastCons30way/hg38.phastCons30way.bw
<br>
phastCons_vertebrate: https://hgdownload.soe.ucsc.edu/goldenPath/hg38/phastCons100way/hg38.phastCons100way.bw

### PhyloP
For PhyloP is no further preparation necessary. Just make sure that it is correctly specified in the YAML configuration file.

GRCh38:
<br>
phyloP_primate: https://hgdownload.soe.ucsc.edu/goldenPath/hg38/phyloP17way/hg38.phyloP17way.bw
<br>
phyloP_mammal: https://hgdownload.soe.ucsc.edu/goldenPath/hg38/phyloP30way/hg38.phyloP30way.bw
<br>
phyloP_vertebrate: https://hgdownload.soe.ucsc.edu/goldenPath/hg38/phyloP100way/hg38.phyloP100way.bw


### Segment Duplication
GRCh38:
<br>
https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/genomicSuperDups.txt.gz

```
wget -c -O hg38.genomicSuperDups.txt.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/genomicSuperDups.txt.gz

zcat hg38.genomicSuperDups.txt.gz | cut -f2,3,4,27 > grch38_segmentDuplication_unsort.bed

grep -v '#' grch38_segmentDuplication_unsort.bed | sort -k1,1 -k2,2n -k3,3n -t $'\t' > grch38_segmentDuplication.bed

rm hg38.genomicSuperDups.txt.gz
rm grch38_segmentDuplication_unsort.bed
```


### SimpleRepeat
GRCh38:
<br>
https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/simpleRepeat.txt.gz

```
wget -c -O hg38.simpleRepeat.txt.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/simpleRepeat.txt.gz

zcat hg38.simpleRepeat.txt.gz | cut -f2,3,4,11 > grch38_simpleRepeat_unsort.bed

grep -v '#' grch38_simpleRepeat_unsort.bed | sort -k1,1 -k2,2n -k3,3n -t '	' > grch38_simpleRepeat.bed

rm hg38.simpleRepeat.txt.gz
rm grch38_simpleRepeat_unsort.bed
```


### RepeatMasker
GRCh38:
<br>
https://www.repeatmasker.org/genomes/hg38/rmsk4.0.5_rb20140131/hg38.fa.out.gz

```
wget -c https://www.repeatmasker.org/genomes/hg38/rmsk4.0.5_rb20140131/hg38.fa.out.gz

python3 prepare_RepeatMasker_bed.py hg38.fa.out.gz grch38_repeatmasker.bed

rm hg38.fa.out.gz
```


### REVEL
GRCh38:
<br>
https://zenodo.org/records/7072866/files/revel-v1.3_all_chromosomes.zip

```
wget -c https://zenodo.org/records/7072866/files/revel-v1.3_all_chromosomes.zip

unzip revel-v1.3_all_chromosomes.zip

python3 prepare_REVEL_vcf.py revel_with_transcript_ids grch37_revel_v13_unsort.vcf grch38_revel_v13_unsort.vcf


# Make sure that the created VCF is sorted
cat grch38_revel_v13_unsort.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > grch38_revel_v13.vcf

bgzip grch38_revel_v13.vcf
tabix -p vcf grch38_revel_v13.vcf.gz

rm revel-v1.3_all_chromosomes.zip
rm revel_with_transcript_ids
rm grch38_revel_v13_unsort.vcf
```


### SpliceAI
GRCh38
<br>
!!! Account needed !!!
<br>
<br>
Manually download the following files `spliceai_scores.masked.snv.hg38.vcf.gz` and `spliceai_scores.masked.indel.hg38.vcf.gz`
<br>
[https://basespace.illumina.com/s/otSPW8hnhaZR](https://basespace.illumina.com/s/otSPW8hnhaZR)

```
tabix -p vcf spliceai_scores.masked.snv.hg38.vcf.gz
tabix -p vcf spliceai_scores.masked.indel.hg38.vcf.gz
```


### AlphaMissense
GRCh38
<br>
!!! Google Account needed !!!
<br>
<br>
Manually donwload the following files `AlphaMissense_hg38.tsv.gz`
<br>
https://console.cloud.google.com/storage/browser/dm_alphamissense

```
tabix -p vcf AlphaMissense_hg38.tsv.gz
```
