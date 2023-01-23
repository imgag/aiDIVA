# Preparation Of Annotation Resources
This document provides links to download the annoation sources used in AIdiva. If it is necessary to further prepare these files it is explained and shown in the respective section with a coode snippet.
\
\
For some of the annotation sources there are no GRCh38/hg38 files available in these cases we lifted the GRCh37/hg19 sources manually to the new assembly using the tool [CrossMap](https://github.com/liguowang/CrossMap). For this manual liftover we also provide code snippets.


## Reference Genome Assembly
Note hg19 and GRCh37 differ in the mitochondrial DNA, but since mitchondrial variants are not supported by AIdiva it does not matter which of the two assemblies is used as reference.

GRCh37:
<br>
[http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz)
<br>
<br>
GRCh38:
<br>
[http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz](http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz)

```
gunzip hg19.fa.gz
bgzip hg19.fa

gunzip hg38.fa.gz
bgzip hg38.fa

samtools dict hg19.fa.gz
samtools dict hg38.fa.gz

samtools faidx hg19.fa.gz
samtools faidx hg38.fa.gz
```

## Annotation Databases
INFO: Please be adviced that some of the annotation sources AIdiva uses are only free for non-commercial and academic use!

If for some reason you choose to exclude one or more of the Annotation sources shown here you have to make sure to modify the source code accordingly. If an annotation source is used as a feature for the random forest and you exclude it in your analysis you have to train a new model before you can use Aidiva.


### CADD
GRCh37:
<br>
https://kircherlab.bihealth.org/download/CADD/v1.6/GRCh37/whole_genome_SNVs.tsv.gz

```
wget -c -O grch37_whole_genome_SNVs.tsv.gz https://kircherlab.bihealth.org/download/CADD/v1.6/GRCh37/whole_genome_SNVs.tsv.gz

python3 prepare_CADD_vcf.py grch37_whole_genome_SNVs.tsv.gz grch37_CADD_snvs_v16.vcf

bgzip grch37_CADD_snvs_v16.vcf

tabix -p vcf grch37_CADD_snvs_v16.vcf.gz

rm grch37_whole_genome_SNVs.tsv.gz
```

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
GRCh37:
<br>
https://zenodo.org/record/3928295/files/capice_v1.0_build37_snvs.tsv.gz

```
wget -c https://zenodo.org/record/3928295/files/capice_v1.0_build37_snvs.tsv.gz

python3 prepare_Capice_vcf.py capice_v1.0_build37_snvs.tsv.gz grch37_capice_v1_snvs.vcf

bgzip grch37_capice_v1_snvs.vcf

tabix -p vcf grch37_capice_v1_snvs.vcf.gz

rm capice_v1.0_build37_snvs.tsv.gz
```


GRCh38:
<br>
Needs manual liftover of the previously generated grch37_capice_v1_snvs.vcf.gz
https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz

```
CrossMap.py vcf hg19ToHg38.over.chain.gz grch37_capice_v1_snvs.vcf.gz hg38.fa.gz grch38_capice_v1_snvs.vcf

bgzip grch38_capice_v1_snvs.vcf

tabix -p vcf grch38_capice_v1_snvs.vcf.gz
```


### ClinVar
GRCh37:
<br>
https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz

```
wget -c https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz

python3 prepare_ClinVar_vcf.py clinvar.vcf.gz grc37_clinvar.vcf

bgzip grch37_clinvar.vcf

tabix -p vcf grch37_clinvar.vcf.gz

rm clinvar.vcf.gz
```

GRCh38:
<br>
https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz

```
wget -c https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz

python3 prepare_ClinVar_vcf.py clinvar.vcf.gz grc38_clinvar.vcf

bgzip grch38_clinvar.vcf

tabix -p vcf grch38_clinvar.vcf.gz

rm clinvar.vcf.gz
```


### Condel
Please first download the FannsDB database from here: [https://bbglab.irbbarcelona.org/fannsdb/](https://bbglab.irbbarcelona.org/fannsdb/) 
(a free registration is required).
\
\
Afterwards the `prepare_Condel_vcf.py` script found in the _scripts_ folder can be used to create the VCF annotation file. For the annotation the VCF file needs to be _bgzipped_ and _indexed_

The following code snippet shows the necessary commands

GRCh37:

```
python3 prepare_Condel_vcf.py fannsdb.tsv.gz grch37_fannsDB_Condel.vcf

bgzip grch37_fannsDB_Condel.vcf

tabix -p vcf grch37_fannsDB_Condel.vcf.gz

rm fannsdb.tsv.gz
```


GRCh38:
\
Needs manual liftover of the previously generated grch37_fannsDB_Condel.vcf.gz
https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz

```
CrossMap.py vcf hg19ToHg38.over.chain.gz grch37_fannsDB_Condel.vcf.gz hg38.fa.gz grch38_fannsDB_Condel.vcf

bgzip grch38_fannsDB_Condel.vcf

tabix -p vcf grch38_fannsDB_Condel.vcf.gz
```


### Eigen phred
GRCh37:
\
https://zenodo.org/record/4046258/files/Eigen_hg19_coding_annot_04092016.tab.bgz

```
wget -c https://zenodo.org/record/4046258/files/Eigen_hg19_coding_annot_04092016.tab.bgz

python3 prepare_EigenPhred_vcf.py Eigen_hg19_coding_annot_04092016.tab.bgz grch37_eigen_phred_coding.vcf

bgzip grch37_eigen_phred_coding.vcf

tabix grch37_eigen_phred_coding.vcf.gz

rm Eigen_hg19_coding_annot_04092016.tab.bgz
```


GRCh38:
\
https://zenodo.org/record/3256671/files/Eigen_hg38_coding_annot_04092016.tab.gz

```
wget -c https://zenodo.org/record/3256671/files/Eigen_hg38_coding_annot_04092016.tab.gz

python3 prepare_EigenPhred_vcf.py Eigen_hg38_coding_annot_04092016.tab.bgz grch38_eigen_phred_coding.vcf

bgzip grch38_eigen_phred_coding.vcf

tabix grch38_eigen_phred_coding.vcf.gz

rm Eigen_hg38_coding_annot_04092016.tab.bgz
```


### FATHMM XF
GRCh37:
<br>
https://fathmm.biocompute.org.uk/fathmm-xf/fathmm_xf_coding.vcf.gz

```
wget -c https://fathmm.biocompute.org.uk/fathmm-xf/fathmm_xf_coding.vcf.gz

python3 prepare_FathmmXF_vcf.py fathmm_xf_coding.vcf.gz grch37_fathmm_xf_coding.vcf

bgzip grch37_fathmm_xf_coding.vcf

tabix -p vcf grch37_fathmm_xf_coding.vcf.gz

rm fathmm_xf_coding.vcf.gz
```

GRCh38:
<br>
https://fathmm.biocompute.org.uk/fathmm-xf/fathmm_xf_coding_hg38.vcf.gz

```
wget -c https://fathmm.biocompute.org.uk/fathmm-xf/fathmm_xf_coding_hg38.vcf.gz

python3 prepare_Fathmm_XF_vcf.py fathmm_xf_coding_hg38.vcf.gz grch38_fathmm_xf_coding.vcf

bgzip grch38_fathmm_xf_coding.vcf

tabix -p vcf grch38_fathmm_xf_coding.vcf.gz

rm fathmm_xf_coding_hg38.vcf.gz
```


### GnomAD (oe_lof and homAF)
The oe_lof score is gene based, therefor this resource is the same for GRCh37/hg19 and GRCh38/hg38.
\
https://storage.googleapis.com/gnomad-public/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz

```
wget -c https://storage.googleapis.com/gnomad-public/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz

zcat gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz | sed -n '1d;p' | awk -v OFS="\t" -F "\t" '{ print $$75,$$76,$$77,$$24}' > gnomad_OE.bed

cat gnomad_OE.bed | sort -k1,1 -k2,2n -k3,3n -t '	' | sed '/NA/s/\bNA//g' > gnomAD_OE_sorted.bed

rm gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz
```


If disk storage is a problem it is also possible to use the exome dataset instead of the genome one.
\
\
GRCh37:
\
https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.sites.vcf.bgz
\
\
(https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.vcf.bgz)
\
\
From the complete gnomAD database we are only interested in a few entries from the VCF files, therefor we suggest to prepare the VCF and remove the unnecessary stuff from the gnomAD annotation VCF. The following code snippet uses the `prepare_gnomAD_vcf.py` script to reduce it to a minimal set of INFO entries that we need for the AIdiva annotations:

```
wget -c https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.vcf.bgz

python3 prepare_gnomAD_vcf.py gnomad.genomes.r2.1.1.sites.vcf.bgz grch37_gnomAD_genomes_r211.vcf

bgzip grch37_gnomAD_genomes_r211.vcf

tabix -p vcf grch37_gnomAD_genomes_r211.vcf.gz
```

GRCh38:
\
https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/liftover_grch38/vcf/genomes/gnomad.genomes.r2.1.1.sites.liftover_grch38.vcf.bgz
\
\
(https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/liftover_grch38/vcf/exomes/gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz)

```
wget -c https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/liftover_grch38/vcf/genomes/gnomad.genomes.r2.1.1.sites.liftover_grch38.vcf.bgz

python3 prepare_gnomAD_vcf.py gnomad.genomes.r2.1.1.sites.liftover_grch38.vcf.bgz grch38_gnomAD_genomes_r211.vcf

bgzip grch38_gnomAD_genomes_r211.vcf

tabix -p vcf grch38_gnomAD_genomes_r211.vcf.gz
```


### \[optional\] HGMD (needs license)
The possibility to use HGMD in AIdiva is optional due to the fact that you need a license for it. If you choose to include HGMD in the AIdiva analysis you can use the public/professional version of HGMD.


### \[optional\] OMIM (needs license)
The possibility to use OMIM in AIdiva is optional due to the fact that you need a license for it.
**Currently not implemented!**


### Low Confidence Regions
GRCh37:
<br>
https://github.com/imgag/megSAP/raw/GRCh37/data/misc/low_conf_regions.bed

```
wget -c -O grch37_low_conf_region.bed https://github.com/imgag/megSAP/raw/GRCh37/data/misc/low_conf_regions.bed
```

GRCh38:
<br>
https://github.com/imgag/megSAP/raw/master/data/misc/low_conf_regions.bed

```
wget -c -O grch38_low_conf_region.bed https://github.com/imgag/megSAP/raw/master/data/misc/low_conf_regions.bed
```


### MutationAssessor
GRCh37:
<br>
http://mutationassessor.org/r3/MA_scores_rel3_hg19_full.tar.gz

```
wget -c http://mutationassessor.org/r3/MA_scores_rel3_hg19_full.tar.gz

tar -xzf MA_scores_rel3_hg19_full.tar.gz

python3 prepare_MutationAssessor_vcf.py $(ls -m MA_scores_rel3_hg19_full/*.csv | tr -d '[:space:]') grch37_precomputed_MutationAssessor_unsort.vcf

# use VcfSort from ngs-bits to make sure that the created VCF is sorted
ngs-bits/VcfSort -in grch37_precomputed_MutationAssessor_unsort.vcf -out grch37_precomputed_MutationAssessor.vcf

bgzip grch37_precomputed_MutationAssessor.vcf

tabix -p vcf grch37_precomputed_MutationAssessor.vcf.gz

rm MA_scores_rel3_hg19_full.tar.gz
rm -r MA_scores_rel3_hg19_full/
rm grch37_precomputed_MutationAssessor_unsort.vcf
```

GRCh38:
<br>
Needs manual liftover
https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz

```
CrossMap.py vcf hg19ToHg38.over.chain.gz grch37_precomputed_MutationAssessor.vcf.gz hg38.fa.gz grch38_precomputed_MutationAssessor.vcf

bgzip grch38_precomputed_MutationAssessor.vcf

tabix -p vcf grch38_precomputed_MutationAssessor.vcf.gz
```


### PhastCons
For PhastCons is no further preparation necessary. Just make sure that it is correctly specified in the YAML configuration file.

GRCH37:
\
phastCons_primate: https://hgdownload.soe.ucsc.edu/goldenPath/hg19/phastCons46way/primates.phastCons46way.bw

phastCons_mammal: https://hgdownload.soe.ucsc.edu/goldenPath/hg19/phastCons46way/placentalMammals.phastCons46way.bw

phastCons_vertebrate: https://hgdownload.soe.ucsc.edu/goldenPath/hg19/phastCons46way/vertebrate.phastCons46way.bw
\
\
GRCH38:
\
phastCons_primate: https://hgdownload.soe.ucsc.edu/goldenPath/hg38/phastCons17way/hg38.phastCons17way.bw

phastCons_mammal: https://hgdownload.soe.ucsc.edu/goldenPath/hg38/phastCons30way/hg38.phastCons30way.bw

phastCons_vertebrate: https://hgdownload.soe.ucsc.edu/goldenPath/hg38/phastCons100way/hg38.phastCons100way.bw

### PhyloP
For PhyloP is no further preparation necessary. Just make sure that it is correctly specified in the YAML configuration file.
\
\
GRCh37:
\
phyloP_primate: https://hgdownload.soe.ucsc.edu/goldenPath/hg19/phyloP46way/primates.phyloP46way.bw

phyloP_mammal: https://hgdownload.soe.ucsc.edu/goldenPath/hg19/phyloP46way/placentalMammals.phyloP46way.bw

phyloP_vertebrate: https://hgdownload.soe.ucsc.edu/goldenPath/hg19/phyloP46way/vertebrate.phyloP46way.bw

GRCh38:
\
phyloP_primate: https://hgdownload.soe.ucsc.edu/goldenPath/hg38/phyloP17way/hg38.phyloP17way.bw

phyloP_mammal: https://hgdownload.soe.ucsc.edu/goldenPath/hg38/phyloP30way/hg38.phyloP30way.bw

phyloP_vertebrate: https://hgdownload.soe.ucsc.edu/goldenPath/hg38/phyloP100way/hg38.phyloP100way.bw


### Segment Duplication
GRCh37:
\
https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/genomicSuperDups.txt.gz

```
wget -c -O hg19.genomicSuperDups.txt.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/genomicSuperDups.txt.gz

zcat hg19.genomicSuperDups.txt.gz | cut -f2,3,4,27 > grch37_genomicSuperDups_unsort.bed

grep -v '#' grch37_genomicSuperDups_unsort.bed | sort -k1,1 -k2,2n -k3,3n -t $'\t' > grch37_genomicSuperDups.bed

rm genomicSuperDups.txt.gz
rm grch37_genomicSuperDups_unsort.bed
```


GRCh38:
\
https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/genomicSuperDups.txt.gz

```
wget -c -O hg38.genomicSuperDups.txt.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/genomicSuperDups.txt.gz

zcat hg38.genomicSuperDups.txt.gz | cut -f2,3,4,27 > grch38_segmentDuplication_unsort.bed

grep -v '#' grch38_segmentDuplication_unsort.bed | sort -k1,1 -k2,2n -k3,3n -t $'\t' > grch38_segmentDuplication.bed

rm hg38.genomicSuperDups.txt.gz
rm grch38_segmentDuplication_unsort.bed
```


### SimpleRepeat
GRCh37
\
https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/simpleRepeat.txt.gz

```
wget -c -O hg19.simpleRepeat.txt.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/simpleRepeat.txt.gz

zcat hg19.simpleRepeat.txt.gz | cut -f2,3,4,11 > grch37_simpleRepeat_unsort.bed

grep -v '#' grch37_simpleRepeat_unsort.bed | sort -k1,1 -k2,2n -k3,3n -t '	' > grch37_simpleRepeat.bed

rm hg19.simpleRepeat.txt.gz
rm grch37_simpleRepeat_unsort.bed
```


GRCh38:
\
https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/simpleRepeat.txt.gz

```
wget -c -O hg38.simpleRepeat.txt.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/simpleRepeat.txt.gz

zcat hg19.simpleRepeat.txt.gz | cut -f2,3,4,11 > grch38_simpleRepeat_unsort.bed

grep -v '#' grch38_simpleRepeat_unsort.bed | sort -k1,1 -k2,2n -k3,3n -t '	' > grch38_simpleRepeat.bed

rm hg38.simpleRepeat.txt.gz
rm grch38_simpleRepeat_unsort.bed
```


### RepeatMasker
GRCh37:
\
http://www.repeatmasker.org/genomes/hg19/RepeatMasker-rm405-db20140131/hg19.fa.out.gz

```
wget -c http://www.repeatmasker.org/genomes/hg19/RepeatMasker-rm405-db20140131/hg19.fa.out.gz

python3 prepare_RepeatMasker_bed.py hg19.fa.out.gz grch37_repeatmasker.bed

rm hg19.fa.out.gz
```


GRCh38:
\
http://www.repeatmasker.org/genomes/hg38/RepeatMasker-rm405-db20140131/hg38.fa.out.gz

```
wget -c http://www.repeatmasker.org/genomes/hg38/RepeatMasker-rm405-db20140131/hg38.fa.out.gz

python3 prepare_RepeatMasker_bed.py hg38.fa.out.gz grch38_repeatmasker.bed

rm hg38.fa.out.gz
```


### REVEL
GRCh37 and GRCh38:
\
https://rothsj06.u.hpc.mssm.edu/revel-v1.3_all_chromosomes.zip

```
wget -c https://rothsj06.u.hpc.mssm.edu/revel-v1.3_all_chromosomes.zip

unzip revel-v1.3_all_chromosomes.zip

python3 prepare_REVEL_vcf.py revel_with_transcript_ids grch37_revel_v13_unsort.vcf grch38_revel_v13_unsort.vcf

ngs-bits/VcfSort -in grch37_revel_v13_unsort.vcf -out grch37_revel_v13.vcf
ngs-bits/VcfSort -in grch38_revel_v13_unsort.vcf -out grch38_revel_v13.vcf

bgzip grch37_revel_v13.vcf
bgzip grch38_revel_v13.vcf

tabix -p vcf grch37_revel_v13.vcf.gz
tabix -p vcf grch38_revel_v13.vcf.gz

rm revel-v1.3_all_chromosomes.zip
rm revel_with_transcript_ids
rm grch37_revel_v13_unsort.vcf
rm grch38_revel_v13_unsort.vcf
```


### dbscSNV
GRCh37 and GRCh38:
\
ftp://dbnsfp:dbnsfp@dbnsfp.softgenetics.com/dbscSNV1.1.zip

```
wget -c ftp://dbnsfp:dbnsfp@dbnsfp.softgenetics.com/dbscSNV1.1.zip

mkdir dbscSNV_files

unzip dbscSNV1.1.zip -d dbscSNV_files

python3 prepare_dbscSNV_vcf.py $(ls -m dbscSNV_files/dbscSNV1.1.* | tr -d '[:space:]') grch37_dbscSNV_scores_unsort.vcf grch38_dbscSNV_scores_unsort.vcf

# use VcfSort from ngs-bits to make sure that the created VCF is sorted
ngs-bits/VcfSort -in grch37_dbscSNV_scores_unsort.vcf -out grch37_dbscSNV_scores.vcf
ngs-bits/VcfSort -in grch38_dbscSNV_scores_unsort.vcf -out grch38_dbscSNV_scores.vcf

bgzip grch37_dbscSNV_scores.vcf
bgzip grch38_dbscSNV_scores.vcf

tabix -p vcf grch37_dbscSNV_scores_unsort.vcf.gz
tabix -p vcf grch38_dbscSNV_scores_unsort.vcf.gz

rm dbscSNV1.1.zip
rm -r dbscSNV_files/
rm grch37_dbscSNV_scores_unsort.vcf
rm grch38_dbscSNV_scores_unsort.vcf
```

