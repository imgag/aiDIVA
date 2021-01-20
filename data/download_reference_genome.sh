#!/bin/bash

folder=`pwd`
genomes_folder=$folder/genomes/
mkdir -p $genomes_folder

# Download GRCh37 reference genome
cd $genomes_folder
wget wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
zcat hs37d5.fa.gz | sed -r 's/>/>chr/g' > GRCh37.fa
rm hs37d5.fa.gz
$folder/tools/samtools-1.10/samtools faidx GRCh37.fa

# Download GRCh38 reference genome
cd $genomes_folder
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz
zcat GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz | sed -r 's/>chrM/>chrMT/g' > GRCh38.fa
rm GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz
$folder/tools/samtools-1.10/samtools faidx GRCh38.fa
