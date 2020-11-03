#!/bin/bash

folder=`pwd`
genomes_folder=$folder/genomes/

cd $genomes_folder
wget wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
gunzip hs37d5.fa.gz |  sed -r 's/>/>chr/g' > GRCh37.fa
rm hs37d5.fa.gz
$folder/tools/samtools-1.10/samtools faidx GRCh37.fa

## TODO: add GRCh38
