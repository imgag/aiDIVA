#!/usr/bin/env python
# downloaded from here http://compbio.charite.de/jenkins/job/hpo.annotations.monthly/lastStableBuild/artifact/annotation/ALL_SOURCES_TYPICAL_FEATURES_phenotype_to_genes.txt

# process the file to save a dictionary with
# HPO terms as keys
# genes as values.

## Goal to use this dictionary within eDiVA to mark all genes related with HPO terms specified

## TODO: Merge all the scripts for resource generation into one

import pickle
import argparse


parser = argparse.ArgumentParser(description='Generate eDiVA HPO - gene association database from HPO database download')
parser.add_argument('--db_file','-db',dest='db', type=str, help='HPO database downloaded from http://compbio.charite.de/jenkins/job/hpo.annotations.monthly/lastStableBuild/artifact/annotation/ALL_SOURCES_TYPICAL_FEATURES_phenotype_to_genes.txt')
args = parser.parse_args()

HPO_association = dict()
with open(args.db) as rd:
    for line in rd:
        if line.startswith('#'):
            pass
        else:
            ff = line.strip().split('\t')
            key = ff[0]
            value = ff[-1]
            try:
                HPO_association[key].append(value)
            except:
                HPO_association[key] = [value]

print('Output will be written to "HPO_gene_association.p". Ensure to put this file in the eDiVA Prioritize folder')
pickle.dump( HPO_association, open( "HPO_gene_assiciation.p", "wb" ) )
