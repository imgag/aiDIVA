SHELL:=/bin/bash

help:
	@cat Makefile_HPO.mk

download:
	wget -c http://purl.obolibrary.org/obo/hp/hpoa/phenotype_to_genes.txt
	wget -c http://purl.obolibrary.org/obo/hp.obo
	wget -c http://purl.obolibrary.org/obo/hp/hpoa/phenotype_annotation.tab
	wget -c http://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/hgnc_complete_set.txt
	wget -c https://stringdb-static.org/download/protein.links.detailed.v11.0/9606.protein.links.detailed.v11.0.txt.gz
	wget -c https://string-db.org/mapping_files/STRING_display_names/human.name_2_string.tsv.gz
	wget -O grch38_ensembl_transcript_length_and_strand.tsv 'https://www.ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "1" count = "" datasetConfigVersion = "0.6" ><Dataset name = "hsapiens_gene_ensembl" interface = "default" ><Attribute name = "ensembl_gene_id" /><Attribute name = "ensembl_transcript_id" /><Attribute name = "cds_length" /><Attribute name = "transcript_length" /><Attribute name = "strand" /></Dataset></Query>'
	wget -O grch37_ensembl_transcript_length_and_strand.tsv 'https://grch37.ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "1" count = "" datasetConfigVersion = "0.6" ><Dataset name = "hsapiens_gene_ensembl" interface = "default" ><Attribute name = "ensembl_gene_id" /><Attribute name = "ensembl_transcript_id" /><Attribute name = "cds_length" /><Attribute name = "transcript_length" /><Attribute name = "strand" /></Dataset></Query>'
	
convert:
	awk -F '	'  '{print $$5}' < phenotype_annotation.tab | sort  | uniq -c | awk '{print $$2 "	" $$1}' > HPO_counts.txt
	python3 generate_HPO_resources.py --hpo_ontology hp.obo --gene_phenotype phenotype_to_genes.txt --gene_hpo gene2hpo.json --hpo_counts HPO_counts.txt --hpo_graph hpo_graph.gexf --hpo_replacements hpo2replacement.json --hgnc_symbols hgnc_complete_set.txt --hgnc_gene hgnc2gene.json --string_links 9606.protein.links.detailed.v11.0.txt.gz --string_mapping human.name_2_string.tsv.gz --gene_interacting gene2interacting.json --transcript_length_mapping transcript2length.json --transcript_information grch37_ensembl_transcript_length_and_strand.tsv
	
	rm phenotype_to_genes.txt
	rm hp.obo
	rm phenotype_annotation.tab
	rm hgnc_complete_set.txt
	rm 9606.protein.links.detailed.v11.0.txt.gz
	rm human.name_2_string.tsv.gz
	rm HPO_counts.txt
	rm grch37_ensembl_transcript_length_and_strand.tsv
	rm grch38_ensembl_transcript_length_and_strand.tsv
	
	mv gene2hpo.json ../hpo_resources/
	mv gene2interacting.json ../hpo_resources/
	mv hgnc2gene.json ../hpo_resources/
	mv hpo_graph.gexf ../hpo_resources/
	mv hpo2replacement.json ../hpo_resources/
	mv transcript2length.json ../hpo_resources/