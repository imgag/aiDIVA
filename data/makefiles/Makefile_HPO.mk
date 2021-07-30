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

convert:
	awk -F '	'  '{print $$5}' < phenotype_annotation.tab | sort  | uniq -c | awk '{print $$2 "	" $$1}' > HPO_counts.txt
	python3 generate_HPO_resources.py --hpo_ontology hp.obo --gene_phenotype phenotype_to_genes.txt --gene_hpo gene2hpo.pkl --hpo_edges hpo_edges.pkl --hpo_counts HPO_counts.txt --hpo_graph hpo_graph.pkl --hgnc_symbols hgnc_complete_set.txt --hgnc_gene hgnc2gene.pkl --string_links 9606.protein.links.detailed.v11.0.txt.gz --string_mapping human.name_2_string.tsv.gz --gene_interacting gene2interacting.pkl
	
	rm phenotype_to_genes.txt
	rm hp.obo
	rm phenotype_annotation.tab
	rm hgnc_complete_set.txt
	rm 9606.protein.links.detailed.v11.0.txt.gz
	rm human.name_2_string.tsv.gz
	rm hpo_edges.pkl
	rm HPO_counts.txt
	
	mv gene2hpo.pkl ../hpo_resources/
	mv gene2interacting.pkl ../hpo_resources/
	mv hgnc2gene.pkl ../hpo_resources/
	mv hpo_graph.pkl ../hpo_resources/