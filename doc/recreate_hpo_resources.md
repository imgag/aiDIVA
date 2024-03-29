# Recreation Of The Included HPO Resources
It is recommended to regularly recreate the HPO resources that are included in the repository (data/hpo_resources/). To include the newest changes from the underlying databases.

## Preparation And Download Of Needed Resources
HPO graph:
```
wget -c wget http://purl.obolibrary.org/obo/hp.obo
wget -c http://purl.obolibrary.org/obo/hp/hpoa/phenotype_annotation.tab

awk -F '\t'  '{print $5}' < phenotype_annotation.tab | sort  | uniq -c | awk '{print $2 "\t" $1}' > HPO_counts.txt
```

Phenotype to gene mapping:
```
wget -c http://purl.obolibrary.org/obo/hp/hpoa/phenotype_to_genes.txt
```

HGNC ID mapping:
```
wget -c ftp://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/hgnc_complete_set.txt
```

Protein interaction mapping:
```
wget -c https://stringdb-static.org/download/protein.links.detailed.v11.0/9606.protein.links.detailed.v11.0.txt.gz
wget -c https://string-db.org/mapping_files/STRING_display_names/human.name_2_string.tsv.gz
```

Transcript length mapping:
```
wget -c -O grch38_ensembl_transcript_length_and_strand.tsv 'https://www.ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "1" count = "" datasetConfigVersion = "0.6" ><Dataset name = "hsapiens_gene_ensembl" interface = "default" ><Attribute name = "ensembl_gene_id" /><Attribute name = "ensembl_transcript_id" /><Attribute name = "cds_length" /><Attribute name = "transcript_length" /><Attribute name = "strand" /></Dataset></Query>'
wget -c -O grch37_ensembl_transcript_length_and_strand.tsv 'https://grch37.ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "1" count = "" datasetConfigVersion = "0.6" ><Dataset name = "hsapiens_gene_ensembl" interface = "default" ><Attribute name = "ensembl_gene_id" /><Attribute name = "ensembl_transcript_id" /><Attribute name = "cds_length" /><Attribute name = "transcript_length" /><Attribute name = "strand" /></Dataset></Query>'
```

## Creation Of Resources
After you successfully downlaoded and prepared the above mentioned resources you can sue the `prepare_HPO_resources.py` script in the _scripts_ folder to create the needed HPO resources for AIdiva. The following code snippets show this process in more detail.

__hpo_graph.gexf__ and __hpo2replacement.json__
```
python3 prepare_HPO_resources.py --hpo_ontology hp.obo --hpo_counts HPO_counts.txt --hpo_graph hpo_graph.gexf --hpo_replacements hpo2replacement.json
```

__gene2hpo.json__
```
python3 prepare_HPO_resources.py --gene_phenotype phenotype_to_genes.txt --gene_hpo gene2hpo.json
```

__hgnc2gene.json__
```
python3 prepare_HPO_resources.py --hgnc_symbols hgnc_complete_set.txt --hgnc_gene hgnc2gene.json
```

__gene2interacting.json__
```
python3 prepare_HPO_resources.py --string_mapping human.name_2_string.tsv.gz --string_links 9606.protein.links.detailed.v11.0.txt.gz --gene_interacting gene2interacting.json
```

__grch37transcript2length.json__ and __grch38transcript2length.json__
```
python3 prepare_HPO_resources.py 
```

Note: You can also use all of the above parameters in one single command and create all resources at once with a single command call. Here they are splitted in separate call for better readability.