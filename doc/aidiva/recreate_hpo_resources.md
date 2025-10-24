# Recreation of the Included HPO Resources
It is recommended to regularly recreate the HPO resources that are included in the repository (data/hpo_resources/). To include the newest changes from the underlying databases.

## Preparation and Download of Needed Resources
HPO graph:
```
wget -c http://purl.obolibrary.org/obo/hp.obo
wget -c http://purl.obolibrary.org/obo/hp/hpoa/phenotype.hpoa

awk -F '\t'  '{print $4}' < phenotype.hpoa | grep -v '^#' | sort  | uniq -c | awk '{print $2 "\t" $1}' > HPO_counts.txt
```

Phenotype to gene mapping:
```
wget -c http://purl.obolibrary.org/obo/hp/hpoa/phenotype_to_genes.txt
```

HGNC ID mapping:
```
wget -c https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt
```

Protein interaction mapping (v11.0b):
```
wget -c https://stringdb-static.org/download/protein.links.detailed.v11.0/9606.protein.links.detailed.v11.0.txt.gz
wget -c https://version-11-0b.string-db.org/mapping_files/STRING_display_names/human.name_2_string.tsv.gz
```

Transcript length mapping:
```
wget -c -O grch38_ensembl_transcript_length_and_strand.tsv 'https://www.ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "1" count = "" datasetConfigVersion = "0.6" ><Dataset name = "hsapiens_gene_ensembl" interface = "default" ><Attribute name = "ensembl_gene_id" /><Attribute name = "ensembl_transcript_id" /><Attribute name = "cds_length" /><Attribute name = "transcript_length" /><Attribute name = "strand" /></Dataset></Query>'
wget -c -O grch37_ensembl_transcript_length_and_strand.tsv 'https://grch37.ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "1" count = "" datasetConfigVersion = "0.6" ><Dataset name = "hsapiens_gene_ensembl" interface = "default" ><Attribute name = "ensembl_gene_id" /><Attribute name = "ensembl_transcript_id" /><Attribute name = "cds_length" /><Attribute name = "transcript_length" /><Attribute name = "strand" /></Dataset></Query>'
```

## Creation of Resources
After you successfully downloaded and prepared the above mentioned resources you can use the `prepare_HPO_resources.py` script in the _scripts_ folder to create the needed HPO resources for aiDIVA. The following code snippets show this process in more detail.

**hpo_graph.gexf** and **hpo2replacement.json**
```
python3 prepare_HPO_resources.py --hpo_ontology hp.obo --hpo_counts HPO_counts.txt --hpo_graph hpo_graph.gexf --hpo_replacements hpo2replacement.json
```

**gene2hpo.json**
```
python3 prepare_HPO_resources.py --gene_phenotype phenotype_to_genes.txt --gene_hpo gene2hpo.json
```

**hpo2name.json**
```
python3 prepare_HPO_resources.py --hpo_ontology hp.obo --hpo_name hpo2name.json
```

**hgnc2gene.json**
```
python3 prepare_HPO_resources.py --hgnc_symbols hgnc_complete_set.txt --hgnc_gene hgnc2gene.json
```

**gene2interacting.json**
```
python3 prepare_HPO_resources.py --string_mapping human.name_2_string.tsv.gz --string_links 9606.protein.links.detailed.v11.0.txt.gz --gene_interacting gene2interacting.json
```

**grch37transcript2length.json** and **grch38transcript2length.json**
```
python3 prepare_HPO_resources.py --transcript_information grch37_ensembl_transcript_length_and_strand.tsv --transcript_mapping grch37transcript2length.json
python3 prepare_HPO_resources.py --transcript_information grch38_ensembl_transcript_length_and_strand.tsv --transcript_mapping grch38transcript2length.json
```

Note: You can also use all of the above parameters in one single command and create all resources at once with a single command call. Here they are splitted in separate calls for better readability.