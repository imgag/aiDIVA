# Recreation of the Included HPO Resources
Before you can successfully run aiDIVA you need to download some resources and make them available to the tool by either specifying the path in the configuration file or putting the resource file in the dedicated folder in the installation directory (`<aidiva-installation-folder>/data/hpo_resources`). It is recommended to regularly update the resource files. To include the newest changes from the underlying databases. 

## Preparation and Download of Needed Resources
HPO graph:
<br>
THe given links should point to the latest version of the resources
```
wget -c http://purl.obolibrary.org/obo/hp.obo
wget -c http://purl.obolibrary.org/obo/hp/hpoa/phenotype.hpoa
```

Phenotype to gene mapping:
```
wget -c http://purl.obolibrary.org/obo/hp/hpoa/phenotype_to_genes.txt
```

HGNC ID mapping:
```
wget -c https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt
```

Protein interaction mapping (v12.0):
<br>
Version 12.0 is the newest version of the String-DB at the time of writing this documentation
```
wget https://stringdb-downloads.org/download/protein.links.detailed.v12.0/9606.protein.links.detailed.v12.0.txt.gz
wget https://stringdb-downloads.org/download/protein.aliases.v12.0/9606.protein.aliases.v12.0.txt.gz
```

Transcript length mapping:
```
wget -c -O grch38_ensembl_transcript_length_and_strand.tsv 'https://www.ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "1" count = "" datasetConfigVersion = "0.6" ><Dataset name = "hsapiens_gene_ensembl" interface = "default" ><Attribute name = "ensembl_gene_id" /><Attribute name = "ensembl_transcript_id" /><Attribute name = "cds_length" /><Attribute name = "transcript_length" /><Attribute name = "strand" /></Dataset></Query>'
wget -c -O grch37_ensembl_transcript_length_and_strand.tsv 'https://grch37.ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "1" count = "" datasetConfigVersion = "0.6" ><Dataset name = "hsapiens_gene_ensembl" interface = "default" ><Attribute name = "ensembl_gene_id" /><Attribute name = "ensembl_transcript_id" /><Attribute name = "cds_length" /><Attribute name = "transcript_length" /><Attribute name = "strand" /></Dataset></Query>'
```

## Creation of Resources
After you successfully downloaded the above mentioned resources you place the files into the folder `data/hpo_resources` so that aiDIVA can find the files. Alternatively you can change the path to the files in the config file by changing the path of the following keys in the Internal-Parameter section of the config file (see the `example_configuration_aiDIVA.yaml`) for aiDIVA: 

- hpo-ontology: \<your-full-path\>/hp.obo
- phenotype-information: \<your-full-path\>/phenotype.hpoa
- phenotype-to-genes: \<your-full-path\>/phenotype_to_genes.txt
- transcript-information: \<your-full-path\>/grch38_ensembl_transcript_length_and_strand.tsv
- hgnc-infromation: \<your-full-path\>/hgnc_complete_set.txt
- string-db-information: \<your-full-path\>/9606.protein.links.detailed.v12.0.txt.gz
- string-db-aliases: \<your-full-path\>/9606.protein.aliases.v12.0.txt.gz