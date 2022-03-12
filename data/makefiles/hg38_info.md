If you want to use the GRCh38 genome assembly you need to use the hg38 annotation sources.
Some of them are already present and you only need to download the correct file. The following list shows the features that can be directly downloaded:
	- simpleRepeat
	- segmentDuplication
	- phastCons
	- phyloP
	- FATHMM_XF


Unfortunately for the following features there are no new resources for hg38 available so these have to be lifted over from the hg19 resource files.
 	- ABB-SCORE (score is not used any more)
	- Condel
	- Eigen
	- MutationAssessor


The tool we used to accomplish this is called CrossMap (http://crossmap.sourceforge.net).
CrossMap needs a chain file and a reference fasta file in case you want to convert VCF files:
	- chain file (ftp://ftp.ensembl.org/pub/assembly_mapping/homo_sapiens/GRCh37_to_GRCh38.chain.gz)
	- target reference fasta (ftp://ftp.ensembl.org/pub/release-76/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz)
