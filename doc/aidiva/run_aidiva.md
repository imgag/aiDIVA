# Running aiDIVA

aiDIVA expects a TAB separted table as input (see [here](https://github.com/imgag/aiDIVA/blob/master/doc/aidiva/run_aidiva.md#overview-of-necessary-and-optional-columns-in-the-input-table) for an overview of the necessary and optional columns) acompanied by a configuration file in the YAML format. Have a look in the `data` folder to get an example configuration file with placeholders. 

It is possible to run only parts of our software. In total we have three different modi in which our software can be run: *aiDIVA-RF*, *aiDIVA-meta*, *aiDIVA-meta-RF*.

Please be adviced that for the *aiDIVA-meta* mode you need the two evidence-based files *eb_dom.GSvar* and *eb_rec.GSvar*. These files can be created using the [VariantRanking](https://github.com/imgag/ngs-bits/blob/master/doc/tools/VariantRanking/index.md) tool which is part of the [ngs-bits](https://github.com/imgag/ngs-bits) tool collection. The VariantRanking tool needs as input a GSvar file which can be created if you process your VCF file with the [megSAP](https://github.com/imgag/megSAP) pipeline. Please head over to these repositories to have detailled instructions on how to use the tools.


## Modes for Running aiDIVA

### aiDIVA-RF:

```
python3 run_aidiva-rf.py --config configuration_aiDIVA.yaml --in_data input.tsv --out_prefix output_folder/aidiva_result [--workdir aidiva_workdir/] [--hpo_list HP:xxxx,HP:xxxx] [--gene_exclusion gene_exclusion.txt] [--family_file family.txt] [--family_type SINGLE] [--skip_db_check] [--rare_disease] [--only_top_results] [--top_rank 25] [--threads 1] [--log_file output_path/logs/aidiva_log.txt] [--log_level INFO]
```
mandatory parameters:

+ *config* -- YAML configuration file (in the `data` folder is an example configuration file)
+ *in\_data* -- TAB separated input table with the annotated variants
+ *out\_prefix* -- A prefix for the resulting output files the output folder can also be specified with that parameter

optional parameters:

+ *workdir* -- Working directory, where all temporary files are created and saved \[optional\]
+ *hpo\_list* -- Comma-separated list of HPO terms observed with the patient \[optional\]
+ *gene\_exclusion* -- TXT file containing genes that should be excluded during the analysis of the HPO relatedness \[optional\]
+ *family_file* -- TXT file containing the sample information if run on multisample VCF files \[optional\]
+ *family_type* -- Type of the family relation \[SINGLE, TRIO\] (default: SINGLE) \[optional\]
+ *skip\_db\_check* -- Skip the database checkup for existing entries in ClinVar (and HGMD) \[optional\]
+ *rare\_disease* -- Adds initial allele frequence filter to only keep variants with allele frequence of less than 2%\[optional\]
+ *only\_top\_results* -- Restrict the results to only report the top X variants (default: 25) \[optional\]
+ *top\_rank* -- Rank that should be used as maximum for the top ranking results to report (default: 25) \[optional\]
+ *threads* -- Number of threads that should be used (default: 1) \[optional\]
+ *log\_file* -- Specify a custom log file to store the log messages from the tool \[optional\]
+ *log\_level* -- Define logging level \[DEBUG, INFO\] (default: INFO) \[optional\]


### aiDIVA-meta-RF:

```
python3 run_aidiva-meta-rf.py --config configuration_aiDIVA.yaml --in_data input.tsv --out_prefix output_folder/aidiva_result --sample_id NA12878 [--workdir aidiva_workdir/] [--hpo_list HP:xxxx,HP:xxxx] [--gene_exclusion gene_exclusion.txt] [--family_file family.txt] [--family_type SINGLE] [--skip_db_check] [--rare_disease] [--only_top_results] [--top_rank 25] [--gender male/female] [--age 42] [--evidence_based] [--threads 1] [--log_file output_path/logs/aidiva_log.txt] [--log_level INFO]
```
mandatory parameters:

+ *config* -- YAML configuration file (in the `data` folder is an example configuration file)
+ *in\_data* -- TAB separated input table with the annotated variants
+ *out\_prefix* -- A prefix for the resulting output files the output folder can also be specified with that parameter
+ *sample\_id* -- Sample ID this is used to extract the genotype in the VCF file

optional parameters:

+ *workdir* -- Working directory, where all temporary files are created and saved \[optional\]
+ *hpo\_list* -- Comma-separated list of HPO terms observed with the patient
+ *gene\_exclusion* -- TXT file containing genes that should be excluded during the analysis of the HPO relatedness \[optional\]
+ *family\_file* -- TXT file containing the sample information if run on multisample VCF files \[optional\]
+ *family\_type* -- Type of the family relation \[SINGLE, TRIO\] (default: SINGLE) \[optional\]
+ *skip\_db\_check* -- Skip the database checkup for existing entries in ClinVar (and HGMD) \[optional\]
+ *rare\_disease* -- Adds initial allele frequence filter to only keep variants with allele frequence of less than 2%\[optional\]
+ *only\_top\_results* -- Restrict the results to only report the top X variants (default: 25) \[optional\]
+ *top\_rank* -- Rank that should be used as maximum for the top ranking results to report (default: 25) \[optional\]
+ *gender* -- Gender of the patient if known \[optional\]
+ *age* -- Age of the patient if known \[optional\]
+ *threads* -- Number of threads that should be used (default: 1) \[optional\]
+ *log\_file* -- Specify a custom log file to store the log messages from the tool \[optional\]
+ *log\_level* -- Define logging level \[DEBUG, INFO\] (default: INFO) \[optional\]


### aiDIVA-meta:

```
python3 run_aidiva-meta.py --config configuration_aiDIVA.yaml --in_data input.tsv --in_eb_dom in_eb_dom.GSvar --in_eb_rec in_eb_rec.GSvar --out_prefix output_path/aidiva_result --sample_id NA12878 [--workdir aidiva_workdir/] [--hpo_list HP:xxxx,HP:xxxx] [--gene_exclusion gene_exclusion.txt] [--family_file family.txt] [--family_type SINGLE] [--skip_db_check] [--rare_disease] [--only_top_results] [--top_rank 25] [--gender male/female] [--age 42] [--evidence_based] [--threads 1] [--log_file output_path/logs/aidiva_log.txt] [--log_level INFO]
```
mandatory parameters:

+ *config* -- YAML configuration file (in the `data` folder is an example configuration file)
+ *in\_data* -- TAB separated input table with the annotated variants
+ *in\_eb\_dom* -- GSvar file containing the evidence-based ranking results of the sample using the dominant mode of the algorithm
+ *in\_eb\_rec* -- GSvar file containing the evidence-based ranking results of the sample using the recessive mode of the algorithm
+ *out\_prefix* -- A prefix for the resulting output files the output folder can also be specified with that parameter
+ *sample\_id* -- Sample ID this is used to access the genotype in the input table and the GSvar files

optional parameters:

+ *workdir* -- Working directory, where all temporary files are created and saved \[optional\]
+ *hpo\_list* -- Comma-separated list of HPO terms observed with the patient \[optional\]
+ *gene\_exclusion* -- TXT file containing genes that should be excluded during the analysis of the HPO relatedness \[optional\]
+ *family\_file* -- TXT file containing the sample information if run on multisample VCF files \[optional\]
+ *family\_type* -- Type of the family relation \[SINGLE, TRIO\] (default: SINGLE) \[optional\]
+ *skip\_db\_check* -- Skip the database checkup for existing entries in ClinVar (and HGMD) \[optional\]
+ *rare\_disease* -- Adds initial allele frequence filter to only keep variants with allele frequence of less than 2%\[optional\]
+ *only\_top\_results* -- Restrict the results to only report the top X variants (default: 25) \[optional\]
+ *top\_rank* -- Rank that should be used as maximum for the top ranking results to report (default: 25) \[optional\]
+ *gender* -- Gender of the patient if known \[optional\]
+ *age* -- Age of the patient if known \[optional\]
+ *threads* -- Number of threads that should be used (default: 1) \[optional\]
+ *log\_file* -- Specify a custom log file to store the log messages from the tool \[optional\]
+ *log\_level* -- Define logging level \[DEBUG, INFO\] (default: INFO) \[optional\]


## aiDIVA Results

aiDIVA will produce multiple different output files. The following lists all possible result files. Depending on your chosen mode (*aiDIVA-RF*, *aiDIVA-meta*, *aiDIVA-meta-RF*) to run aiDIVA you will only get a subset of these result files.

+ *\<your-result-prefix\>\_aidiva\_result.tsv* -- The unfiltered result table (aiDIVA-RF).
+ *\<your-result-prefix\>\_aidiva\_result\_filtered.tsv* -- The filtered result table (aiDIVA-RF) this table is also used for the subsequent analysis steps.
+ *\<your-result-prefix\>\_aidiva\_random\_forest\_based\_llm\_result.tsv* -- The LLM results based on the random forest-based ranking.
+ *\<your-result-prefix\>\_aidiva\_evidence\_dominant\_based\_llm\_result.tsv* -- The LLM results based on the evidence-based ranking (dominant model).
+ *\<your-result-prefix\>\_aidiva\_evidence\_recessive\_based\_llm\_result.tsv* -- The LLM results based on the evidence-based ranking (recessive model).
+ *\<your-result-prefix\>\_aidiva\_metascore\_results.tsv* -- The final result table with the variant ranking based on the meta model (aiDIVA-meta).


## Overview of Necessary and Optional Columns in the Input Table

The following shows all columns that need to be present in the annotated input table for the software to work.

The column names should exactly match the column name specified in the table if not otherwise specified in the description.

### Strictly necessary columns

These columns give the basic information for each variant in the table.

+ *#CHROM* -- Chromosome identifier
+ *POS* -- Variant position 
+ *REF* -- Reference allele
+ *ALT* -- Observed alternative allele
+ *FILTER* -- Matches the FILTER column of a VCF file (checked terms: off-target, low_conf_region)
+ *SYMBOL* -- Approved gene symbol
+ *IMPACT* -- Impact of the variant according to ensembl (HIGH, MODERATE, LOW, MODIFIER)


### Necessary feature columns

These column names must match the feature-list specified in the configuration file. The following shows the column names for the feature-list specified in the example configuration given in the `data` folder.

+ *SIFT*
+ *PolyPhen*
+ *CADD_PHRED*
+ *REVEL*
+ *MAX_AF*
+ *EIGEN_PHRED*
+ *CONDEL*
+ *FATHMM_XF*
+ *MutationAssessor*
+ *phastCons_mammal*
+ *phastCons_primate*
+ *phastCons_vertebrate*
+ *phyloP_mammal*
+ *phyloP_primate*
+ *phyloP_vertebrate*
+ *oe_lof*
+ *homAF*
+ *CAPICE*
+ *ALPHA_MISSENSE_SCORE*
+ *HIGH_IMPACT*
+ *IS_INDEL*


### Necessary allele frequency columns

These columns are necessary if the MAX_AF column from the feature list is not present in the table. Not restricted to the columns shown below. column names must match the entries in the allele-frequency-list specified in the configuration file.

+ *gnomAD_AFR_AF*
+ *gnomAD_AMR_AF*
+ *gnomAD_EAS_AF*
+ *gnomAD_NFE_AF*
+ *gnomAD_SAS_AF*


### Necessary splicing features

Scores used to handle splicing variants.

+ *SpliceAI* -- Single SpliceAI score per variant (we use the maximum of the four scores)


### Necessary variant consequence information

Variant consequence terms annotated with ensembl VEP. The terms are shown in the overview on the ensembl [website](https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html). This site was also used as reference to determine the severity of the different terms.

+ *Consequence* -- Consequence terms annotated by VEP
+ *MOST_SEVERE_CONSEQUENCE* -- Single consequence term per variant (if overlapping terms were present we chose the most severe one)


### Necessary Sample Information
Information on the genotype of the variant. \<sample-id\> must match the one you pass with the *--sample\_id* parameter.

+ *GT_\<sample-id\>* -- Genotype information of the sample (\<sample-id\> must match the id of the sample you are currently analyzing)


### Optional columns

The following columns are not necessarily needed to run aiDIVA.

+ *CLINVAR_DETAILS* -- ClinVar classification of the variant (needed if the parameter --skip_db_check is not given)
+ *HGMD_CLASS* -- HGMD classification of the variant (needed if the parameter --skip_db_check is not given)

<br>

+ *segmentDuplication* -- Segment duplciation annotation if the variant lies in a segment duplicaiton region.
+ *simpleRepeat* -- Simple repeat annotation if the variant lies inside a tandem repeat region.

<br>

+ *REPEATMASKER* -- RepeatMasker annotation (is used as additional filter in the prioritization step of aiDIVA-RF)

<br>


+ *HGNC_ID* -- HGNC gene ids to update outdated gene symbols in the HPO resources (can be annotated with VEP)

<br>

+ *CDS_position* -- Position of the variant on the transcript (can be annotated with VEP)
+ *Feature* -- Ensembl transcript ID used to match with the IDs in the canonical-transcript-file specified in the configuration file
