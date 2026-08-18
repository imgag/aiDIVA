# Running aiDIVA

aiDIVA expects a TAB separted table as input (see [here](https://github.com/imgag/aiDIVA/blob/master/doc/aidiva/run_aidiva.md#overview-of-necessary-and-optional-columns-in-the-input-table) for an overview of the necessary and optional columns) acompanied by a configuration file in the YAML format. Have a look in the `data` folder to get an example configuration file with placeholders. 

It is possible to run only parts of our software. In total we have three different modi in which our software can be run: *aiDIVA-RF*, *aiDIVA-meta*, *aiDIVA-meta-RF*.

Please be adviced that for the *aiDIVA-meta* mode you need the two evidence-based files *eb_dom.GSvar* and *eb_rec.GSvar*. These files can be created using the [VariantRanking](https://github.com/imgag/ngs-bits/blob/master/doc/tools/VariantRanking/index.md) tool which is part of the [ngs-bits](https://github.com/imgag/ngs-bits) tool collection. The VariantRanking tool needs as input a GSvar file which can be created if you process your VCF file with the [megSAP](https://github.com/imgag/megSAP) pipeline. Please head over to these repositories to have detailled instructions on how to use the tools.

To test aiDIVA-meta without the need to install the whole [megSAP](https://github.com/imgag/megSAP) pipeline we included a conversion script to convert the annotated table into a basic GSvar file that includes all necessary information needed to run the [VariantRanking](https://github.com/imgag/ngs-bits/blob/master/doc/tools/VariantRanking/index.md) tool (please be adviced that not all optional information used by [VariantRanking](https://github.com/imgag/ngs-bits/blob/master/doc/tools/VariantRanking/index.md) are included in htis converted file many of these informations are missing). This converison script is meant for initial tests, but is highly recommended to use the whole [megSAP](https://github.com/imgag/megSAP) pipeline to generate the GSvar files to use the full potential of the [VariantRanking](https://github.com/imgag/ngs-bits/blob/master/doc/tools/VariantRanking/index.md) tool. The [VariantRanking](https://github.com/imgag/ngs-bits/blob/master/doc/tools/VariantRanking/index.md) tool was written specifically with the [megSAP](https://github.com/imgag/megSAP) annotation in mind.

Head over to the [documentation](https://github.com/imgag/aiDIVA/blob/master/doc/annotation/run_annotation.md#convert-annotated-table-to-gsvar) of the annotion part to see how to use the conversion script.

Please be aware that the annotation from the [megSAP](https://github.com/imgag/megSAP) pipeline and the annotation script shipped in this repository are not identical. 


## Modes for Running aiDIVA

### aiDIVA-Score:

```
python3 run_aidiva-score.py --config configuration_aiDIVA.yaml --in_data input.tsv --out_prefix output_folder/aidiva_result [--workdir aidiva_workdir/] [--threads 1] [--log_file output_path/logs/aidiva_log.txt] [--log_level INFO]
```
mandatory parameters:

+ *config* -- YAML configuration file (in the `data` folder is an example configuration file)
+ *in\_data* -- TAB separated input table with the annotated variants
+ *out\_prefix* -- A prefix for the resulting output files the output folder can also be specified with that parameter

optional parameters:

+ *workdir* -- Working directory, where all temporary files are created and saved \[optional\]
+ *threads* -- Number of threads that should be used (default: 1) \[optional\]
+ *log\_file* -- Specify a custom log file to store the log messages from the tool \[optional\]
+ *log\_level* -- Define logging level \[DEBUG, INFO\] (default: INFO) \[optional\]


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

+ *\<your-result-prefix\>\_result\_aidiva-score.tsv* -- The result table with the predicted pathogenicity score (aiDIVA-Score) this table can be used as input for the other modes to skip the prediction part.
+ *\<your-result-prefix\>\_result\_aidiva-rf.tsv* -- The unfiltered result table (aiDIVA-RF).
+ *\<your-result-prefix\>\_result\_filtered\_aidiva-rf.tsv* -- The filtered result table (aiDIVA-RF) this table is also used for the subsequent analysis steps.
+ *\<your-result-prefix\>\_aidiva-rf\_based\_llm\_results.tsv* -- The LLM results based on the random forest-based ranking.
+ *\<your-result-prefix\>\_aidiva-eb-dom\_based\_llm\_result.tsv* -- The LLM results based on the evidence-based ranking (dominant model).
+ *\<your-result-prefix\>\_aidiva-eb-rec\_based\_llm\_result.tsv* -- The LLM results based on the evidence-based ranking (recessive model).
+ *\<your-result-prefix\>\_metascore\_results\_aidiva-meta-rf.tsv* -- The final result table with the variant ranking based on the meta-rf model (aiDIVA-meta-RF).
+ *\<your-result-prefix\>\_metascore\_results\_aidiva-meta.tsv* -- The final result table with the variant ranking based on the meta model (aiDIVA-meta).


## Overview of Necessary and Optional Columns in the Input Table

The following shows all columns that need to be present in the annotated input table for the software to work.

The column names should exactly match the column name specified in the table if not otherwise specified in the description.

### Strictly Necessary Columns

These columns give the basic information for each variant in the table.

+ *#CHROM* -- Chromosome identifier
+ *POS* -- Variant position 
+ *REF* -- Reference allele
+ *ALT* -- Observed alternative allele
+ *FILTER* -- Matches the FILTER column of a VCF file (checked terms: off-target, low_conf_region)
+ *SYMBOL* -- Approved gene symbol
+ *IMPACT* -- Impact of the variant according to ensembl (HIGH, MODERATE, LOW, MODIFIER)


### Necessary Feature Columns

These column names must match the feature-list specified in the configuration file. The following shows the column names for the feature-list specified in the example configuration given in the `data` folder.

+ *SIFT* -- SIFT score (float, minimum: 0.0, maximum 1.0)
+ *PolyPhen* -- PolyPhen2 score (float, minimum: 0.0, maximum 1.0)
+ *CADD_PHRED*  -- phred scaled CADD score (float, minimum: 0.0, maximum: ~40, no fixed upper limit)
+ *REVEL* -- REVEL score (float, minimum: 0.0, maximum: 1.0)
+ *MAX_AF* -- maximum allele frequency (float, minimum: 0.0, maximum: 1.0)
+ *EIGEN_PHRED* -- phred scaled Eigen score (float, minimum: 0.0, maximum: ~30, no fixed upper limit)
+ *CONDEL* -- CONDEL score (float, minimum: 0.0, maximum: 1.0)
+ *FATHMM_XF* -- FATHMM_XF score (float, minimum: 0.0, maximum: 1.0)
+ *MutationAssessor* -- MutationAssessor score (float, minimum: ~-5.5, maximum: ~6.0, no fixed range)
+ *phastCons_mammal* -- phastCons conservation score (float, minimum: 0.0, maximum: 1.0)
+ *phastCons_primate* -- phastCons conservation score (float, minimum: 0.0, maximum: 1.0)
+ *phastCons_vertebrate* -- phastCons conservation score (float, minimum: 0.0, maximum: 1.0)
+ *phyloP_mammal* -- phyloP conservation score (float, minimum: ~-10, maximum: ~10, no fixed range)
+ *phyloP_primate* -- phyloP conservation score (float, minimum: ~-10, maximum: ~10, no fixed range)
+ *phyloP_vertebrate* -- phyloP conservation score (float, minimum: ~-10, maximum: ~10, no fixed range)
+ *oe_lof* -- observed/expected loss-of-function ratio (float, minimum: 0.0, maximum: ~1.0, no fixed upper limit)
+ *homAF* -- homozygous allele frequency (float, minimum: 0.0, maximum: 1.0)
+ *CAPICE* -- CAPICE score (float, minimum: 0.0, maximum: 1.0)
+ *ALPHA_MISSENSE_SCORE* -- AlphaMissense score (float, minimum: 0.0, maximum: 1.0)
+ *HIGH_IMPACT* -- binary score indicating if a variants impact is HIGH or not (0=False, 1=True)
+ *IS_INDEL* --  binary score indicating if a variant is an indel (0=False, 1=True)


### Necessary Allele Frequency Columns

These columns are necessary if the MAX_AF column from the feature list is not present in the table. Not restricted to the columns shown below. column names must match the entries in the allele-frequency-list specified in the configuration file.

+ *gnomAD_AFR_AF* -- allele frequency in the african subpopulation (float, minimum: 0.0, maximum: 1.0)
+ *gnomAD_AMR_AF* -- allele frequency in the american subpopulation (float, minimum: 0.0, maximum: 1.0)
+ *gnomAD_EAS_AF* -- allele frequency in the east asian subpopulation (float, minimum: 0.0, maximum: 1.0)
+ *gnomAD_NFE_AF* -- allele frequency in the non-finnish european subpopulation (float, minimum: 0.0, maximum: 1.0)
+ *gnomAD_SAS_AF* -- allele frequency in the south asian subpopulation (float, minimum: 0.0, maximum: 1.0)


### Necessary Splicing Features

Scores used to handle splicing variants.

+ *SpliceAI* -- Single SpliceAI score per variant (we use the maximum of the four scores)


### Necessary Variant Consequence Information

Variant consequence terms annotated with ensembl VEP. The terms are shown in the overview on the ensembl [website](https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html). This site was also used as reference to determine the severity of the different terms.

+ *Consequence* -- Consequence terms annotated by VEP
+ *MOST_SEVERE_CONSEQUENCE* -- Single consequence term per variant (if overlapping terms were present we chose the most severe one)


### Necessary Sample Information
Information on the genotype of the variant. \<sample-id\> must match the one you pass with the *--sample\_id* parameter.

+ *GT_\<sample-id\>* -- Genotype information of the sample (\<sample-id\> must match the id of the sample you are currently analyzing)


### Optional Columns

The following columns are not necessarily needed to run aiDIVA.

+ *CLINVAR_DETAILS* -- ClinVar classification of the variant (needed if the parameter --skip_db_check is not given)
+ *HGMD_CLASS* -- HGMD classification of the variant (needed if the parameter --skip_db_check is not given)

<br>

+ *segmentDuplication* -- Segment duplciation annotation if the variant lies in a segment duplicaiton region.
+ *simpleRepeat* -- Simple repeat annotation if the variant lies inside a tandem repeat region.

<br>

+ *REPEATMASKER* -- RepeatMasker annotation (is used as additional filter in the prioritization step of aiDIVA-RF)

<br>

+ *low_conf_region* -- low confidence region annotation (is used as additional filter in the prioritization step of aiDIVA-RF)

<br>

+ *HGNC_ID* -- HGNC gene ids to update outdated gene symbols in the HPO resources (can be annotated with VEP)

<br>

+ *CDS_position* -- Position of the variant on the transcript (can be annotated with VEP)
+ *Feature* -- Ensembl transcript ID used to match with the IDs in the canonical-transcript-file specified in the configuration file
