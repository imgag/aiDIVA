# aiDIVA - augmented intelligence-based DIsease Variant Analysis

aiDIVA is an analysis pipeline that combines a pathogenicity-based approach and an optional evidence-based approach with state-of-the-art large language models (LLMs) to identify potential disease causing variants in a given rare disease sample.

aiDIVA comprises the following steps:

1. A pathogenicity-based approach that utilizes the predictions of a random forest model trained on ClinVar data, supplemented with phenotype information given as HPO terms.

2. An optional evidence-based approach includes the ranks and scores that can be obtained from the [VariantRanking](https://github.com/imgag/ngs-bits/blob/master/doc/tools/VariantRanking/index.md) tool included in [ngs-bits](https://github.com/imgag/ngs-bits/).

3. State-of-the-art LLMs are utilized to refine the ranking results from the previous two approaches.

4. A final meta model is used to combine all preliminary results and create a final ranking of the variants. For this model we also use a random forest.


## System Requirements
The program is written in Python 3

Tested version: 3.12.3

The following additional libraries need to be installed in order to use the program:

+ networkx (v3.4.2)
+ numpy (v1.26.4)
+ openai (v1.60.2)
+ pandas (v2.2.3)
+ pysam (v0.22.1)
+ pyyaml (v6.0.2)
+ scipy (v1.15.1)
+ scikit-learn (v1.3.2)

For easy package installation in your Python virtual environment we included a requirements.txt just run `pip install -r requiremnts.txt` to install all necessary packages.

If a newer scikit-learn version is used it is advised to create a new model with the newer scikit version.


## Annotation Resources and Tools
If the annotation is not made beforehand make sure that the necessary database resources and tools are present on your system and the paths in the configuration file are set correctly.

For detailed instructions on how to download and prepare the annotation resources please head over to the respective [documentation](https://github.com/imgag/aiDIVA/blob/master/doc/prepare_annotation_resources.md) (prepare\_annotation\_resources.md) found in the _doc_ folder.

For detailed instructions on how to download and install the annotation tools please head over to the respective [documentation](https://github.com/imgag/aiDIVA/blob/master/doc/install_additional_tools.md) (install\_additional\_tools.md) found in the _doc_ folder.


## HPO Resources
The HPO resources needed for the prioritization step can be found in the `data` folder. The path to the files is specified in the configuration file make sure that it leads to the correct location. Although we used networkx v3.4 to create the HPO graph it should also be possible to use networkx in version 1 or 2. We included some workarounds to still support the older versions.

To recreate the HPO resources please head over to the detailed [instructions](https://github.com/imgag/aiDIVA/blob/master/doc/recreate_hpo_resources.md) found in the _doc_ folder.

Protein interactions based on String-DB v11.0b

Last update of HPO resources: 22nd May, 2023


## Pathogenicity Prediction
There is one random forest model that is used in aiDIVA to predict the pathogenicity of a given variant. It is a combined model for SNV and inframe InDel variants. The training data of the model consists of variants from Clinvar.

The scripts used to train the model can be found in the following GitHub repository: [aiDIVA-Training](https://github.com/imgag/aiDIVA-Training)

_Frameshift_ variants will get a default score of 0.9, whereas _synonymous_ variants always get the lowest score 0.0

A pretrained random forest model (_aidiva-rf_) using our current feature set can be found [here](https://download.imgag.de/aidiva/aidiva_pretrained_models/). The latest model was trained using scikit-learn v1.3.2. The trained models of scikit-learn are version dependent.


## LLM Usage
aiDIVA supports the use of the official [OpenAI API](https://platform.openai.com/docs/api-reference/introduction) to send the requests to GPT-4o or GPT-4.1 for example. To use the OpenAI API you need an account and an API-Key that needs to be specified in the configuration file.
Alternatively it is possible to set up your own local LLM (eg., LLama-8b, Mistral-12b, ...) and provide it locally as a Webservice. For an easy deployment you could use the NVIDIA NIM Containers see [here](https://build.nvidia.com/meta/llama-3_1-8b-instruct/deploy) for more details on how to do that. These local LLMs use the same python package for inference you just have to specify the port and URL where to find the local model in the configuration file.


## Meta Model
You can download the pretrained meta models (_aidiva-meta & aidiva-meta-rf_) [here](https://download.imgag.de/aidiva/aidiva_pretrained_models/). For these two models we used a random forest model that takes as features the ranking position and scores from the initial rankings (pathogenicity-based and evidence-based) plus the ranking result from the LLMs and the inheritance mode used in the evidence-based model.


## Running aiDIVA
aiDIVA can be run either on already annotated VCF files or unannotated VCF files. In both cases a configuration file in the YAML format is required. Furthermore it is possible to run only segments of the whole pipeline. In total we have six different modi in which our pipeline can be run.

Make sure to specify the trained model in the configuration file `aiDIVA_configuration_annotated.yaml` or `aiDIVA_configuration_with_annotation.yaml` respectively.

In the *data* folder you can find example configuration files with placeholders.


### Running Only Annotation:

```
python3 run_annotation.py --config aiDIVA_configuration_with_annotation.yaml --vcf input.vcf --out_folder output_path/aidiva_annotation [--filtered] [--filtered_folder output_path/aidiva_filtered/] [--inhouse_sample] [--output_table] [--threads 1] [--log_file output_path/logs/aidiva_log.txt] [--log_level INFO]
```

+ *config* -- YAML configuration file (in the `data` folder there are example configuration files)
+ *vcf* -- Input VCF file containing all variants of the given patient
+ *out_folder* -- Folder where the resulting output files should be saved
+ *filtered* -- Skip filtering step if it was already done [optional]
+ *filtered_folder* -- Folder where to find the already filtered VCF files [optional]
+ *inhouse_sample* -- The input VCF was prepared inhouse (skips some additional preparation steps to prevent possible problems) [optional]
+ *output_table* -- Save the annotations additionally as a TAB separated file [optional]
+ *threads* -- Number of threads that should be used (default: 1) [optional]
+ *log_file* -- Specify a custom log file to store the log messages from the tool [optional]
+ *log_level* -- Define logging level \[DEBUG, INFO\] (default: INFO) [optional]



### Running aiDIVA Without Meta Ranking on Already Annotated Data:

```
python3 run_aiDIVA.py --config aiDIVA_configuration_annotated.yaml --snp_vcf annotated_snp.vcf --indel_vcf annotated_indel.vcf --expanded_indel_vcf annotated_expanded_indel.vcf --out_prefix aidiva_result --workdir aidiva_workdir/ [--hpo_list HP:xxxx,HP:xxxx] [--gene_exclusion gene_exclusion.txt] [--family_file family.txt] [--family_type SINGLE] [--skip_db_check] [--only_top_results] [--top_rank 25] [--threads 1] [--log_file output_path/logs/aidiva_log.txt] [--log_level INFO] [--save_as_vcf]
```

+ *config* -- YAML configuration file (in the `data` folder there are example configuration files)
+ *snp_vcf* -- Input VCF file containing all annotated SNP variants of the given patient
+ *indel_vcf* -- Input VCF file containing all InDel variants with basic annotation of the given patient
+ *expanded\_indel\_vcf* -- Input VCF file containing all annotated expanded InDel variants of the given patient
+ *out_prefix* -- A prefix for the resulting output files the output folder can also be specified with that parameter
+ *workdir* -- Working directory, where all temporary files are created and saved [optional]
+ *hpo_list* -- Comma-separated list of HPO terms observed with the patient [optional]
+ *gene_exclusion* -- TXT file containing genes that should be excluded during the analysis of the HPO relatedness [optional]
+ *family_file* -- TXT file containing the sample information if run on multisample VCF files [optional]
+ *family_type* -- Type of the family relation \[SINGLE, TRIO\] (default: SINGLE) [optional]
+ *skip\_db\_check* -- Skip the database checkup for existing entries in ClinVar (and HGMD) [optional]
+ *only\_top\_results* -- Restrict the results to only report the top X variants (default: 25) [optional]
+ *top_rank* -- Rank that should be used as maximum for the top ranking results to report (default: 25) [optional]
+ *threads* -- Number of threads that should be used (default: 1) [optional]
+ *log_file* -- Specify a custom log file to store the log messages from the tool [optional]
+ *log_level* -- Define logging level \[DEBUG, INFO\] (default: INFO) [optional]


### Running aiDIVA  Without Meta Ranking and Perform the Annotation:

```
python3 run_annotation_and_aiDIVA.py --config aiDIVA_configuration_with_annotation.yaml --vcf input.vcf --out_prefix output_path/aidiva_result [--workdir aidiva_workdir/] [--hpo_list HP:xxxx,HP:xxxx] [--gene_exclusion gene_exclusion.txt] [--family_file family.txt] [--family_type SINGLE] [--skip_db_check] [--only_top_results] [--top_rank 25] [--inhouse_sample] [--threads 1] [--log_file output_path/logs/aidiva_log.txt] [--log_level INFO] [--save_as_vcf]
```

+ *config* -- YAML configuration file (in the `data` folder there are example configuration files)
+ *vcf* -- Input VCF file containing all variants of the given patient
+ *out_prefix* -- A prefix for the resulting output files the output folder can also be specified with that parameter
+ *workdir* -- Working directory, where all temporary files are created and saved [optional]
+ *hpo_list* -- Comma-separated list of HPO terms observed with the patient [optional]
+ *gene_exclusion* -- TXT file containing genes that should be excluded during the analysis of the HPO relatedness [optional]
+ *family_file* -- TXT file containing the sample information if run on multisample VCF files [optional]
+ *family_type* -- Type of the family relation \[SINGLE, TRIO\] (default: SINGLE) [optional]
+ *skip\_db\_check* -- Skip the database checkup for existing entries in ClinVar (and HGMD) [optional]
+ *only\_top\_results* -- Restrict the results to only report the top X variants (default: 25) [optional]
+ *top_rank* -- Rank that should be used as maximum for the top ranking results to report (default: 25) [optional]
+ *inhouse_sample* -- The input VCF was prepared inhouse (skips some additional preparation steps to prevent possible problems) [optional]
+ *threads* -- Number of threads that should be used (default: 1) [optional]
+ *log_file* -- Specify a custom log file to store the log messages from the tool [optional]
+ *log_level* -- Define logging level \[DEBUG, INFO\] (default: INFO) [optional]
+ *save\_as\_vcf* -- Save results additionally in VCF format [optional]


### Running Only the Metascore Ranking:

```
python3 run_metamodel.py --config aiDIVA_configuration_with_annotation.yaml --in_rf in_rf.tsv --out_prefix output_path/aidiva_result --sample_id NA12878 [--in_eb_dom in_eb_dom.GSvar] [--in_eb_rec in_eb_rec.GSvar] [--workdir aidiva_workdir/] --hpo_list HP:xxxx,HP:xxxx [--gender male/female] [--age 42] [--evidence_based] [--threads 1] [--log_file output_path/logs/aidiva_log.txt] [--log_level INFO]
```

+ *config* -- YAML configuration file (in the `data` folder there are example configuration files)
+ *in_rf* -- TAB separated file with the random forest-based ranking results from aiDIVA without using the metamodel
+ *out_prefix* -- A prefix for the resulting output files the output folder can also be specified with that parameter
+ *sample_id* -- Sample ID this is used to extract the genotype in the VCF file
+ *in\_eb\_dom* -- GSvar file containing the evidence-based ranking results of the patients variants using the dominant mode of the algorithm [optional]
+ *in\_eb\_rec* -- GSvar file containing the evidence-based ranking results of the patients variants using the recessive mode of the algorithm [optional]
+ *workdir* -- Working directory, where all temporary files are created and saved [optional]
+ *hpo_list* -- Comma-separated list of HPO terms observed with the patient [optional]
+ *gender* -- Gender of the patient if known [optional]
+ *age* -- Age of the patient if known [optional]
+ *evidence_based* -- Include evidence-based ranking results [optional]
+ *threads* -- Number of threads that should be used (default: 1) [optional]
+ *log_file* -- Specify a custom log file to store the log messages from the tool [optional]
+ *log_level* -- Define logging level \[DEBUG, INFO\] (default: INFO) [optional]


### Running aiDIVA on Already Annotated Data With Metascore Ranking:

```
python3 run_aiDIVA_with_metamodel.py --config aiDIVA_configuration_with_annotation.yaml --snp_vcf annotated_snp.vcf --indel_vcf annotated_indel.vcf --expanded_indel_vcf annotated_expanded_indel.vcf [--in_eb_dom in_eb_dom.GSvar] [--in_eb_rec in_eb_rec.GSvar] --out_prefix output_path/aidiva_result --sample_id NA12878 [--workdir aidiva_workdir/] --hpo_list HP:xxxx,HP:xxxx [--gene_exclusion gene_exclusion.txt] [--family_file family.txt] [--family_type SINGLE] [--skip_db_check] [--only_top_results] [--top_rank 25] [--gender male/female] [--age 42] [--evidence_based] [--threads 1] [--log_file output_path/logs/aidiva_log.txt] [--log_level INFO] [--save_as_vcf]
```

+ *config* -- YAML configuration file (in the `data` folder there are example configuration files)
+ *snp_vcf* -- Input VCF file containing all annotated SNP variants of the given patient
+ *indel_vcf* -- Input VCF file containing all InDel variants with basic annotation of the given patient
+ *expanded\_indel\_vcf* -- Input VCF file containing all annotated expanded InDel variants of the given patient
+ *in\_eb\_dom* -- GSvar file containing the evidence-based ranking results of the patients variants using the dominant mode of the algorithm [optional]
+ *in\_eb\_rec* -- GSvar file containing the evidence-based ranking results of the patients variants using the recessive mode of the algorithm [optional]
+ *out_prefix* -- A prefix for the resulting output files the output folder can also be specified with that parameter
+ *sample_id* -- Sample ID this is used to extract the genotype in the VCF file
+ *workdir* -- Working directory, where all temporary files are created and saved [optional]
+ *hpo_list* -- Comma-separated list of HPO terms observed with the patient
+ *gene_exclusion* -- TXT file containing genes that should be excluded during the analysis of the HPO relatedness [optional]
+ *family_file* -- TXT file containing the sample information if run on multisample VCF files [optional]
+ *family_type* -- Type of the family relation \[SINGLE, TRIO\] (default: SINGLE) [optional]
+ *skip\_db\_check* -- Skip the database checkup for existing entries in ClinVar (and HGMD) [optional]
+ *only\_top\_results* -- Restrict the results to only report the top X variants (default: 25) [optional]
+ *top_rank* -- Rank that should be used as maximum for the top ranking results to report (default: 25) [optional]
+ *gender* -- Gender of the patient if known [optional]
+ *age* -- Age of the patient if known [optional]
+ *evidence_based* -- Include evidence-based ranking results [optional]
+ *threads* -- Number of threads that should be used (default: 1) [optional]
+ *log_file* -- Specify a custom log file to store the log messages from the tool [optional]
+ *log_level* -- Define logging level \[DEBUG, INFO\] (default: INFO) [optional]
+ *save\_as\_vcf* -- Save results additionally in VCF format [optional]


### Running aiDIVA With Metascore Ranking and Perform the Annotation:

```
python3 run_annotation_and_aiDIVA_with_metamodel.py --config aiDIVA_configuration_with_annotation.yaml --vcf input.vcf [--in_eb_dom in_eb_dom.GSvar] [--in_eb_rec in_eb_rec.GSvar] --out_prefix output_path/aidiva_result --sample_id NA12878 [--workdir aidiva_workdir/] --hpo_list HP:xxxx,HP:xxxx [--gene_exclusion gene_exclusion.txt] [--family_file family.txt] [--family_type SINGLE] [--skip_db_check] [--only_top_results] [--top_rank 25] [--inhouse_sample] [--gender male/female] [--age 42] [--evidence_based] [--threads 1] [--log_file output_path/logs/aidiva_log.txt] [--log_level INFO] [--save_as_vcf]
```

+ *config* -- YAML configuration file (in the `data` folder there are example configuration files)
+ *vcf* -- Input VCF file containing all variants of the given patient
+ *in\_eb\_dom* -- GSvar file containing the evidence-based ranking results of the patients variants using the dominant mode of the algorithm [optional]
+ *in\_eb\_rec* -- GSvar file containing the evidence-based ranking results of the patients variants using the recessive mode of the algorithm [optional]
+ *out_prefix* -- A prefix for the resulting output files the output folder can also be specified with that parameter
+ *sample_id* -- Sample ID this is used to extract the genotype in the VCF file
+ *workdir* -- Working directory, where all temporary files are created and saved [optional]
+ *hpo_list* -- Comma-separated list of HPO terms observed with the patient [optional]
+ *gene_exclusion* -- TXT file containing genes that should be excluded during the analysis of the HPO relatedness [optional]
+ *family_file* -- TXT file containing the sample information if run on multisample VCF files [optional]
+ *family_type* -- Type of the family relation \[SINGLE, TRIO\] (default: SINGLE) [optional]
+ *skip\_db\_check* -- Skip the database checkup for existing entries in ClinVar (and HGMD) [optional]
+ *only\_top\_results* -- Restrict the results to only report the top X variants (default: 25) [optional]
+ *top_rank* -- Rank that should be used as maximum for the top ranking results to report (default: 25) [optional]
+ *inhouse_sample* -- The input VCF was prepared inhouse (skips some additional preparation steps to prevent possible problems) [optional]
+ *gender* -- Gender of the patient if known [optional]
+ *age* -- Age of the patient if known [optional]
+ *evidence_based* -- Include evidence-based ranking results [optional]
+ *threads* -- Number of threads that should be used (default: 1) [optional]
+ *log_file* -- Specify a custom log file to store the log messages from the tool [optional]
+ *log_level* -- Define logging level \[DEBUG, INFO\] (default: INFO) [optional]
+ *save\_as\_vcf* -- Save results additionally in VCF format [optional]


## aiDIVA Results
aiDIVA will produce multiple different output files. The following lists all possible result files. Depending on your chosen mode to run aiDIVA you will only get a subset of these result files.

+ \<your-result-prefix\>\_aidiva\_result.tsv -- The unfiltered result table for the random forest-based ranking (aiDIVA-RF).
+ \<your-result-prefix\>\_aidiva\_result\_filtered.tsv -- The filtered result table for the random forest-based ranking (aiDIVA-RF) this table is internally used for the subsequent analysis steps (if they are performed).
+ \<your-result-prefix\>\_aidiva\_random\_forest\_based\_llm\_result.tsv -- The LLM results based on the random forest-based top ranking variants.
+ \<your-result-prefix\>\_aidiva\_evidence\_dominant\_based\_llm\_result.tsv -- The LLM results based on the evidence-based top ranking variants (dominant model).
+ \<your-result-prefix\>\_aidiva\_evidence\_recessive\_based\_llm\_result.tsv -- The LLM results based on the evidence-based top ranking variants (recessive model).
+ \<your-result-prefix\>\_aidiva\_metascore\_results.tsv -- The final result table with the variant ranking based on the meta model (aiDIVA-meta).

If aiDIVA is run in annotation mode you will only get the following annotated files:

+ \<your-result-prefix\>\_snp\_annotated.vcf
+ \<your-result-prefix\>\_indel\_annotated.vcf
+ \<your-result-prefix\>\_indelExpanded\_annotated.vcf
