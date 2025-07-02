# aiDIVA - augmented intelligence-based Disease Variant Analysis

aiDIVA is an analysis pipeline that predicts the pathogenicity of given variants and prioritizes them according to the predicted pathogenicity to identify potential disease causing variants of a given sample.

The pathogenicity prediction is based on a random forest (RF) model. That covers both SNV and inframe InDel variants (frameshift variants are not covered, they are handled separately).

Later we combine this prediction-based ranking with an optional evidence-based approach and refine it with an LLM.

In the end we create a metascore-based ranking, where we use a random forest with the ranking results as features.


## System requirements
The program is written in Python 3 (v3.12.3)

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


## Annotation resources and tools
If the annotation is not made beforehand make sure that the necessary database resources and tools are present on your system and the paths in the configuration file are set correctly.

For detailed instructions on how to download and prepare the annotation resources please head over to the respective [documentation](https://github.com/imgag/aiDIVA/blob/master/doc/prepare_annotation_resources.md) (prepare_annotation_resources.md) found in the _doc_ folder.

For detailed instructions on how to download and install the annotation tools please head over to the respective [documentation](https://github.com/imgag/aiDIVA/blob/master/doc/install_additional_tools.md) (isntall_additional_tools.md) found in the _doc_ folder.


## HPO resources
The HPO resources needed for the prioritization step can be found in the `data` folder. The path to the files is specified in the configuration file make sure that it leads to the correct location. Although we used networkx v3.4 to create the HPO graph it should also be possible to use networkx in version 1 or 2. We included some workarounds to still support the older versions.

To recreate the HPO resources please head over to the detailed [instructions](https://github.com/imgag/aiDIVA/blob/master/doc/recreate_hpo_resources.md) found in the _doc_ folder.
<br>
<br>

Protein interactions based on String-DB v11.0b

Last update of HPO resources: 22nd May, 2023


## LLM setup
aiDIVA supports the use of the official [OpenAI API](https://platform.openai.com/docs/api-reference/introduction) to send the requests to GPT-4o or GPT-4.1 for example. For that you need an account and an API-Key that needs to be specified in the configuration file.
Alternatively it is possible to set up your own local LLM (eg., LLama-8b, Mistral-12b, ...) and serve it as a local Webservice. For an easy deployment you could use the NVIDIA NIM Containers see [here](https://build.nvidia.com/meta/llama-3_1-8b-instruct/deploy) for more details on how to do that. These local LLMs use the same python package for inference you just have to specify the port and URL where to find the local model in the configuration file.


## Pathogenicity prediction
There is one random forest model that is used in aiDIVA to predict the pathogenicity of a given variant. It is a combined model for SNV and inframe InDel variants. The training data of the model consists of variants from Clinvar.

The scripts used to train the model can be found in the following GitHub repository: [aiDIVA-Training](https://github.com/imgag/aiDIVA-Training)

_Frameshift_ variants will get a default score of 0.9, whereas _synonymous_ variants always get the lowest score 0.0

A pretrained random forest model (aiDIVA-RF) using our current feature set can be found [here](https://download.imgag.de/aidiva/aidiva_pretrained_models/). The latest model was trained using scikit-learn v1.3.2. The trained models of scikit-learn are version dependent.


## Meta model
You can download the pretrained meta models (aiDIVA-meta & aiDIVA-meta-RF) [here](https://download.imgag.de/aidiva/aidiva_pretrained_models/). For these two models we used a random forest model that takes as features the ranking position and scores from the initial rankings (random forest-based and evidence-based) plus the ranking result from the LLMs and the inheritance mode used in the evidence-based model.


## Running aiDIVA
aiDIVA can be run either on already annotated VCF files or unannotated VCF files. In both cases a configuration file in the YAML format is required. Furthermore it is possible to run only segments of the whole pipeline. In total we have six different modi in which our pipeline can be run.

Make sure to specify the trained model in the config file `aiDIVA_configuration_annotated.yaml`/`aiDIVA_configuration_with_annotation.yaml`.


### Running only annotation:

```
python3 run_annotation.py --config aiDIVA_configuration_with_annotation.yaml --vcf input.vcf --out_folder output_path/aidiva_annotation [--filtered] [--filtered_folder output_path/aidiva_filtered/] [--inhouse_sample] [--output_table] [--threads 1] [--log_file output_path/logs/aidiva_log.txt] [--log_level INFO]
```

+ _config_ -- YAML configuration file (in the `data` folder there are example configuration files)
+ _vcf_ -- Input VCF file containing all variants of the given patient
+ _out_folder_ -- Folder where the resulting output files should be saved
+ _filtered_ -- Skip filtering step if it was already done [optional]
+ _filtered_folder_ -- Folder where to find the already filtered VCF files [optional]
+ _inhouse_sample_ -- The input VCF was prepared inhouse (skips some additional preparation steps to prevent possible problems) [optional]
+ _output_table_ -- Save the annotations additionally as a TAB separated file [optional]
+ _threads_ -- Number of threads that should be used (default: 1) [optional]
+ _log_file_ -- Specify a custom log file to store the log messages from the tool [optional]
+ _log_level_ -- Define logging level [DEBUG, INFO] (default: INFO) [optional]



### Running aiDIVA without meta ranking on already annotated data:

```
python3 run_aiDIVA.py --config aiDIVA_configuration_annotated.yaml --snp_vcf annotated_snp.vcf --indel_vcf annotated_indel.vcf --expanded_indel_vcf annotated_expanded_indel.vcf --out_prefix aidiva_result --workdir aidiva_workdir/ [--hpo_list HP:xxxx,HP:xxxx] [--gene_exclusion gene_exclusion.txt] [--family_file family.txt] [--family_type SINGLE] [--skip_db_check] [--only_top_results] [--top_rank 25] [--threads 1] [--log_file output_path/logs/aidiva_log.txt] [--log_level INFO] [--save_as_vcf]
```

+ _config_ -- YAML configuration file (in the `data` folder there are example configuration files)
+ _snp_vcf_ -- Input VCF file containing all annotated SNP variants of the given patient
+ _indel_vcf_ -- Input VCF file containing all InDel variants with basic annotation of the given patient
+ _expanded_indel_vcf_ -- Input VCF file containing all annotated expanded InDel variants of the given patient
+ _out_prefix_ -- A prefix for the resulting output files the output folder can also be specified with that parameter
+ _workdir_ -- Working directory, where all temporary files are created and saved [optional]
+ _hpo_list_ -- Comma-separated list of HPO terms observed with the patient [optional]
+ _gene_exclusion_ -- TXT file containing genes that should be excluded during the analysis of the HPO relatedness [optional]
+ _family_file_ -- TXT file containing the sample information if run on multisample VCF files [optional]
+ _family_type_ -- Type of the family relation [SINGLE, TRIO] (default: SINGLE) [optional]
+ _skip_db_check_ -- Skip the database checkup for existing entries in ClinVar (and HGMD) [optional]
+ _only_top_results_ -- Restrict the results to only report the top X variants (default: 25) [optional]
+ _top_rank_ -- Rank that should be used as maximum for the top ranking results to report (default: 25) [optional]
+ _threads_ -- Number of threads that should be used (default: 1) [optional]
+ _log_file_ -- Specify a custom log file to store the log messages from the tool [optional]
+ _log_level_ -- Define logging level [DEBUG, INFO] (default: INFO) [optional]


### Running aiDIVA  without meta ranking and perform the annotation:

```
python3 run_annotation_and_aiDIVA.py --config aiDIVA_configuration_with_annotation.yaml --vcf input.vcf --out_prefix output_path/aidiva_result [--workdir aidiva_workdir/] [--hpo_list HP:xxxx,HP:xxxx] [--gene_exclusion gene_exclusion.txt] [--family_file family.txt] [--family_type SINGLE] [--skip_db_check] [--only_top_results] [--top_rank 25] [--inhouse_sample] [--threads 1] [--log_file output_path/logs/aidiva_log.txt] [--log_level INFO] [--save_as_vcf]
```

+ _config_ -- YAML configuration file (in the `data` folder there are example configuration files)
+ _vcf_ -- Input VCF file containing all variants of the given patient
+ _out_prefix_ -- A prefix for the resulting output files the output folder can also be specified with that parameter
+ _workdir_ -- Working directory, where all temporary files are created and saved [optional]
+ _hpo_list_ -- Comma-separated list of HPO terms observed with the patient [optional]
+ _gene_exclusion_ -- TXT file containing genes that should be excluded during the analysis of the HPO relatedness [optional]
+ _family_file_ -- TXT file containing the sample information if run on multisample VCF files [optional]
+ _family_type_ -- Type of the family relation [SINGLE, TRIO] (default: SINGLE) [optional]
+ _skip_db_check_ -- Skip the database checkup for existing entries in ClinVar (and HGMD) [optional]
+ _only_top_results_ -- Restrict the results to only report the top X variants (default: 25) [optional]
+ _top_rank_ -- Rank that should be used as maximum for the top ranking results to report (default: 25) [optional]
+ _inhouse_sample_ -- The input VCF was prepared inhouse (skips some additional preparation steps to prevent possible problems) [optional]
+ _threads_ -- Number of threads that should be used (default: 1) [optional]
+ _log_file_ -- Specify a custom log file to store the log messages from the tool [optional]
+ _log_level_ -- Define logging level [DEBUG, INFO] (default: INFO) [optional]
+ _save_as_vcf_ -- Save results additionally in VCF format [optional]


### Running only the metascore ranking:

```
python3 run_metamodel.py --config aiDIVA_configuration_with_annotation.yaml --in_rf in_rf.tsv --out_prefix output_path/aidiva_result --sample_id NA12878 [--in_eb_dom in_eb_dom.GSvar] [--in_eb_rec in_eb_rec.GSvar] [--workdir aidiva_workdir/] --hpo_list HP:xxxx,HP:xxxx [--gender male/female] [--age 42] [--evidence_based] [--threads 1] [--log_file output_path/logs/aidiva_log.txt] [--log_level INFO]
```

+ _config_ -- YAML configuration file (in the `data` folder there are example configuration files)
+ _in_rf_ -- TAB separated file with the random forest-based ranking results from aiDIVA without using the metamodel
+ _out_prefix_ -- A prefix for the resulting output files the output folder can also be specified with that parameter
+ _sample_id_ -- Sample ID this is used to extract the genotype in the VCF file
+ _in_eb_dom_ -- GSvar file containing the evidence-based ranking results of the patients variants using the dominant mode of the algorithm [optional]
+ _in_eb_rec_ -- GSvar file containing the evidence-based ranking results of the patients variants using the recessive mode of the algorithm [optional]
+ _workdir_ -- Working directory, where all temporary files are created and saved [optional]
+ _hpo_list_ -- Comma-separated list of HPO terms observed with the patient [optional]
+ _gender_ -- Gender of the patient if known [optional]
+ _age_ -- Age of the patient if known [optional]
+ _evidence_based_ -- Include evidence-based ranking results [optional]
+ _threads_ -- Number of threads that should be used (default: 1) [optional]
+ _log_file_ -- Specify a custom log file to store the log messages from the tool [optional]
+ _log_level_ -- Define logging level [DEBUG, INFO] (default: INFO) [optional]


### Running aiDIVA on already annotated data with metascore ranking:

```
python3 run_aiDIVA_with_metamodel.py --config aiDIVA_configuration_with_annotation.yaml --snp_vcf annotated_snp.vcf --indel_vcf annotated_indel.vcf --expanded_indel_vcf annotated_expanded_indel.vcf [--in_eb_dom in_eb_dom.GSvar] [--in_eb_rec in_eb_rec.GSvar] --out_prefix output_path/aidiva_result --sample_id NA12878 [--workdir aidiva_workdir/] --hpo_list HP:xxxx,HP:xxxx [--gene_exclusion gene_exclusion.txt] [--family_file family.txt] [--family_type SINGLE] [--skip_db_check] [--only_top_results] [--top_rank 25] [--gender male/female] [--age 42] [--evidence_based] [--threads 1] [--log_file output_path/logs/aidiva_log.txt] [--log_level INFO] [--save_as_vcf]
```

+ _config_ -- YAML configuration file (in the `data` folder there are example configuration files)
+ _snp_vcf_ -- Input VCF file containing all annotated SNP variants of the given patient
+ _indel_vcf_ -- Input VCF file containing all InDel variants with basic annotation of the given patient
+ _expanded_indel_vcf_ -- Input VCF file containing all annotated expanded InDel variants of the given patient
+ _in_eb_dom_ -- GSvar file containing the evidence-based ranking results of the patients variants using the dominant mode of the algorithm [optional]
+ _in_eb_rec_ -- GSvar file containing the evidence-based ranking results of the patients variants using the recessive mode of the algorithm [optional]
+ _out_prefix_ -- A prefix for the resulting output files the output folder can also be specified with that parameter
+ _sample_id_ -- Sample ID this is used to extract the genotype in the VCF file
+ _workdir_ -- Working directory, where all temporary files are created and saved [optional]
+ _hpo_list_ -- Comma-separated list of HPO terms observed with the patient
+ _gene_exclusion_ -- TXT file containing genes that should be excluded during the analysis of the HPO relatedness [optional]
+ _family_file_ -- TXT file containing the sample information if run on multisample VCF files [optional]
+ _family_type_ -- Type of the family relation [SINGLE, TRIO] (default: SINGLE) [optional]
+ _skip_db_check_ -- Skip the database checkup for existing entries in ClinVar (and HGMD) [optional]
+ _only_top_results_ -- Restrict the results to only report the top X variants (default: 25) [optional]
+ _top_rank_ -- Rank that should be used as maximum for the top ranking results to report (default: 25) [optional]
+ _gender_ -- Gender of the patient if known [optional]
+ _age_ -- Age of the patient if known [optional]
+ _evidence_based_ -- Include evidence-based ranking results [optional]
+ _threads_ -- Number of threads that should be used (default: 1) [optional]
+ _log_file_ -- Specify a custom log file to store the log messages from the tool [optional]
+ _log_level_ -- Define logging level [DEBUG, INFO] (default: INFO) [optional]
+ _save_as_vcf_ -- Save results additionally in VCF format [optional]


### Running aiDIVA with metascore ranking and perform the annotation:

```
python3 run_annotation_and_aiDIVA_with_metamodel.py --config aiDIVA_configuration_with_annotation.yaml --vcf input.vcf [--in_eb_dom in_eb_dom.GSvar] [--in_eb_rec in_eb_rec.GSvar] --out_prefix output_path/aidiva_result --sample_id NA12878 [--workdir aidiva_workdir/] --hpo_list HP:xxxx,HP:xxxx [--gene_exclusion gene_exclusion.txt] [--family_file family.txt] [--family_type SINGLE] [--skip_db_check] [--only_top_results] [--top_rank 25] [--inhouse_sample] [--gender male/female] [--age 42] [--evidence_based] [--threads 1] [--log_file output_path/logs/aidiva_log.txt] [--log_level INFO] [--save_as_vcf]
```

+ _config_ -- YAML configuration file (in the `data` folder there are example configuration files)
+ _vcf_ -- Input VCF file containing all variants of the given patient
+ _in_eb_dom_ -- GSvar file containing the evidence-based ranking results of the patients variants using the dominant mode of the algorithm [optional]
+ _in_eb_rec_ -- GSvar file containing the evidence-based ranking results of the patients variants using the recessive mode of the algorithm [optional]
+ _out_prefix_ -- A prefix for the resulting output files the output folder can also be specified with that parameter
+ _sample_id_ -- Sample ID this is used to extract the genotype in the VCF file
+ _workdir_ -- Working directory, where all temporary files are created and saved [optional]
+ _hpo_list_ -- Comma-separated list of HPO terms observed with the patient [optional]
+ _gene_exclusion_ -- TXT file containing genes that should be excluded during the analysis of the HPO relatedness [optional]
+ _family_file_ -- TXT file containing the sample information if run on multisample VCF files [optional]
+ _family_type_ -- Type of the family relation [SINGLE, TRIO] (default: SINGLE) [optional]
+ _skip_db_check_ -- Skip the database checkup for existing entries in ClinVar (and HGMD) [optional]
+ _only_top_results_ -- Restrict the results to only report the top X variants (default: 25) [optional]
+ _top_rank_ -- Rank that should be used as maximum for the top ranking results to report (default: 25) [optional]
+ _inhouse_sample_ -- The input VCF was prepared inhouse (skips some additional preparation steps to prevent possible problems) [optional]
+ _gender_ -- Gender of the patient if known [optional]
+ _age_ -- Age of the patient if known [optional]
+ _evidence_based_ -- Include evidence-based ranking results [optional]
+ _threads_ -- Number of threads that should be used (default: 1) [optional]
+ _log_file_ -- Specify a custom log file to store the log messages from the tool [optional]
+ _log_level_ -- Define logging level [DEBUG, INFO] (default: INFO) [optional]
+ _save_as_vcf_ -- Save results additionally in VCF format [optional]


## aiDIVA results
aiDIVA will produce multiple different output files. The following lists all possible result files. Depending on your chosen mode to run aiDIVA you will only get a subset of these result files.

+ <your-result-prefix>_aidiva_result.tsv -- The unfiltered result table for the random forest-based ranking (aiDIVA-RF).
+ <your-result-prefix>_aidiva_result_filtered.tsv -- The filtered result table for the random forest-based ranking (aiDIVA-RF) this table is internally used for the subsequent analysis steps (if they are performed).
+ <your-result-prefix>_aidiva_random_forest_based_llm_result.tsv -- The LLM results based on the random forest-based top ranking variants.
+ <your-result-prefix>_aidiva_evidence_dominant_based_llm_result.tsv -- The LLM results based on the evidence-based top ranking variants (dominant model).
+ <your-result-prefix>_aidiva_evidence_recessive_based_llm_result.tsv -- The LLM results based on the evidence-based top ranking variants (recessive model).
+ <your-result-prefix>_aidiva_metascore_results.tsv -- The final result table with the variant ranking based on the meta model (aiDIVA-meta).

If aiDIVA is run in annotation mode only you will get the following annotated files:

+ <your-result-prefix>_snp_annotated.vcf
+ <your-result-prefix>_indel_annotated.vcf
+ <your-result-prefix>_indelExpanded_annotated.vcf
