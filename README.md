# aiDIVA - Augmented Intelligence Disease Variant Analysis

aiDIVA is an analysis pipeline that predicts the pathogenicity of given variants and prioritizes them according to the predicted pathogenicity to identify potential disease causing variants of a given sample.

The pathogenicity prediction is based on two random forest (RF) models. One is covering SNP variants, whereas the other one covers inframe InDel variants (frameshift variants are not covered).


## System requirements
The program is written in Python 3 (>= v3.6.9)

The following additional libraries need to be installed in order to use the program:

+ networkx (>= v1.11)
+ numpy (>= v1.13.3)
+ pandas (>= v0.22.0)
+ pysam (>= v0.14)
+ pyyaml (>= v3.12)
+ scipy (>= v0.19.1)
+ scikit-learn (>= v0.19.1)

If a newer scikit-learn version is used the models provided should still work (they were created using the version v0.19.1 and the import was tested with v0.21.3 and v0.22.2).


## Annotation resources and tools
If the annotation is not made beforehand make sure that the necessary database resources and tools are present on your system and the paths in the configuration file are set correctly.

For detailed instructions on how to download and prepare the annotation resources please head over to the respective [documentation](https://github.com/imgag/AIdiva/blob/master/doc/prepare_annotation_resources.md) (prepare_annotation_resources.md) found in the _doc_ folder.

For detailed instructions on how to download and install the annotation tools please head over to the respective [documentation](https://github.com/imgag/AIdiva/blob/master/doc/install_additional_tools.md) (isntall_additional_tools.md) found in the _doc_ folder.


## HPO resources
The HPO resources needed for the prioritization step can be found in the `data` folder. The path to the files is specified in the configuration file make sure that it leads to the correct location.For compatibility reasons the HPO graph resources were generated using networkx v1. Version 2 can still import these resources, but the graph generated with version 2 is not compatible with version 1.

To recreate the HPO resources please head over to the detailed [instructions](https://github.com/imgag/AIdiva/blob/master/doc/recreate_hpo_resources.md) found in the _doc_ folder.
<br>
<br>

Protein interactions based on String-DB v11.0

Last update of HPO resources: 30th July, 2021



## Pathogenicity prediction
There are two random forest models that are used in aiDIVA to predict the pathogenicity of a given variant. One for SNP variants and the other for inframe InDel variants. The training data of the two models consists of variants from Clinvar combined with additional variants from HGMD that are not present in Clinvar.

The scripts used to train the models can be found in the following GitHub repository: [aiDIVA-Training](https://github.com/imgag/aiDIVA-Training)

_Frameshift_ variants will get the no score, whereas _synonymous_ variants always get the lowest score 0.0

Pretrained random forest models using our current feature set can be found [here](https://download.imgag.de/ahboced1/AIdiva_pretrained_models/). The models were trained using scikit-learn v0.19.1. The trained models of scikit-learn are version dependent, but during our tests it also worked to load the 0.19.1 model with newer versions of scikit-learn (only the other way round it didn't work).

## Running aiDIVA
aiDIVA can be run either on already annotated VCF files or unannotated VCF files. In both cases a configuration file in the YAML format is required. When the variant annotation with VEP should also be performed with aiDIVA this is the only required command line argument. In the other case when an already annotated file is given there are a few more arguments that needs to be passed instead of being specified in the configuration file. The reason for this different set of parameters is due to the fact that it makes it more convenient to include aiDIVA in another existing pipeline.

Make sure to put the trained models in the data folder and make sure that the filename in the `aiDIVA_configuration_annotated.yaml` is correct.

### Running aiDIVA on already annotated data:

```
python run_AIdiva.py --config AIdiva_configuration_annotated.yaml --snp_vcf annotated_snp.vcf --indel_vcf annotated_indel.vcf --expanded_indel_vcf annotated_expanded_indel.vcf --out_prefix aidiva_result --workdir aidiva_workdir/ [--hpo_list hpo_terms.txt] [--gene_exclusion gene_exclusion.txt] [--family_file family.txt] [--family_type SINGLE] [--skip_db_check] [--only_top_results] [--threads 1] [--log_level INFO]
```

+ _config_ -- YAML configuration file (in the `data` folder there are example configuration files for each of the two modes)
+ _snp_vcf_ -- Input VCF file containing all annotated SNP variants of the given patient
+ _indel_vcf_ -- Input VCF file containing all InDel variants with basic annotation of the given patient
+ _expanded_indel_vcf_ -- Input VCF file containing all annotated expanded InDel of the given patient
+ _out_prefix_ -- A prefix for the resulting output files
+ _workdir_ -- Working directory, where all temporary files are created and saved (the results will also be stored here)
+ _hpo_list_ -- TXT file containing all the HPO terms observed with the patient [optional]
+ _gene_exclusion_ -- TXT file containing genes that should be excluded during the analysis of the HPO relatedness [optional]
+ _family_file_ -- TXT file containing the sample information if run on multisample VCF files [optional]
+ _family_type_ -- Type of the family relation [SINGLE, TRIO, FAMILY] (default: SINGLE) [optional]
+ _skip_db_check_ -- Skip the database checkup for existing entries in ClinVar (and HGMD) [optional]
+ _only_top_results_ -- Restrict the results to only report the top 25 variants [optional]
+ _threads_ -- Number of threads that should be used (default: 1) [optional]
+ _log_level_ -- Define logging level [DEBUG, INFO, WARN, ERROR, CRITICAL] (default: INFO) [optional]

### Running aiDIVA and perform the annotation:

```
python run_annotation_and_AIdiva.py --config AIdiva_configuration_with_annotation.yaml --vcf input.vcf --workdir aidiva_workdir/ [--hpo_list hpo_terms.txt] [--gene_exclusion gene_exclusion.txt] [--family_file family.txt] [--family_type SINGLE] [--skip_db_check] [--only_top_results] [--threads 1] [--log_level INFO]
```


## aiDIVA results
aiDIVA will produce three different output files two CSV files with the annotation and the results from the prediction and the prioritization. One of these two contains all variants, whereas the other one only contains the variants that passed the internal filtering step. The third file is a VCF file with the results in the INFO field, there are nine different fields in the INFO field. Four indicate a possible inheritance if the given VCF was a multisample VCF (These are missing if only a single sample VCF was given). One indicates if all internal filters were passed (AIDIVA_FILTER).
