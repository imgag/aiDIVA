# Annotate VCF file for aiDIVA

The `doc/annotation` folder contains the documentation for the following required steps to successfully run the annotation: 

1) Installation of all necessary software dependencies ([*install\_additional\_tools.md*](https://github.com/imgag/aiDIVA/blob/master/doc/annotation/install_additional_tools.md))

2) Preparation of annotation resources required for running the pipeline ([*prepare\_annotation\_resources.md*](https://github.com/imgag/aiDIVA/blob/master/doc/annotation/prepare_annotation_resources.md))

## Running Annotation:

```
python3 run_annotation.py --vcf input.vcf --config configuration_annotation.yaml --out_folder output_path/ [--filtered] [--filtered_folder output_path/aidiva_filtered/] [--inhouse_sample] [--compress] [--threads 1] [--log_file annotation_log.txt] [--log_level INFO]
```
mandatory parameters:

+ *vcf* -- Input VCF file containing all variants of the given sample
+ *config* -- YAML configuration file (in the `data` folder is an example configuration file)
+ *out\_folder* -- Folder where the resulting output file should be saved

optional parameters:

+ *filtered* -- Skip filtering step if it was already done \[optional\]
+ *filtered\_folder* -- Folder where to find the already filtered VCF files \[optional\]
+ *inhouse\_sample* -- The input VCF was prepared inhouse (skips some additional preparation steps to prevent possible problems) \[optional\]
+ *compress* -- If this flag is set the output table will be compressed with gzip \[optional\]
+ *threads* -- Number of threads that should be used (default: 1) \[optional\]
+ *log\_file* -- Specify a custom log file to store the log messages from the tool \[optional\]
+ *log\_level* -- Define logging level \[DEBUG, INFO\] (default: INFO) \[optional\]

## Annotation Result

After running the annotation script you will get the following annotated file that is ready to be used with the aiDIVA software:

+ *\<your-output-folder\>/\<input-filename\>\_annotated.tsv*

### Convert Annotated Table to GSvar
To convert the annotated table *\<input-filename\>\_annotated.tsv* to the GSvar format needed to run the VariantRanking tool you can use the `convert_annotated_table_to_gsvar.py` script found in the `scripts` folder.

```
python3 convert_annotated_table_to_gsvar.py --in_file NA12878_annotated.tsv --out_file NA12878_annotated.GSvar --sample_id NA12878
```
