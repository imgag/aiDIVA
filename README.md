# AIdiva - Augmented Intelligence Disease Variant Analysis

AIdiva is an analysis pipeline that predicts the pathogenicity of given variants and prioritizes them according to the predicted pathogenicity to identify potential disease causing variants of a given sample.

The pathogenicity prediction is based on two random forest (RF) models. One is covering SNP variants, whereas the other on covers inframe InDel variants (frameshift variants are always scored with the highest pathogenicity).

## System Requirements
The program is written in Python v3.7

The following additional libraries need to be installed in order to use the program:

+ argparse (v1.4.0)
+ biopython (v1.76)
+ networkx (v2.4)
+ numpy (v1.17.2)
+ pandas (v0.25.1)
+ pybigwig (v0.3.17)
+ pyyaml (v5.1.2)
+ scipy (v1.3.1)
+ scikit-learn (v0.21.3)


## Mandatory Commandline Arguments

To run the whole AIdiva pipeline one argument providing the configuration file is required.
This configuration file is encoded in the YAML format. In the `data` folder an example file can be found.


## Running AIdiva

The following commands show how to call the program from the commandline

`python run_AIdiva.py --config AIdiva_configuration.yaml`
