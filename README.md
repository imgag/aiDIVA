# aiDIVA - augmented intelligence-based DIsease Variant Analysis

aiDIVA is an analysis pipeline that combines a pathogenicity-based approach and an optional evidence-based approach with state-of-the-art large language models (LLMs) to identify potential disease causing variants in a given rare disease sample.

aiDIVA comprises the following steps:

1. A pathogenicity-based approach that utilizes the predictions of a random forest model trained on ClinVar data, supplemented with phenotype information given as HPO terms.

2. An optional evidence-based approach includes the ranks and scores that can be obtained from the [VariantRanking](https://github.com/imgag/ngs-bits/blob/master/doc/tools/VariantRanking/index.md) tool included in [ngs-bits](https://github.com/imgag/ngs-bits/).

3. State-of-the-art LLMs are utilized to refine the ranking results from the previous two approaches.

4. A final meta model is used to combine all preliminary results and create a final ranking of the variants. For this model we also use a random forest.


## Citing

If you use aiDIVA in your work please cite our preprint: 

**aiDIVA - Diagnostics of Rare Genetic Diseases Using Large Language Models** ([link](https://www.medrxiv.org/content/10.1101/2025.09.04.25335099v1))

<!--
<br><br>

Additionally you may use the following Zenodo-Records to point to specific versions of the tool:

+ 1.0.0: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.16966353.svg)](https://doi.org/10.5281/zenodo.16966353)
+ 1.0.1: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17277130.svg)](https://doi.org/10.5281/zenodo.17277130)
-->


## Support

Please report any issues or questions to the [aiDIVA Issue Tracker](https://github.com/imgag/aiDIVA/issues).


## System Requirements

The program is written in Python 3

Latest used version: 3.12.3

The following additional libraries need to be installed in order to use the program (latest used version):

+ networkx (v3.4.2)
+ numpy (v1.26.4)
+ openai (v1.60.2)
+ pandas (v2.2.3)
+ pysam (v0.22.1)
+ pyyaml (v6.0.2)
+ scipy (v1.15.1)
+ scikit-learn (v1.3.2)

For easy package installation in your Python virtual environment we included a requirements.txt just run `pip install -r requirements.txt` to install all necessary packages (the versions in the *requirements.txt* match our own setup at that time).

If a newer scikit-learn version is used it is advised to create a new model with the newer scikit-learn version.


## Run the aiDIVA software

To run the aiDIVA software you need a TAB separated file containing the annotation information for every variant present in your sample file.

Detailed instructions on how to run the software and what columns need to be present in the input table can be found [here](https://github.com/imgag/aiDIVA/blob/master/doc/aidiva/run_aidiva.md).


## Annotation Resources and Tools

If you don't have an annotated table with the necessary columns mentioned before. You can use the _run_annotation_ script provided in the _annotation_ folder to create a table with the necessary information.
Before you use this annotation script make sure that the necessary database resources and tools are present on your system and the paths in the configuration file are set correctly (IMPORTANT: use the correct configuration file, it differs between the annotation and aiDIVA!).

Instructions on how to use the annotation script and prepare the annotation resources and tools can be found [here](https://github.com/imgag/aiDIVA/blob/master/doc/annotation/run_annotation.md).


## HPO Resources

The HPO resources needed for the prioritization step can be found in the `data` folder. The path to the files is specified in the configuration file make sure that it leads to the correct location. Although we used networkx v3.4 to create the HPO graph it should also be possible to use networkx in version 1 or 2. We included some workarounds to still support the older versions.

To recreate the HPO resources please head over to the detailed [instructions](https://github.com/imgag/aiDIVA/blob/master/doc/aidiva/recreate_hpo_resources.md) found in the `doc` folder.

Protein interactions based on String-DB v11.0b

Last update of HPO and other resources: 22nd May, 2023


## Pathogenicity Prediction

There is one random forest model that is used in aiDIVA to predict the pathogenicity of a given variant. It is a combined model for SNV and inframe indel variants. The training data of the model consists of variants from Clinvar.

The scripts used to train the model can be found in the following GitHub repository: [aiDIVA-Training](https://github.com/imgag/aiDIVA-Training)

*Frameshift* variants will get a default score of 0.9, whereas *synonymous* variants always get the lowest score 0.0

A pretrained random forest model (*aidiva-rf*) using our current feature set can be found [here](https://download.imgag.de/aidiva/aidiva_pretrained_models/). The latest model was trained using scikit-learn v1.3.2. The trained models of scikit-learn are version dependent.


## LLM Usage

aiDIVA supports the use of the official [OpenAI API](https://platform.openai.com/docs/api-reference/introduction) to send the requests to GPT-4o or GPT-4.1 for example. To use the OpenAI API you need an account and an API-Key that needs to be specified in the configuration file.
Alternatively it is possible to set up your own local LLM (eg., LLama-8b, Mistral-12b, ...) and provide it locally as a Webservice. For an easy deployment you could use the NVIDIA NIM Containers see [here](https://build.nvidia.com/meta/llama-3_1-8b-instruct/deploy) for more details on how to do that. These local LLMs use the same python package for inference you just have to specify the port and URL where to find the local model in the configuration file.


## Meta Model

You can download the pretrained meta models (*aidiva-meta* & *aidiva-meta-rf*) [here](https://download.imgag.de/aidiva/aidiva_pretrained_models/). For these two models we used a random forest model that takes as features the ranking position and scores from the initial rankings (pathogenicity-based and evidence-based) plus the ranking result from the LLMs and the inheritance mode used in the evidence-based model.

