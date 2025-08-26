import argparse
import logging
import os
import numpy as np
import pandas as pd
import tempfile
import time
import meta_model.get_top_ranking_genes as top_ranking
import meta_model.llm_handling as llm_handler
import meta_model.metascore_handling as meta_handler
import yaml

#from mistralai import Mistral
from openai import OpenAI


if __name__=="__main__":
    parser = argparse.ArgumentParser(description = "aiDIVA -- augmented intelligence-based DIsease Variant Analysis")
    parser.add_argument("--in_rf", type=str, dest="in_rf", metavar="in_rf.tsv", required=True, help="TSV file containing the filtered results of aiDIVAs random forest prediction an ranking [required]")
    parser.add_argument("--in_eb_dom", type=str, dest="in_eb_dom", metavar="in_eb_dom.GSvar", required=False, help="GSvar file containing the evidence dominant based ranking [optional]")
    parser.add_argument("--in_eb_rec", type=str, dest="in_eb_rec", metavar="in_eb_rec.GSvar", required=False, help="GSvar file containing the evidence recessive based ranking [optional]")
    parser.add_argument("--out_prefix", type=str, dest="out_prefix", metavar="/output_path/aidiva_result", required=True, help="Prefix that is used to save the results [required]")
    parser.add_argument("--sample_id", type=str, dest="sample_id", metavar="NA12878_01", required=True, help="Sample ID that was used in previous annotations to store the genotype [required]")
    parser.add_argument("--config", type=str, dest="config", metavar="config.yaml", required=True, help="Config file specifying the parameters for aiDIVA [required]")
    parser.add_argument("--workdir", type=str, dest="workdir", metavar="/tmp/aidiva_workdir/", required=False, help="Path to the working directory, here all intermediate files are saved (if not specified a temporary folder will be created and used)")
    parser.add_argument("--hpo_list", type=str, dest="hpo_list", metavar="hpo.txt", required=False, help="TXT file containing the HPO terms reported for the current patient")
    parser.add_argument("--sex", type=str, dest="sex", metavar="male/female", required=False, help="Sex of the current patient")
    parser.add_argument("--age", type=str, dest="age", metavar="42", required=False, help="Age of the current patient")
    parser.add_argument("--evidence_based", dest="evidence_based", action="store_true", required=False, help="Flag to include evidence based rankings")
    parser.add_argument("--threads", type=int, dest="threads", metavar="1", required=False, help="Number of threads to use (default: 1)")
    parser.add_argument("--log_file", type=str, dest="log_file", metavar="/output_path/logs/aidiva_log.txt", required=False, help="Path plus name of the log file to be saved, if not specified the log file is saved in the working directory")
    parser.add_argument("--log_level", type=str, dest="log_level", metavar="INFO", required=False, help="Define logging level, if unsure just leave the default [DEBUG, INFO] (default: INFO)")
    args = parser.parse_args()

    # parse input files
    if (args.in_rf is not None):
        in_rf = args.in_rf
        result_data_random_forest = pd.read_csv(in_rf, sep="\t", low_memory=False)

    else:
        raise SystemExit("ERROR: The random forest based input TSV file was not specified!")

    evidence_based = args.evidence_based

    if evidence_based:
        if (args.in_eb_dom is not None) and (args.in_eb_rec is not None):
            in_eb_dom = args.in_eb_dom
            in_eb_rec = args.in_eb_rec

        else:
            raise SystemExit("ERROR: The evidence based input TSV files were not specified!")

    # parse configuration file
    if args.config is not None:
        try:
            config_file = open(args.config, "r")

        except OSError:
            raise SystemExit("ERROR: The given config file could not be opened!")

        else:
            configuration = yaml.load(config_file, Loader=yaml.SafeLoader)
            config_file.close()

    else:
        raise SystemExit("ERROR: The parameter 'config' was not specified!")

    # parse output files
    if args.out_prefix is not None:
        output_filename = args.out_prefix

    else:
        raise SystemExit("ERROR: The parameter 'out_prefix' was not specified!")

    # parse output files
    if args.sample_id is not None:
        sample_id = args.sample_id

    else:
        raise SystemExit("ERROR: The parameter 'sample_id' was not specified!")

    if args.workdir is not None:
        working_directory = args.workdir

        if not working_directory.endswith("/"):
            working_directory = working_directory + "/"

    else:
        workdir = tempfile.TemporaryDirectory(suffix="", prefix="aidiva_workdir_", dir=None)
        working_directory = workdir.name

        if not working_directory.endswith("/"):
            working_directory = working_directory + "/"

    # parse disease and patient information
    if args.hpo_list is not None:
        hpo_terms = args.hpo_list

    else:
        hpo_terms = None

    if args.sex is not None:
        sex = args.sex

    else:
        sex = ""

    if args.age is not None:
        age = args.age

    else:
        age = np.nan

    if args.threads is not None:
        num_cores = int(args.threads)

    else:
        num_cores = 1

    # use log level INFO as default
    if args.log_level is not None:
        if args.log_level == "DEBUG":
            log_level = logging.DEBUG
            log_format = "%(asctime)s -- %(name)s - %(levelname)s - %(message)s"

        elif args.log_level == "INFO":
            log_level = logging.INFO
            log_format = "%(asctime)s -- %(levelname)s - %(message)s"

        elif args.log_level == "WARNING":
            log_level = logging.WARNING
            log_format = "%(asctime)s -- %(levelname)s - %(message)s"

        elif args.log_level == "ERROR":
            log_level = logging.ERROR
            log_format = "%(asctime)s -- %(levelname)s - %(message)s"

        elif args.log_level == "CRITICAL":
            log_level = logging.CRITICAL
            log_format = "%(asctime)s -- %(levelname)s - %(message)s"

        else:
            log_level = logging.INFO
            log_format = "%(asctime)s -- %(levelname)s - %(message)s"

    else:
        log_level = logging.INFO
        log_format = "%(asctime)s -- %(levelname)s - %(message)s"

    if args.log_file is not None:
        log_file = args.log_file

    else:
        timestamp = time.strftime("%Y%m%d-%H%M%S")
        log_file = str(working_directory + "/" + "aidiva_" + timestamp + ".txt")

    # set up logger
    logging.basicConfig(filename=log_file,
                        filemode="a",
                        format=log_format,
                        datefmt="%H:%M:%S",
                        level=log_level)
    logger = logging.getLogger()

    logger.info("Running aiDIVAs meta model on predicted data")
    logger.info("Start program")
    logger.info(f"Working directory: {working_directory}")

    # load internal parameters
    internal_parameter_dict = configuration["Internal-Parameters"]

    # load meta ML model
    if evidence_based:
        meta_model = configuration["LLM-Input"]["meta-model"]

    else:
        meta_model = configuration["LLM-Input"]["meta-model-rf"]

    # specify LLM API to use
    llm_api = configuration["LLM-Input"]["llm-api"]

    # specify URL and port for locally hosted LLM
    if llm_api == "LOCAL":
        llm_api_url = configuration["LLM-Input"]["llm-api-url"]
        llm_api_port = configuration["LLM-Input"]["llm-api-port"]

    # load LLM API keyfile
    key_file = configuration["LLM-Input"]["llm-api-key"]

    # load LLM API keyfile
    llm_model = configuration["LLM-Input"]["llm-model"]

    # load LLM system instructions
    llm_instructions = configuration["LLM-Input"]["llm-instruction"]

    if (not result_data_random_forest.dropna(how='all').empty):
        # makes sure that the AIDIVA_RANK column exists, otherwise it is created (in older versions that column did not exist)
        if "AIDIVA_RANK" not in result_data_random_forest.columns:
            result_data_random_forest["AIDIVA_RANK"] = result_data_random_forest["FINAL_AIDIVA_SCORE"].rank(method='min', ascending=False)

        # extract top 10 gene list from random forest ranking
        top_ranking_genes_random_forest = top_ranking.extract_top_ranking_entries_random_forest_based(sample_id, result_data_random_forest, 10)

        # extract top10 list from evidence ranking (dominant and recessive)
        if evidence_based:
            top_ranking_genes_evidence_dominant = top_ranking.extract_top_ranking_entries_evidence_based(sample_id, in_eb_dom, 10) #result_data_evidence_dominant, 10)
            top_ranking_genes_evidence_recessive = top_ranking.extract_top_ranking_entries_evidence_based(sample_id, in_eb_rec, 10) #result_data_evidence_recessive, 10)

        if llm_api == "LOCAL":
            ## use this client definition to use a local LLM (make sure that the service serving the model is accessible)
            client = OpenAI(base_url=f"{llm_api_url}:{llm_api_port}/v1/", api_key="UNKNOWN")

        elif llm_api == "OPENAI":
            # Get your OpenAI API key
            with open(key_file, "r") as keyfile:
                llm_api_key = keyfile.readline().rstrip()

            # Initialize the OpenAI API client
            client = OpenAI(api_key=llm_api_key)

        #elif llm_api == "MISTRALAI"
        #    # Get your MistralAI API key
        #    with open(llm_api_key_file, "r") as keyfile:
        #        llm_api_key = keyfile.readline().rstrip()

        #    # Initialize the MistralAI API client
        #    client = Mistral(api_key=llm_api_key)

        else:
            raise SystemExit("You need to specify a valid LLM api (OPENAI, LOCAL)!")

        llm_prompt_random_forest = llm_handler.create_llm_prompt(sex, age, hpo_terms, top_ranking_genes_random_forest, "rf", internal_parameter_dict)
        llm_refined_results_random_forest = llm_handler.call_llm_api(client, llm_prompt_random_forest, llm_instructions, llm_model, llm_api)

        llm_refined_results_random_forest.to_csv(str(output_filename + "_aidiva_random_forest_based_llm_results.tsv"), sep="\t", index=False)

        if evidence_based:
            llm_prompt_evidence_dominant = llm_handler.create_llm_prompt(sex, age, hpo_terms, top_ranking_genes_evidence_dominant, "eb_dom", internal_parameter_dict)
            llm_prompt_evidence_recessive = llm_handler.create_llm_prompt(sex, age, hpo_terms, top_ranking_genes_evidence_recessive, "eb_rec", internal_parameter_dict)

            llm_refined_results_evidence_dominant = llm_handler.call_llm_api(client, llm_prompt_evidence_dominant, llm_instructions, llm_model, llm_api)
            llm_refined_results_evidence_recessive = llm_handler.call_llm_api(client, llm_prompt_evidence_recessive, llm_instructions, llm_model, llm_api)

            llm_refined_results_evidence_dominant.to_csv(str(output_filename + "_aidiva_evidence_dominant_based_llm_results.tsv"), sep="\t", index=False)
            llm_refined_results_evidence_recessive.to_csv(str(output_filename + "_aidiva_evidence_recessive_based_llm_results.tsv"), sep="\t", index=False)

        # create metascore table
        metascore_table_random_forest = meta_handler.create_table_rf_based(llm_refined_results_random_forest, top_ranking_genes_random_forest, False)

        if evidence_based:
            metascore_table_random_forest_and_evidence = meta_handler.create_table_rf_and_evidence_based(llm_refined_results_random_forest, top_ranking_genes_random_forest, llm_refined_results_evidence_dominant, top_ranking_genes_evidence_dominant, llm_refined_results_evidence_recessive, top_ranking_genes_evidence_recessive, False)

        # use metascore model to predict new values and create final ranking
        if evidence_based:
            metascore_table_random_forest_and_evidence_predicted = meta_handler.meta_scoring(metascore_table_random_forest_and_evidence, meta_model, False)

        else:
           metascore_table_random_forest_and_evidence_predicted = meta_handler.meta_scoring(metascore_table_random_forest, meta_model, True)

        # extract final ranking and create result file
        metascore_table_random_forest_and_evidence_predicted.to_csv(str(output_filename + "_aidiva_metascore_results.tsv"), sep="\t", index=False)

        logger.info("Pipeline successfully finsished!")

    else:
        logger.warning("The given input files were empty!")
