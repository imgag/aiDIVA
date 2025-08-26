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


def perform_llm_refinement(row, evidence_based, llm_api, llm_api_key, llm_api_url, llm_api_port, llm_model, llm_instructions, meta_ml_model):
    sex = row["sex"]
    age = row["age"]
    hpo_terms = row["hpo"]
    causal_gene = row["causal_gene"]

    top_ranking_genes_random_forest = row["top_10_rf"]

    if evidence_based:
        top_ranking_genes_evidence_dominant = row["top_10_eb_dom"]
        top_ranking_genes_evidence_recessive = row["top_10_eb_rec"]

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

    if evidence_based:
        llm_prompt_evidence_dominant = llm_handler.create_llm_prompt(sex, age, hpo_terms, top_ranking_genes_evidence_dominant, "eb_dom", internal_parameter_dict)
        llm_prompt_evidence_recessive = llm_handler.create_llm_prompt(sex, age, hpo_terms, top_ranking_genes_evidence_recessive, "eb_rec", internal_parameter_dict)

        llm_refined_results_evidence_dominant = llm_handler.call_llm_api(client, llm_prompt_evidence_dominant, llm_instructions, llm_model, llm_api)
        llm_refined_results_evidence_recessive = llm_handler.call_llm_api(client, llm_prompt_evidence_recessive, llm_instructions, llm_model, llm_api)

        print(llm_prompt_evidence_dominant, "\n\n")
        print(llm_refined_results_evidence_dominant, "\n\n")
        print(llm_prompt_evidence_recessive, "\n\n")
        print(llm_refined_results_evidence_recessive, "\n\n")

    # create metascore table
    metascore_table_random_forest = meta_handler.create_table_rf_based(llm_refined_results_random_forest, top_ranking_genes_random_forest, True)

    if evidence_based:
        metascore_table_random_forest_and_evidence = meta_handler.create_table_rf_and_evidence_based(llm_refined_results_random_forest, top_ranking_genes_random_forest, llm_refined_results_evidence_dominant, top_ranking_genes_evidence_dominant, llm_refined_results_evidence_recessive, top_ranking_genes_evidence_recessive, True)

    # use metascore model to predict new values and create final ranking
    if evidence_based:
        metascore_table_random_forest_and_evidence_predicted = meta_handler.meta_scoring(metascore_table_random_forest_and_evidence, meta_ml_model, False)

    else:
       metascore_table_random_forest_and_evidence_predicted = meta_handler.meta_scoring(metascore_table_random_forest, meta_ml_model, True)

    llm_rank_rf = -1
    llm_rank_eb = -1
    rank_metascore = -1
    metascore = np.nan

    for entry in metascore_table_random_forest_and_evidence_predicted.itertuples():
        if entry.gene_name in causal_gene:
            if entry.rf_rank_llm < 4:
                llm_rank_rf = entry.rf_rank_llm

            if evidence_based:
                if entry.eb_rank_llm < 4:
                    llm_rank_eb = entry.eb_rank_llm

            rank_metascore = entry.meta_rank
            metascore = entry.metascore_prediction

            break

    if evidence_based:
        return llm_rank_rf, llm_rank_eb, rank_metascore, metascore

    else:
        return llm_rank_rf, rank_metascore, metascore


if __name__=="__main__":
    parser = argparse.ArgumentParser(description = "aiDIVA -- augmented intelligence-based DIsease Variant Analysis")
    parser.add_argument("--in_data", type=str, dest="in_data", metavar="in_data.tsv", required=True, help="TSV file containing the filtered results of aiDIVAs random forest prediction and ranking as well as the evidence based results [required]")
    parser.add_argument("--out_data", type=str, dest="out_data", metavar="out_data.tsv", required=True, help="Prefix that is used to save the results [required]")
    parser.add_argument("--config", type=str, dest="config", metavar="config.yaml", required=True, help="Config file specifying the parameters for aiDIVA [required]")
    parser.add_argument("--evidence_based", dest="evidence_based", action="store_true", required=False, help="Flag to include evidence based rankings")
    parser.add_argument("--log_file", type=str, dest="log_file", metavar="/output_path/logs/aidiva_log.txt", required=False, help="Path plus name of the log file to be saved, if not specified the log file is saved in the working directory")
    parser.add_argument("--log_level", type=str, dest="log_level", metavar="INFO", required=False, help="Define logging level, if unsure just leave the default [DEBUG, INFO] (default: INFO)")
    args = parser.parse_args()

    # parse input files
    if (args.in_data is not None):
        in_data = args.in_data
        result_data = pd.read_csv(in_data, sep="\t", low_memory=False)

    else:
        raise SystemExit("ERROR: The random forest based input TSV file was not specified!")

    evidence_based = args.evidence_based

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
    if args.out_data is not None:
        output_filename = args.out_data

    else:
        raise SystemExit("ERROR: The parameter 'out_prefix' was not specified!")

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

    workdir = tempfile.TemporaryDirectory(suffix="", prefix="aidiva_workdir_", dir=None)
    working_directory = workdir.name

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

    logger.info("Running aiDIVAs meta model on anonymized benchmark data")
    logger.info("Start program")

    # load internal parameters
    internal_parameter_dict = configuration["Internal-Parameters"]

    # load meta ML model
    if evidence_based:
        meta_ml_model = configuration["LLM-Input"]["meta-model"]

    else:
        meta_ml_model = configuration["LLM-Input"]["meta-model-rf"]

    # specify LLM API to use
    llm_api = configuration["LLM-Input"]["llm-api"]

    # specify URL and port for locally hosted LLM
    if llm_api == "LOCAL":
        llm_api_url = configuration["LLM-Input"]["llm-api-url"]
        llm_api_port = configuration["LLM-Input"]["llm-api-port"]

    # load LLM API keyfile
    llm_api_key = configuration["LLM-Input"]["llm-api-key"]

    # load LLM API keyfile
    llm_model = configuration["LLM-Input"]["llm-model"]

    # load LLM system instructions
    llm_instructions = configuration["LLM-Input"]["llm-instruction"]

    if (not result_data.dropna(how='all').empty):
        # add anonymous patient id to handle intermediate result files (e.g. patient_01) if id is missing (this can also be just the row index)
        # afterwards use apply to create and add all the results
        #(anonymized_id) | causal gene | (causal variant) | age | sex | hpo | transl. hpo | top 10 rf | top 10 eb dom | top 10 eb rec | rf rank | rf score | eb dom rank | eb dom score | eb rec rank | eb rec score | exomiser rank | AIM rank 
        # add directly the result for the causal gene
        if evidence_based:
            result_data[["llm_rank_rf", "llm_rank_eb", "rank_metascore", "metascore"]] = result_data.apply(lambda row: pd.Series(perform_llm_refinement(row, True, llm_api, llm_api_key, llm_api_url, llm_api_port, llm_model, llm_instructions, meta_ml_model)), axis=1)

        else:
            result_data[["llm_rank_rf", "rank_metascore", "metascore"]] = result_data.apply(lambda row: pd.Series(perform_llm_refinement(row, False, llm_api, llm_api_key, llm_api_url, llm_api_port, llm_model, llm_instructions, meta_ml_model)), axis=1)

        # save updated table
        result_data.to_csv(output_filename, sep="\t", index=False)

        logger.info("Pipeline successfully finsished!")

    else:
        logger.warning("The given input files were empty!")
