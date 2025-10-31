import argparse
import logging
import os
import pandas as pd
import tempfile
import time
import variant_scoring.score_variants as predict
import variant_prioritization.prioritize_variants as prio
import yaml
import time


if __name__=="__main__":
    parser = argparse.ArgumentParser(description = "aiDIVA-RF")
    parser.add_argument("--config", type=str, dest="config", metavar="config.yaml", required=True, help="Config file specifying the parameters for aiDIVA [required]")
    parser.add_argument("--in_data", type=str, dest="in_data", metavar="input.tsv", required=True, help="TSV file with the annotated variants [required]")
    parser.add_argument("--out_prefix", type=str, dest="out_prefix", metavar="output_path/aidiva_result", required=True, help="Prefix that is used to save the results [required]")
    parser.add_argument("--workdir", type=str, dest="workdir", metavar="/tmp/aidiva_workdir/", required=False, help="Path to the working directory, here all intermediate files are saved (if not specified a temporary folder will be created and used)")
    parser.add_argument("--hpo_list", type=str, dest="hpo_list", metavar="hpo.txt", required=False, help="TXT file containing the HPO terms reported for the current patient")
    parser.add_argument("--gene_exclusion", type=str, dest="gene_exclusion", metavar="gene_exclusion.txt", required=False, help="Tab separated file containing the genes to exclude in the analysis. Genes are assumed to be in the first column.")
    parser.add_argument("--family_file", type=str, dest="family_file", metavar="family.txt", required=False, help="TXT file showing the sample relations of the current data")
    parser.add_argument("--family_type", type=str, dest="family_type", metavar="SINGLE", required=False, help="In case of multisample data the kind of sample relation [SINGLE, TRIO, MULTI]")
    parser.add_argument("--skip_db_check", dest="skip_db_check", action="store_true", required=False, help="Flag to skip database (ClinVar, HGMD) lookup")
    parser.add_argument("--rare_disease", dest="rare_disease", action="store_true", required=False, help="Set rare disease mode: Activate initial filtering to remove all variants with a maximum allele frequency of more than 2%. This setting is meant to speed up the software in a rare disease setting.")
    parser.add_argument("--only_top_results", dest="only_top_results", action="store_true", required=False, help="Report only the top ranking variants as result. The desired rank can be given as parameter with '--top_rank' (default: 25)")
    parser.add_argument("--top_rank", type=str, dest="top_rank", metavar="25", required=False, help="Rank parameter for '--only_top_results' (default: 25)")
    parser.add_argument("--threads", type=int, dest="threads", metavar="1", required=False, help="Number of threads to use (default: 1)")
    parser.add_argument("--log_file", type=str, dest="log_file", metavar="output_path/logs/aidiva_log.txt", required=False, help="Path plus name of the log file to be saved, if not specified the log file is saved in the working directory")
    parser.add_argument("--log_level", type=str, dest="log_level", metavar="INFO", required=False, help="Define logging level, if unsure just leave the default [DEBUG, INFO] (default: INFO)")
    args = parser.parse_args()

    # parse input files
    if (args.in_data is not None):
        input_table = args.in_data

    else:
        raise SystemExit("The input VCF files were not specified!")

    # parse configuration file
    if args.config is not None:
        try:
            config_file = open(args.config, "r")

        except OSError:
            raise SystemExit("The given config file could not be opened!")

        else:
            configuration = yaml.load(config_file, Loader=yaml.SafeLoader)
            config_file.close()

    else:
        raise SystemExit("The parameter 'config' was not specified!")

    # parse output files
    if args.out_prefix is not None:
        output_filename = args.out_prefix

    else:
        raise SystemExit("The parameter 'out_prefix' was not specified!")

    if args.workdir is not None:
        working_directory = args.workdir

        if not working_directory.endswith("/"):
            working_directory = working_directory + "/"

    else:
        workdir = tempfile.TemporaryDirectory(suffix="", prefix="aidiva_workdir_", dir=None)
        working_directory = workdir.name

        if not working_directory.endswith("/"):
            working_directory = working_directory + "/"

    # parse disease and inheritance information
    if args.hpo_list is not None:
        hpo_file = args.hpo_list

    else:
        hpo_file = None

    if args.gene_exclusion is not None:
        gene_exclusion_file = args.gene_exclusion

    else:
        gene_exclusion_file = None

    if (args.family_file is not None) and (args.family_type is not None):
        family_file = args.family_file
        family_type = args.family_type

    else:
        family_file = None
        family_type = "SINGLE"

    rare_disease = args.rare_disease
    skip_db_check = args.skip_db_check
    only_top_results = args.only_top_results

    if args.threads is not None:
        num_cores = int(args.threads)

    else:
        num_cores = 1

    if args.top_rank is not None:
        top_rank = int(args.top_rank)

    else:
        top_rank = 25

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
        log_file = str(working_directory + "/" + "aidiva-rf_" + timestamp + ".txt")

    # set up logger
    logging.basicConfig(filename=log_file,
                            filemode="a",
                            format=log_format,
                            datefmt="%H:%M:%S",
                            level=log_level)
    logger = logging.getLogger()

    logger.info("Running aiDIVA-RF")
    logger.info("Start program")
    logger.info(f"Working directory: {working_directory}")

    # load ML model
    scoring_model = configuration["Analysis-Input"]["scoring-model"]

    # load internal parameters
    internal_parameter_dict = configuration["Internal-Parameters"]
    CONSTANT_DICTIONARY = configuration["Internal-Parameters"]["CONSTANTS"]

    allele_frequency_list = configuration["Model-Features"]["allele-frequency-list"]
    feature_list = configuration["Model-Features"]["feature-list"]
    assembly_build = configuration["Assembly-Build"]
    ref_path = configuration["Reference-Genome"]

    #transcript_path = configuration["Canonical-Transcripts"]
    prioritization_weights = configuration["Analysis-Input"]["prioritization-weights"]
    filter_identifiers = configuration["Analysis-Input"]["prioritization-information"]

    # convert splitted input data to vcf and annotate
    if input_table is not None:
        ## TODO: Add possibility to add comments above the header line in the input table
        variant_table = pd.read_csv(input_table, sep="\t", low_memory=False)

        # prefilter input table for rare disease mode
        if rare_disease:
            logger.info("Rare disease mode activated! Filter out all variants with allele frequence higher than 2%")
            if "MAX_AF" in variant_table.columns:
                variant_table = variant_table[variant_table["MAX_AF"] <= 0.02]

            else:
                if allele_frequency_list:
                    for allele_frequency in allele_frequency_list:
                        variant_table[allele_frequency] = variant_table[allele_frequency].fillna(0)
                        variant_table[allele_frequency] = variant_table.apply(lambda row: pd.Series(max([float(frequency) for frequency in str(row[allele_frequency]).split("&")], default=np.nan)), axis=1)

                    variant_table["MAX_AF"] = variant_table.apply(lambda row: pd.Series(max([float(frequency) for frequency in row[allele_frequency_list].tolist()], default=np.nan)), axis=1)

                else:
                    raise SystemExit("Could not identify the allele frequency information in the input table!")

    else:
        variant_table = pd.DataFrame()

    logger.debug(f"Condition-Check: {variant_table.dropna(how='all').empty}")
    logger.debug(f"Condition: {(not variant_table.dropna(how='all').empty)}")

    if (not variant_table.dropna(how='all').empty):
        # predict pathogenicity score
        logger.info("Score variants ...")
        variant_table_predicted = predict.perform_pathogenicity_score_prediction(scoring_model, variant_table, allele_frequency_list, feature_list, CONSTANT_DICTIONARY, num_cores)

        # prioritize and filter variants
        logger.info("Prioritize variants and finalize score ...")
        variant_table_prioritized = prio.prioritize_variants(variant_table_predicted, internal_parameter_dict, prioritization_weights, filter_identifiers, ref_path, num_cores, assembly_build, feature_list, CONSTANT_DICTIONARY, skip_db_check, family_file, family_type, hpo_file, gene_exclusion_file)

        variant_table_prioritized_filtered = variant_table_prioritized[variant_table_prioritized["FILTER_PASSED"] == 1].copy(deep=True)
        variant_table_prioritized_filtered["AIDIVA_RANK"] = variant_table_prioritized_filtered["FINAL_AIDIVA_SCORE"].rank(method='min', ascending=False)

        if only_top_results:
            variant_table_prioritized_filtered[variant_table_prioritized_filtered["AIDIVA_RANK"] <= top_rank].to_csv(str(output_filename + "_result_filtered_aidiva-rf.tsv"), sep="\t", index=False)
            logger.info(f"Only top ranking variants are reported as result (rank: {top_rank})!")

        else:
            variant_table_prioritized = variant_table_prioritized.rename(columns={"CHROM": "#CHROM"})
            variant_table_prioritized.to_csv(str(output_filename + "_result_aidiva-rf.tsv"), sep="\t", index=False)
            variant_table_prioritized_filtered.to_csv(str(output_filename + "_result_filtered_aidiva-rf.tsv"), sep="\t", index=False)

            ##TODO add additional condition for type FAMILY (this type is currently not supported)
            if family_type == "TRIO":
                variant_table_prioritized_filtered[(variant_table_prioritized_filtered["COMPOUND"] == 1)].to_csv(str(output_filename + "_result_filtered_compound_aidiva-rf.tsv"), sep="\t", index=False)
                variant_table_prioritized_filtered[(variant_table_prioritized_filtered["DOMINANT_DENOVO"] == 1)].to_csv(str(output_filename + "_result_filtered_dominant_denovo_aidiva-rf.tsv"), sep="\t", index=False)
                variant_table_prioritized_filtered[(variant_table_prioritized_filtered["DOMINANT"] == 1)].to_csv(str(output_filename + "_result_filtered_dominant_aidiva-rf.tsv"), sep="\t", index=False)
                variant_table_prioritized_filtered[(variant_table_prioritized_filtered["XLINKED"] == 1)].to_csv(str(output_filename + "_result_filtered_xlinked_aidiva-rf.tsv"), sep="\t", index=False)
                variant_table_prioritized_filtered[(variant_table_prioritized_filtered["RECESSIVE"] == 1)].to_csv(str(output_filename + "_result_filtered_recessive_aidiva-rf.tsv"), sep="\t", index=False)

        logger.info("Pipeline successfully finsished!")

    else:
        logger.warning("The given input file was empty!")
