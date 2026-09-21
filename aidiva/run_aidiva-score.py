import argparse
import logging
import pandas as pd
import tempfile
import variant_scoring.score_variants as predict
import yaml
import time


if __name__=="__main__":
    parser = argparse.ArgumentParser(description = "aiDIVA-Score")
    parser.add_argument("--config", type=str, dest="config", metavar="config.yaml", required=True, help="Config file specifying the parameters for aiDIVA [required]")
    parser.add_argument("--in_data", type=str, dest="in_data", metavar="input.tsv", required=True, help="TSV file with the annotated variants [required]")
    parser.add_argument("--out_prefix", type=str, dest="out_prefix", metavar="output_path/aidiva_result", required=True, help="Prefix that is used to save the results [required]")
    parser.add_argument("--workdir", type=str, dest="workdir", metavar="/tmp/aidiva_workdir/", required=False, help="Path to the working directory, here all intermediate files are saved (if not specified a temporary folder will be created and used)")
    parser.add_argument("--threads", type=int, dest="threads", metavar="1", required=False, help="Number of threads to use (default: 1)")
    parser.add_argument("--log_file", type=str, dest="log_file", metavar="output_path/logs/aidiva_log.txt", required=False, help="Path plus name of the log file to be saved, if not specified the log file is saved in the working directory")
    parser.add_argument("--log_level", type=str, dest="log_level", metavar="INFO", required=False, help="Define logging level, if unsure just leave the default [DEBUG, INFO] (default: INFO)")
    args = parser.parse_args()

    # parse input files
    if (args.in_data is not None):
        input_table = args.in_data

    else:
        raise SystemExit("The input table was not specified!")

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
        log_file = str(working_directory + "/" + "aidiva-score_" + timestamp + ".txt")

    # set up logger
    logging.basicConfig(filename=log_file,
                            filemode="a",
                            format=log_format,
                            datefmt="%H:%M:%S",
                            level=log_level)
    logger = logging.getLogger()

    logger.info("Running aiDIVA-score")
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

    # convert splitted input data to vcf and annotate
    if input_table is not None:
        ## TODO: Add possibility to add comments above the header line in the input table
        variant_table = pd.read_csv(input_table, sep="\t", low_memory=False)

    else:
        variant_table = pd.DataFrame()

    logger.debug(f"Condition-Check: {variant_table.dropna(how='all').empty}")
    logger.debug(f"Condition: {(not variant_table.dropna(how='all').empty)}")

    if (not variant_table.dropna(how='all').empty):
        # predict pathogenicity score
        logger.info("Score variants ...")
        variant_table_predicted = predict.perform_pathogenicity_score_prediction(scoring_model, variant_table, allele_frequency_list, feature_list, CONSTANT_DICTIONARY, num_cores)

        variant_table_predicted = variant_table_predicted.rename(columns={"CHROM": "#CHROM"})
        variant_table_predicted.to_csv(str(output_filename + "_result_aidiva-score.tsv"), sep="\t", index=False)

        logger.info("Pipeline successfully finished!")

    else:
        logger.warning("The given input file was empty!")
