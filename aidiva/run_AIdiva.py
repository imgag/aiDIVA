import argparse
import helper_modules.combine_expanded_indels_and_create_csv as combine_expanded_indels
import helper_modules.create_result_vcf as write_result
import helper_modules.convert_vcf_to_csv as convert_vcf
import logging
import os
import pandas as pd
import tempfile
import time
import variant_scoring.score_variants as predict
import variant_prioritization.prioritize_variants as prio
import yaml


if __name__=="__main__":
    parser = argparse.ArgumentParser(description = "AIdiva -- Augmented Intelligence Disease Variant Analysis")
    parser.add_argument("--snp_vcf", type=str, dest="snp_vcf", metavar="snp.vcf", required=True, help="VCF file with the annotated SNP variants [required]")
    parser.add_argument("--indel_vcf", type=str, dest="indel_vcf", metavar="indel.vcf", required=True, help="VCF file with the annotated (only basic annotation) InDel variants [required]")
    parser.add_argument("--expanded_indel_vcf", type=str, dest="expanded_indel_vcf", metavar="expanded_indel.vcf", required=True, help="VCF file with the annotated expanded InDel variants [required]")
    parser.add_argument("--out_prefix", type=str, dest="out_prefix", metavar="/output_path/aidiva_result", required=True, help="Prefix that is used to save the results [required]")
    parser.add_argument("--config", type=str, dest="config", metavar="config.yaml", required=True, help="Config file specifying the parameters for AIdiva [required]")
    parser.add_argument("--workdir", type=str, dest="workdir", metavar="/tmp/aidiva_workdir/", required=False, help="Path to the working directory, here all intermediate files are saved (if not specified a temporary folder will be created and used)")
    parser.add_argument("--hpo_list", type=str, dest="hpo_list", metavar="hpo.txt", required=False, help="TXT file containing the HPO terms reported for the current patient")
    parser.add_argument("--gene_exclusion", type=str, dest="gene_exclusion", metavar="gene_exclusion.txt", required=False, help="Tab separated file containing the genes to exclude in the analysis. Genes are assumed to be in the first column.")
    parser.add_argument("--family_file", type=str, dest="family_file", metavar="family.txt", required=False, help="TXT file showing the sample relations of the current data")
    parser.add_argument("--family_type", type=str, dest="family_type", metavar="SINGLE", required=False, help="In case of multisample data the kind of sample relation [SINGLE, TRIO, MULTI]")
    parser.add_argument("--skip_db_check", dest="skip_db_check", action="store_true", required=False, help="Flag to skip database (ClinVar, HGMD) lookup")
    parser.add_argument("--only_top_results", dest="only_top_results", action="store_true", required=False, help="Report only the top 25 variants as result")
    parser.add_argument("--threads", type=int, dest="threads", metavar="1", required=False, help="Number of threads to use (default: 1)")
    parser.add_argument("--log_path", type=str, dest="log_path", metavar="/output_path/logs/", required=False, help="Path where the log file should be saved, if not specified the log file is saved in the working directory")
    parser.add_argument("--log_level", type=str, dest="log_level", metavar="INFO", required=False, help="Define logging level, if unsure just leave the default [DEBUG, INFO] (default: INFO)")
    parser.add_argument("--save_as_vcf", dest="save_as_vcf", action="store_true", required=False, help="Flag to create additional result files in VCF format.")
    args = parser.parse_args()
    
    # parse input files
    if (args.snp_vcf is not None) and (args.indel_vcf is not None and args.expanded_indel_vcf is not None):
        snp_vcf = args.snp_vcf
        indel_vcf = args.indel_vcf
        expanded_indel_vcf = args.expanded_indel_vcf
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

    skip_db_check = args.skip_db_check
    only_top_results = args.only_top_results
    save_as_vcf = args.save_as_vcf
    
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
        #elif args.log_level == "WARN":
        #    log_level = logging.WARN
        #    log_format = "%(asctime)s -- %(levelname)s - %(message)s"
        #elif args.log_level == "ERROR":
        #    log_level = logging.ERROR
        #    log_format = "%(asctime)s -- %(levelname)s - %(message)s"
        #elif args.log_level == "CRITICAL":
        #    log_level = logging.CRITICAL
        #    log_format = "%(asctime)s -- %(levelname)s - %(message)s"
        else:
            log_level = logging.INFO
            log_format = "%(asctime)s -- %(levelname)s - %(message)s"
    else:
        log_level = logging.INFO
        log_format = "%(asctime)s -- %(levelname)s - %(message)s"

    if args.log_path is not None:
        log_path = args.log_path
    else:
        log_path = working_directory

    # set up logger
    timestamp = time.strftime("%Y%m%d-%H%M%S")
    logging.basicConfig(filename=str(log_path + "/" + "aidiva_" + timestamp + ".log"),
                            filemode="a",
                            format=log_format,
                            datefmt="%H:%M:%S",
                            level=log_level)
    logger = logging.getLogger()

    logger.info("Running AIdiva on annotated data")
    logger.info("Start program")
    logger.info(f"Working directory: {working_directory}")

    # load SNP ML model
    scoring_model_snp = configuration["Analysis-Input"]["scoring-model-snp"]
    
    # load InDel ML model
    scoring_model_indel = configuration["Analysis-Input"]["scoring-model-indel"]

    # load internal parameters
    internal_parameter_dict = configuration["Internal-Parameters"]

    allele_frequency_list = configuration["Model-Features"]["allele-frequency-list"]
    feature_list = configuration["Model-Features"]["feature-list"]
    assembly_build = configuration["Assembly-Build"]
    ref_path = configuration["Analysis-Input"]["ref-path"]

    # convert splitted input data to vcf and annotate
    if snp_vcf is not None:
        input_data_snp = convert_vcf.convert_vcf_to_pandas_dataframe(snp_vcf, False, False, num_cores)
    else:
        input_data_snp = pd.DataFrame()

    if indel_vcf is not None and expanded_indel_vcf is not None:
        input_data_indel = convert_vcf.convert_vcf_to_pandas_dataframe(indel_vcf, True, False, num_cores)
        input_data_expanded_indel = convert_vcf.convert_vcf_to_pandas_dataframe(expanded_indel_vcf, True, True, num_cores)
    else:
        input_data_indel = pd.DataFrame()
        input_data_expanded_indel = pd.DataFrame()

    #logger.debug(f"Condition-Check: {input_data_snp.dropna(how='all').empty}, {input_data_indel.dropna(how='all').empty}, {input_data_expanded_indel.dropna(how='all').empty}")
    #logger.debug(f"Condition: {(not input_data_snp.dropna(how='all').empty) or ((not input_data_indel.dropna(how='all').empty) and (not input_data_expanded_indel.dropna(how='all').empty))}")

    if (not input_data_snp.dropna(how='all').empty) or ((not input_data_indel.dropna(how='all').empty) and (not input_data_expanded_indel.dropna(how='all').empty)):
        if ((not input_data_indel.empty) and (not input_data_expanded_indel.empty)):
            logger.info("Combine InDel variants ...")
            input_data_combined_indel = combine_expanded_indels.parallelized_indel_combination(input_data_indel, input_data_expanded_indel, feature_list, num_cores)
        
        else:
            logger.warn("No InDel variants given move on to SNP processing!")
            input_data_combined_indel = pd.DataFrame()

        # predict pathogenicity score
        logger.info("Score variants ...")
        
        if not input_data_snp.dropna(how='all').empty:
            predicted_data_snp = predict.perform_pathogenicity_score_prediction(scoring_model_snp, input_data_snp, allele_frequency_list, feature_list, num_cores)
        
        else:
            logger.warn("No SNP variants, skip SNP prediction!")
            predicted_data_snp = pd.DataFrame()

        if not input_data_combined_indel.dropna(how='all').empty:
            predicted_data_indel = predict.perform_pathogenicity_score_prediction(scoring_model_indel, input_data_combined_indel, allele_frequency_list, feature_list, num_cores)
        
        else:
            logger.warn("No InDel variants, skip InDel prediction!")
            predicted_data_indel = pd.DataFrame()

        if (not predicted_data_snp.dropna(how='all').empty) and (not predicted_data_indel.dropna(how='all').empty):
            predicted_data = pd.concat([predicted_data_snp, predicted_data_indel])
            predicted_data.sort_values(["CHROM", "POS"], ascending=[True, True], inplace=True)
            predicted_data.reset_index(inplace=True, drop=True)
            predicted_data = predicted_data[predicted_data_snp.columns]

        elif (predicted_data_snp.dropna(how='all').empty) and (not predicted_data_indel.dropna(how='all').empty):
            predicted_data = predicted_data_indel

        elif (predicted_data_indel.dropna(how='all').empty) and (not predicted_data_snp.dropna(how='all').empty):
            predicted_data = predicted_data_snp

        else:
            logger.critical("Something went terribly wrong!")

        # prioritize and filter variants
        logger.info("Prioritize variants and finalize score ...")
        prioritized_data = prio.prioritize_variants(predicted_data, internal_parameter_dict, ref_path, num_cores, assembly_build, feature_list, skip_db_check, family_file, family_type, hpo_file, gene_exclusion_file)

        ## TODO: create additional output files according to the inheritance information (only filtered data)
        if only_top_results:
            prioritized_data[prioritized_data["FILTER_PASSED"] == 1].head(n=25).to_csv(str(output_filename + "_aidiva_result_filtered.tsv"), sep="\t", index=False)
            logger.info("Only 25 best variants are reported as result!")

        else:
            if save_as_vcf:
                write_result.write_result_vcf(prioritized_data, str(output_filename + ".vcf"), assembly_build, bool(family_type == "SINGLE"))
                write_result.write_result_vcf(prioritized_data[prioritized_data["FILTER_PASSED"] == 1], str(output_filename + "_aidiva_result_filtered.vcf"), assembly_build, bool(family_type == "SINGLE"))

            prioritized_data = prioritized_data.rename(columns={"CHROM": "#CHROM"})
            prioritized_data.to_csv(str(output_filename + "_aidiva_result.tsv"), sep="\t", index=False)
            prioritized_data[prioritized_data["FILTER_PASSED"] == 1].to_csv(str(output_filename + "_aidiva_result_filtered.tsv"), sep="\t", index=False)

            if family_type == "TRIO": #TODO add additional condition for type FAMILY (this type is currently not supported)
                prioritized_data[(prioritized_data["FILTER_PASSED"] == 1) & (prioritized_data["COMPOUND"] == 1)].to_csv(str(output_filename + "_aidiva_result_filtered_compound.tsv"), sep="\t", index=False)
                prioritized_data[(prioritized_data["FILTER_PASSED"] == 1) & (prioritized_data["DOMINANT_DENOVO"] == 1)].to_csv(str(output_filename + "_aidiva_result_filtered_dominant_denovo.tsv"), sep="\t", index=False)
                prioritized_data[(prioritized_data["FILTER_PASSED"] == 1) & (prioritized_data["DOMINANT"] == 1)].to_csv(str(output_filename + "_aidiva_result_filtered_dominant.tsv"), sep="\t", index=False)
                prioritized_data[(prioritized_data["FILTER_PASSED"] == 1) & (prioritized_data["XLINKED"] == 1)].to_csv(str(output_filename + "_aidiva_result_filtered_xlinked.tsv"), sep="\t", index=False)
                prioritized_data[(prioritized_data["FILTER_PASSED"] == 1) & (prioritized_data["RECESSIVE"] == 1)].to_csv(str(output_filename + "_aidiva_result_filtered_recessive.tsv"), sep="\t", index=False)

        logger.info("Pipeline successfully finsished!")

    else:
        if save_as_vcf:
            write_result.write_result_vcf(None, str(output_filename + "_aidiva_result.vcf"), assembly_build, bool(family_type == "SINGLE"))
            write_result.write_result_vcf(None, str(output_filename + "_aidiva_result_filtered.vcf"), assembly_build, bool(family_type == "SINGLE"))

        logger.warn("The given input files were empty!")
