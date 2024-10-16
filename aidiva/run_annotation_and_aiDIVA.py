import argparse
import helper_modules.combine_expanded_indels_and_create_csv as combine_expanded_indels
import helper_modules.create_result_vcf as write_result
import helper_modules.convert_indels_to_snps_and_create_vcf as expand_indels_and_create_vcf
import helper_modules.convert_vcf_to_csv as convert_vcf
import helper_modules.filter_vcf as filt_vcf
import helper_modules.split_vcf_in_indel_and_snp_set as split_vcf
import logging
import os
import pandas as pd
import tempfile
import time
import variant_scoring.score_variants as predict
import variant_prioritization.prioritize_variants as prio
import variant_annotation.annotate_with_vep as annotate
import yaml


if __name__=="__main__":
    parser = argparse.ArgumentParser(description = "AIdiva -- Augmented Intelligence Disease Variant Analysis")
    parser.add_argument("--vcf", type=str, dest="vcf", metavar="input.vcf(.gz)", required=True, help="VCF file with the variants to analyze [required]")
    parser.add_argument("--config", type=str, dest="config", metavar="config.yaml", required=True, help="Config file specifying the parameters for aiDIVA [required]")
    parser.add_argument("--out_prefix", type=str, dest="out_prefix", metavar="/output_path/aidiva_result", required=True, help="Prefix that is used to save the results [required]")
    parser.add_argument("--workdir", type=str, dest="workdir", metavar="/tmp/aidiva_workdir/", required=False, help="Path to the working directory, here all intermediate files are saved (if not specified a temporary folder will be created and used)")
    parser.add_argument("--hpo_list", type=str, dest="hpo_list", metavar="hpo.txt", required=False, help="TXT file containing the HPO terms reported for the current patient")
    parser.add_argument("--gene_exclusion", type=str, dest="gene_exclusion", metavar="gene_exclusion.txt", required=False, help="TXT file containing the genes to exclude in the analysis")
    parser.add_argument("--family_file", type=str, dest="family_file", metavar="family.txt", required=False, help="TXT file showing the family relation of the current patient")
    parser.add_argument("--family_type", type=str, dest="family_type", metavar="SINGLE", required=False, help="String indicating the present family type [SINGLE, TRIO]")
    parser.add_argument("--skip_db_check", dest="skip_db_check", action="store_true", required=False, help="Flag to skip DB lookup of variants")
    parser.add_argument("--only_top_results", dest="only_top_results", action="store_true", required=False, help="Report only the top ranking variants as result. The desired rank can be given as parameter with '--top_rank' (default: 25)")
    parser.add_argument("--top_rank", type=str, dest="top_rank", metavar="25", required=False, help="Rank parameter for '--only_top_results' (default: 25)")
    parser.add_argument("--inhouse_sample", dest="inhouse_sample", action="store_true", required=False, help="Flag to indicate that we are annotating an inhouse sample (skips leftNormalize since it is already performed)")
    parser.add_argument("--threads", type=int, dest="threads", metavar="1", required=False, help="Number of threads to use. (default: 1)")
    parser.add_argument("--log_file", type=str, dest="log_file", metavar="/output_path/logs/aidiva_log.txt", required=False, help="Path plus name of the log file to be saved, if not specified the log file is saved in the working directory")
    parser.add_argument("--log_level", type=str, dest="log_level", metavar="INFO", required=False, help="Define logging level, if unsure just leave the default [DEBUG, INFO] (default: INFO)")
    parser.add_argument("--save_as_vcf", dest="save_as_vcf", action="store_true", required=False, help="Flag to create additional result files in VCF format.")
    args = parser.parse_args()

    # parse input files
    if os.path.isfile(args.vcf):
        input_vcf = args.vcf

    else:
        raise SystemExit("The given input file could not be opened!")

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

    if args.family_file is not None:
        family_file = args.family_file

        if args.family_type is not None:
            family_type = args.family_type

        else:
            family_type = "SINGLE"

    else:
        family_file = None
        family_type = "SINGLE"

    save_as_vcf = args.save_as_vcf
    skip_db_check = args.skip_db_check
    only_top_results = args.only_top_results
    inhouse_sample = args.inhouse_sample

    # obtain number of threads to use during computation
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

        #elif args.log_level == "WARNING":
        #    log_level = logging.WARNING
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

    logger.info("Running annotation and AIdiva")
    logger.info("Start program")
    logger.info(f"Working directory: {working_directory}")

    # load ML model
    scoring_model = configuration["Analysis-Input"]["scoring-model"]

    annotation_dict = configuration["Annotation-Resources"]
    prioritization_information_dict = configuration["Analysis-Input"]["prioritization-information"]
    internal_parameter_dict = configuration["Internal-Parameters"]

    allele_frequency_list = configuration["Model-Features"]["allele-frequency-list"]
    feature_list = configuration["Model-Features"]["feature-list"]
    assembly_build = configuration["Assembly-Build"]

    ref_path = configuration["Analysis-Input"]["ref-path"]

    transcript_path = configuration["Canonical-Transcripts"]

    # convert splitted input data to vcf and annotate
    input_file = os.path.splitext(input_vcf)[0]
    input_filename = os.path.basename(input_file)
    input_filename = input_filename.split(".")[0]

    logger.info("Starting VCF preparation...")
    # left normalizing, sorting, and filtering step to remove unsupported variants
    annotate.left_normalize_and_sort_vcf(input_vcf, str(working_directory + "/" + input_filename + "_sorted.vcf"), annotation_dict, ref_path, inhouse_sample)
    annotate.annotate_consequence_information(str(working_directory + "/" + input_filename + "_sorted.vcf"), str(working_directory + "/" + input_filename + "_consequence.vcf"), annotation_dict, assembly_build, num_cores)
    filt_vcf.filter_coding_variants(str(working_directory + "/" + input_filename + "_consequence.vcf"), str(working_directory + "/" + input_filename + "_filtered.vcf"), "CONS")

    # convert input vcf to pandas dataframe
    split_vcf.split_vcf_file_in_indel_and_snps_set(str(working_directory + "/" + input_filename + "_filtered.vcf"), str(working_directory + "/" + input_filename + "_snp.vcf"), str(working_directory + "/" + input_filename + "_indel.vcf"))
    expand_indels_and_create_vcf.convert_indel_vcf_to_expanded_indel_vcf(str(working_directory + "/" + input_filename + "_indel.vcf"), str(working_directory + "/" + input_filename + "_indel_expanded.vcf"), ref_path)

    # Annotation with VEP
    logger.info("Starting VEP annotation...")
    annotate.call_vep_and_annotate_vcf(str(working_directory + "/" + input_filename + "_snp.vcf"), str(working_directory + "/" + input_filename + "_snp_vep.vcf"), annotation_dict, assembly_build, False, False, num_cores)
    annotate.call_vep_and_annotate_vcf(str(working_directory + "/" + input_filename + "_indel.vcf"), str(working_directory + "/" + input_filename + "_indel_vep.vcf"), annotation_dict, assembly_build, True, False, num_cores)
    annotate.call_vep_and_annotate_vcf(str(working_directory + "/" + input_filename + "_indel_expanded.vcf"), str(working_directory + "/" + input_filename + "_indel_expanded_vep.vcf"), annotation_dict, assembly_build, False, True, num_cores)

    # Additional annotation with AnnotateFromVCF (a ngs-bits tool)
    logger.info("Starting AnnotateFromVCF annotation...")
    annotate.annotate_from_vcf(str(working_directory + "/" + input_filename + "_snp_vep.vcf"), str(working_directory + "/" + input_filename + "_snp_vep_annotated.vcf"), annotation_dict, False, False, num_cores)
    annotate.annotate_from_vcf(str(working_directory + "/" + input_filename + "_indel_vep.vcf"), str(working_directory + "/" + input_filename + "_indel_vep_annotated.vcf"), annotation_dict, False, True, num_cores)
    annotate.annotate_from_vcf(str(working_directory + "/" + input_filename + "_indel_expanded_vep.vcf"), str(working_directory + "/" + input_filename + "_indel_expanded_vep_annotated.vcf"), annotation_dict, True, False, num_cores)

    # Additional annotation with AnnotateFromBed (a ngs-bits tool)
    annotate.annotate_from_bed(str(working_directory + "/" + input_filename + "_snp_vep_annotated.vcf"), str(working_directory + "/" + input_filename + "_snp_vep_annotated_bed.vcf"), annotation_dict, num_cores)
    annotate.annotate_from_bed(str(working_directory + "/" + input_filename + "_indel_vep_annotated.vcf"), str(working_directory + "/" + input_filename + "_indel_vep_annotated_bed.vcf"), annotation_dict, num_cores)

    # Additional annotation with AnnotateFromBigWig (a ngs-bits tool)
    annotate.annotate_from_bigwig(str(working_directory + "/" + input_filename + "_snp_vep_annotated_bed.vcf"), str(working_directory + "/" + input_filename + "_snp_vep_annotated_bed_bw.vcf"), annotation_dict, num_cores)
    annotate.annotate_from_bigwig(str(working_directory + "/" + input_filename + "_indel_expanded_vep_annotated.vcf"), str(working_directory + "/" + input_filename + "_indel_expanded_vep_annotated_bw.vcf"), annotation_dict, num_cores)

    # Filter low confidence regions with VariantFilterRegions (a ngs-bits tool)
    annotate.filter_regions(str(working_directory + "/" + input_filename + "_snp_vep_annotated_bed_bw.vcf"), str(working_directory + "/" + input_filename + "_snp_vep_annotated_bed_bw_filtered.vcf"), annotation_dict)
    annotate.filter_regions(str(working_directory + "/" + input_filename + "_indel_vep_annotated_bed.vcf"), str(working_directory + "/" + input_filename + "_indel_vep_annotated_bed_filtered.vcf"), annotation_dict)

    # convert annotated vcfs back to pandas dataframes
    input_data_snp_annotated = convert_vcf.convert_vcf_to_pandas_dataframe(str(working_directory + "/" + input_filename + "_snp_vep_annotated_bed_bw_filtered.vcf"), False, False, transcript_path, num_cores)
    input_data_indel_annotated = convert_vcf.convert_vcf_to_pandas_dataframe(str(working_directory + "/" + input_filename + "_indel_vep_annotated_bed_filtered.vcf"), True, False, transcript_path, num_cores)
    input_data_indel_expanded_annotated = convert_vcf.convert_vcf_to_pandas_dataframe(str(working_directory + "/" + input_filename + "_indel_expanded_vep_annotated_bw.vcf"), True, True, transcript_path, num_cores)

    if (not input_data_snp_annotated.dropna(how='all').empty) or ((not input_data_indel_annotated.dropna(how='all').empty) and (not input_data_indel_expanded_annotated.dropna(how='all').empty)):
        if ((not input_data_indel_annotated.empty) and (not input_data_indel_expanded_annotated.empty)):
            # combine the two indel sets
            logger.info("Combine InDel variants ...")
            input_data_indel_combined_annotated = combine_expanded_indels.parallelized_indel_combination(input_data_indel_annotated, input_data_indel_expanded_annotated, feature_list, num_cores)

        else:
            logger.warning("No InDel variants given move on to SNV processing!")
            input_data_indel_combined_annotated = pd.DataFrame()

        # predict pathogenicity score
        logger.info("Score variants...")

        if not input_data_snp_annotated.dropna(how='all').empty:
            predicted_data_snp = predict.perform_pathogenicity_score_prediction(scoring_model, input_data_snp_annotated, allele_frequency_list, feature_list, num_cores)

        else:
            logger.warning("No SNV variants, skip SNV prediction!")
            predicted_data_snp = pd.DataFrame()

        if not input_data_indel_combined_annotated.dropna(how='all').empty:
            predicted_data_indel = predict.perform_pathogenicity_score_prediction(scoring_model, input_data_indel_combined_annotated, allele_frequency_list, feature_list, num_cores)

        else:
            logger.warning("No InDel variants, skip InDel prediction!")
            predicted_data_indel = pd.DataFrame()

        if (not predicted_data_snp.dropna(how='all').empty) and (not predicted_data_indel.dropna(how='all').empty):
            predicted_data = pd.concat([predicted_data_snp, predicted_data_indel])
            predicted_data.sort_values(["#CHROM", "POS"], ascending=[True, True], inplace=True)
            predicted_data.reset_index(inplace=True, drop=True)
            predicted_data = predicted_data[predicted_data_snp.columns]

        elif (predicted_data_snp.dropna(how='all').empty) and (not predicted_data_indel.dropna(how='all').empty):
            predicted_data = predicted_data_indel

        elif (predicted_data_indel.dropna(how='all').empty) and (not predicted_data_snp.dropna(how='all').empty):
            predicted_data = predicted_data_snp

        else:
            logger.error("Something went terribly wrong!")

        # prioritize and filter variants
        logger.info("Filter variants and finalize score...")
        prioritized_data = prio.prioritize_variants(predicted_data, internal_parameter_dict, ref_path, num_cores, assembly_build, feature_list, skip_db_check, family_file, family_type, hpo_file, gene_exclusion_file)
        prioritized_data["AIDIVA_RANK"] = prioritized_data["FINAL_AIDIVA_SCORE"].rank(method='min', ascending=False)

        ## TODO: create additional output files according to the inheritance information (only filtered data)
        if only_top_results:
            prioritized_data[prioritized_data["FILTER_PASSED"] == 1 & prioritized_data["AIDIVA_RANK"] <= top_rank].to_csv(str(output_filename + "_aidiva_result_filtered.tsv"), sep="\t", index=False)
            logger.info(f"Only top ranking variants are reported as result (rank: {top_rank})!")

        else:
            if save_as_vcf:
                write_result.write_result_vcf(prioritized_data, str(output_filename + "_aidiva_result.vcf"), assembly_build, bool(family_type == "SINGLE"))
                write_result.write_result_vcf(prioritized_data[prioritized_data["FILTER_PASSED"] == 1], str(output_filename + "_aidiva_result_filtered.vcf"), assembly_build, bool(family_type == "SINGLE"))

            #prioritized_data = prioritized_data.rename(columns={"CHROM": "#CHROM"}) # not needed anymore
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

        logger.warning("The given input files were empty!")
