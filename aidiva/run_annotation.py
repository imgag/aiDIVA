import argparse
import helper_modules.convert_indels_to_snps_and_create_vcf as expand_indels_and_create_vcf
import helper_modules.filter_vcf as filt_vcf
import helper_modules.split_vcf_in_indel_and_snp_set as split_vcf
import logging
import os
import tempfile
import time
import variant_annotation.annotate_with_vep as annotate
import yaml


if __name__=="__main__":
    parser = argparse.ArgumentParser(description = "AIdiva -- Augmented Intelligence Disease Variant Analysis")
    parser.add_argument("--vcf", type=str, dest="vcf", metavar="input.vcf(.gz)", required=True, help="VCF file with the variants to annotate [required]")
    parser.add_argument("--config", type=str, dest="config", metavar="config.yaml", required=True, help="Config file specifying the parameters for AIdiva [required]")
    parser.add_argument("--out_folder", type=str, dest="out_folder", metavar="/output_path/aidiva_result", required=True, help="Prefix that is used to save the annotated files [required]")
    parser.add_argument("--threads", type=int, dest="threads", metavar="1", required=False, help="Number of threads to use. (default: 1)")
    parser.add_argument("--log_path", type=str, dest="log_path", metavar="/output_path/logs/", required=False, help="Path where the log file should be saved, if not specified the log file is saved in the working directory")
    parser.add_argument("--log_level", type=str, dest="log_level", metavar="INFO", required=False, help="Define logging level, if unsure just leave the default [DEBUG, INFO] (default: INFO)")
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
    if args.out_folder is not None:
        output_folder = args.out_folder
    else:
        raise SystemExit("The parameter 'out_prefix' was not specified!")

    workdir = tempfile.TemporaryDirectory(suffix="", prefix="aidiva_workdir_", dir=None)
    working_directory = workdir.name

    if not working_directory.endswith("/"):
        working_directory = working_directory + "/"

    # obtain number of threads to use during computation
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

    logger.info("Running annotation")
    logger.info("Start program")
    logger.info(f"Output directory: {output_folder}")

    annotation_dict = configuration["Annotation-Resources"]

    assembly_build = configuration["Assembly-Build"]

    ref_path = configuration["Analysis-Input"]["ref-path"]

    # convert splitted input data to vcf and annotate
    input_file = os.path.splitext(input_vcf)[0]
    input_filename = os.path.basename(input_file)
    input_filename = input_filename.split(".")[0]

    logger.info("Starting VCF preparation...")
    # sorting and filtering step to remove unsupported variants
    annotate.sort_vcf(input_vcf, str(working_directory + "/" + input_filename + "_sorted.vcf"), annotation_dict)
    annotate.annotate_consequence_information(str(working_directory + "/" + input_filename + "_sorted.vcf"), str(working_directory + "/" + input_filename + "_consequence.vcf"), annotation_dict, assembly_build, num_cores)
    filt_vcf.filter_coding_variants(str(working_directory + "/" + input_filename + "_consequence.vcf"), str(output_folder + "/" + input_filename + "_filtered.vcf"), "CONS")

    # convert input vcf to pandas dataframe
    split_vcf.split_vcf_file_in_indel_and_snps_set(str(output_folder + "/" + input_filename + "_filtered.vcf"), str(working_directory + "/" + input_filename + "_snp.vcf"), str(working_directory + "/" + input_filename + "_indel.vcf"))
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
    annotate.annotate_from_bigwig(str(working_directory + "/" + input_filename + "_indel_expanded_vep_annotated.vcf"), str(output_folder + "/" + input_filename + "_indelExpanded_annotated.vcf"), annotation_dict, num_cores)

    # Filter low confidence regions with VariantFilterRegions (a ngs-bits tool)
    annotate.filter_regions(str(working_directory + "/" + input_filename + "_snp_vep_annotated_bed_bw.vcf"), str(output_folder + "/" + input_filename + "_snp_annotated.vcf"), annotation_dict)
    annotate.filter_regions(str(working_directory + "/" + input_filename + "_indel_vep_annotated_bed.vcf"), str(output_folder + "/" + input_filename + "_indel_annotated.vcf"), annotation_dict)

