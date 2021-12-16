import argparse
import os
import pandas as pd
import time
import helper_modules.combine_expanded_indels_and_create_csv as combine_expanded_indels
import helper_modules.create_result_vcf as write_result
import helper_modules.convert_vcf_to_csv as convert_vcf
import variant_scoring.score_variants as predict
import variant_prioritization.prioritize_variants as prio
import yaml
import logging


if __name__=="__main__":
    parser = argparse.ArgumentParser(description = "AIdiva -- Augmented Intelligence Disease Variant Analysis")
    parser.add_argument("--snp_vcf", type=str, dest="snp_vcf", metavar="snp.vcf", required=True, help="VCF file with the annotated SNP variants [required]")
    parser.add_argument("--indel_vcf", type=str, dest="indel_vcf", metavar="indel.vcf", required=True, help="VCF file with the annotated (only basic annotation) InDel variants [required]")
    parser.add_argument("--expanded_indel_vcf", type=str, dest="expanded_indel_vcf", metavar="expanded_indel.vcf", required=True, help="VCF file with the annotated expanded InDel variants [required]")
    parser.add_argument("--out_prefix", type=str, dest="out_prefix", metavar="result", required=True, help="Prefix that is used to save the results [required]")
    parser.add_argument("--workdir", type=str, dest="workdir", metavar="workdir/", required=True, help="Path to the working directory (here all results are saved) [required]")
    parser.add_argument("--hpo_list", type=str, dest="hpo_list", metavar="hpo.txt", required=False, help="TXT file containing the HPO terms reported for the current patient")
    parser.add_argument("--gene_exclusion", type=str, dest="gene_exclusion", metavar="gene_exclusion.txt", required=False, help="Tab separated file containing the genes to exclude in the analysis. Genes are assumed to be in the first column.")
    parser.add_argument("--family_file", type=str, dest="family_file", metavar="family.txt", required=False, help="TXT file showing the sample relations of the current data")
    parser.add_argument("--family_type", type=str, dest="family_type", metavar="SINGLE", required=False, help="In case of multisample data the kind of sample relation [SINGLE, TRIO, MULTI]")
    parser.add_argument("--config", type=str, dest="config", metavar="config.yaml", required=True, help="Config file specifying the parameters for AIdiva [required]")
    parser.add_argument("--reference", type=str, dest="reference", metavar="GRCh37.fa", required=True, help="Reference sequence to use as FASTA [required]")
    parser.add_argument("--threads", type=int, dest="threads", metavar="1", required=False, help="Number of threads to use (default: 1)")
    args = parser.parse_args()
    
    # set up logger
    timestamp = time.strftime("%Y%m%d-%H%M%S")
    logging.basicConfig(filename=str(args.workdir + "/" + args.out_prefix + "_aidiva_" + timestamp + ".log"),
                            filemode='a',
                            format='%(asctime)s -- %(name)s - %(levelname)s - %(message)s',
                            datefmt='%H:%M:%S',
                            level=logging.DEBUG)
    logger = logging.getLogger()
    logger.info("Running AIdiva on annotated data")
    logger.info("Start program")


    if args.threads is not None:
        num_cores = int(args.threads)
    else:
        num_cores = 1

    # parse configuration file
    with open(args.config, "r") as config_file:
        configuration = yaml.load(config_file, Loader=yaml.SafeLoader)
    
    working_directory = args.workdir

    if not working_directory.endswith("/"):
        working_directory = working_directory + "/"

    # parse input files
    snp_vcf = args.snp_vcf
    indel_vcf = args.indel_vcf
    expanded_indel_vcf = args.expanded_indel_vcf

    # get machine learning models
    scoring_model_snp = os.path.dirname(__file__) + "/../data/" + configuration["Analysis-Input"]["scoring-model-snp"]
    scoring_model_indel = os.path.dirname(__file__) + "/../data/" + configuration["Analysis-Input"]["scoring-model-indel"]
    hpo_resources_folder = os.path.dirname(__file__) + "/../data/hpo_resources/"

    # parse output files
    #output_filename = configuration["Analysis-Output"]["out-filename"]
    output_filename = args.out_prefix

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
    
    ref_path = args.reference

    allele_frequency_list = configuration["Model-Features"]["allele-frequency-list"]
    feature_list = configuration["Model-Features"]["feature-list"]

    # convert splitted input data to vcf and annotate
    input_data_snp = convert_vcf.convert_vcf_to_pandas_dataframe(snp_vcf, False, num_cores)
    input_data_indel = convert_vcf.convert_vcf_to_pandas_dataframe(indel_vcf, True, num_cores)
    input_data_expanded_indel = convert_vcf.convert_vcf_to_pandas_dataframe(expanded_indel_vcf, True, num_cores)

    # TODO: make it work if only InDel or only SNP variants are given
    if (not input_data_snp.empty) and (not input_data_indel.empty) and (not input_data_expanded_indel.empty):
        logger.info("Combine InDel variants ...")
        input_data_combined_indel = combine_expanded_indels.parallelized_indel_combination(input_data_indel, input_data_expanded_indel, feature_list, num_cores)

        # predict pathogenicity score
        logger.info("Score variants ...")
        predicted_data_snp = predict.perform_pathogenicity_score_prediction(scoring_model_snp, input_data_snp, allele_frequency_list, feature_list, num_cores)
        predicted_data_indel = predict.perform_pathogenicity_score_prediction(scoring_model_indel, input_data_combined_indel, allele_frequency_list, feature_list, num_cores)

        predicted_data = pd.concat([predicted_data_snp, predicted_data_indel])
        predicted_data.sort_values(["CHROM", "POS"], ascending=[True, True], inplace=True)
        predicted_data.reset_index(inplace=True, drop=True)
        predicted_data = predicted_data[predicted_data_snp.columns]

        # prioritize and filter variants
        logger.info("Prioritize variants and finalize score ...")
        prioritized_data = prio.prioritize_variants(predicted_data, hpo_resources_folder, ref_path, num_cores, family_file, family_type, hpo_file, gene_exclusion_file)

        write_result.write_result_vcf(prioritized_data, str(working_directory + output_filename + ".vcf"), bool(family_type == "SINGLE"))
        write_result.write_result_vcf(prioritized_data[prioritized_data["FILTER_PASSED"] == 1], str(working_directory + output_filename + "_filtered.vcf"), bool(family_type == "SINGLE"))
        prioritized_data = prioritized_data.rename(columns={"CHROM": "#CHROM"})
        prioritized_data.to_csv(str(working_directory + output_filename + ".tsv"), sep="\t", index=False)
        prioritized_data[prioritized_data["FILTER_PASSED"] == 1].to_csv(str(working_directory + output_filename + "_filtered.tsv"), sep="\t", index=False)
        logger.info("Pipeline successfully finsished!")
    else:
        write_result.write_result_vcf(pd.DataFrame(), str(working_directory + output_filename + ".vcf"), bool(family_type == "SINGLE"))
        write_result.write_result_vcf(pd.DataFrame(), str(working_directory + output_filename + "_filtered.vcf"), bool(family_type == "SINGLE"))
        logger.warn("The given input files were empty!")
