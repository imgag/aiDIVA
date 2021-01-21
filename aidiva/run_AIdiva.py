import argparse
import os
import pandas as pd
import tempfile
import time
import helper_modules.combine_expanded_indels_and_create_csv as combine_expanded_indels
import helper_modules.create_result_vcf as write_result
import helper_modules.convert_indels_to_snps_and_create_vcf as expand_indels_and_create_vcf
import helper_modules.convert_vcf_to_csv as convert_vcf
import variant_scoring.score_variants as predict
import variant_prioritization.prioritize_variants as prio
import yaml


if __name__=="__main__":
    parser = argparse.ArgumentParser(description = "AIdiva -- Augmented Intelligence Disease Variant Analysis")
    parser.add_argument("--snp_vcf", type=str, dest="snp_vcf", metavar="snp.vcf", required=True, help="VCF file with the annotated SNP variants [required]")
    parser.add_argument("--indel_vcf", type=str, dest="indel_vcf", metavar="indel.vcf", required=True, help="VCF file with the annotated (only basic annotation) InDel variants [required]")
    parser.add_argument("--expanded_indel_vcf", type=str, dest="expanded_indel_vcf", metavar="expanded_indel.vcf", required=True, help="VCF file with the annotated expanded InDel variants [required]")
    parser.add_argument("--out_prefix", type=str, dest="out_prefix", metavar="result", required=True, help="Prefix that is used to save the results [required]")
    parser.add_argument("--workdir", type=str, dest="workdir", metavar="workdir/", required=True, help="Path to the working directory (here all results are saved) [required]")
    parser.add_argument("--hpo_list", type=str, dest="hpo_list", metavar="hpo.txt", required=False, help="TXT file containing the HPO terms reported for the current patient")
    parser.add_argument("--family_file", type=str, dest="family_file", metavar="family.txt", required=False, help="TXT file showing the family relation of the current patient")
    parser.add_argument("--config", type=str, dest="config", metavar="config.yaml", required=True, help="Config file specifying the parameters for AIdiva [required]")
    parser.add_argument("--threads", type=int, dest="threads", metavar="1", nargs="?", const=1, required=False, help="Number of threads to use.")
    args = parser.parse_args()

    if "threads" in args:
        num_cores = args.threads
    else:
        num_cores = 1

    # parse configuration file
    config_file = open(args.config, "r")
    configuration = yaml.load(config_file, Loader=yaml.SafeLoader)
    config_file.close()

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
    if "hpo_list" in args:
        hpo_file = args.hpo_list
    else:
        hpo_file = None
    gene_exclusion_file = configuration["Analysis-Input"]["prioritization-information"]["gene-exclusion"]

    if "family_file" in args:
        family_file = args.family_file
        family_type = "SINGLE" ## TODO: get correct family type based on the family file
    else:
        family_file = None
        family_type = "SINGLE"

    ## TODO: Choose whether to delete the allele frequency list
    #allele_frequency_list = configuration["Model-Features"]["allele-frequency-list"]
    allele_frequency_list = []
    feature_list = configuration["Model-Features"]["feature-list"]

    # convert splitted input data to vcf and annotate
    input_data_snp = convert_vcf.convert_vcf_to_pandas_dataframe(snp_vcf, False, num_cores)
    input_data_indel = convert_vcf.convert_vcf_to_pandas_dataframe(indel_vcf, True, num_cores)
    input_data_expanded_indel = convert_vcf.convert_vcf_to_pandas_dataframe(expanded_indel_vcf, True, num_cores)

    ## TODO: handle the situation if one or more (but not all) of the input dataframes are empty
    ## TODO: make it work if only InDel or only SNP variants are given
    if (not input_data_snp.empty) & (not input_data_indel.empty) & (not input_data_expanded_indel.empty):
        print("Combine InDel variants ...")
        input_data_combined_indel = combine_expanded_indels.parallelized_indel_combination(input_data_indel, input_data_expanded_indel, feature_list, num_cores)

        # predict pathogenicity score
        print("Score variants ...")
        #predicted_data = predict.perform_pathogenicity_score_prediction(input_data_snp, input_data_combined_indel, scoring_model_snp, scoring_model_indel, allele_frequency_list, feature_list, num_cores)
        predicted_data_snp = predict.perform_pathogenicity_score_prediction(scoring_model_snp, input_data_snp, allele_frequency_list, feature_list, num_cores)
        predicted_data_indel = predict.perform_pathogenicity_score_prediction(scoring_model_indel, input_data_combined_indel, allele_frequency_list, feature_list, num_cores)

        predicted_data = pd.concat([predicted_data_snp, predicted_data_indel])
        predicted_data.sort_values(["CHROM", "POS"], ascending=[True, True], inplace=True)
        predicted_data.reset_index(inplace=True, drop=True)
        predicted_data = predicted_data[predicted_data_snp.columns]

        # prioritize and filter variants
        print("Filter variants and finalize score ...")
        prioritized_data = prio.prioritize_variants(predicted_data, hpo_resources_folder, family_file, family_type, hpo_file, gene_exclusion_file)

        write_result.write_result_vcf(prioritized_data, str(working_directory + output_filename + ".vcf"), bool(family_type == "SINGLE"))
        prioritized_data.to_csv(str(working_directory + output_filename + ".csv"), sep="\t", index=False)
        prioritized_data[prioritized_data["FILTER_PASSED"] == 1].to_csv(str(working_directory + output_filename + "_passed_filters.csv"), sep="\t", index=False)
        print("Pipeline successfully finsished!")
    else:
        print("ERROR: The given input files were empty!")
        #write_result.write_result_vcf(pd.concat([input_data_snp, input_data_indel, input_data_expanded_indel]), str(working_directory + output_filename + ".vcf"), bool(family_type == "SINGLE"))
