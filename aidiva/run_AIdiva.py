import argparse
import os
import ntpath
import pandas as pd
import tempfile
import helper_modules.combine_expanded_indels_and_create_csv as combine_expanded_indels
import helper_modules.create_result_vcf as write_result
import helper_modules.convert_indels_to_snps_and_create_vcf as expand_indels_and_create_vcf
import helper_modules.convert_vcf_to_csv as convert_vcf
import variant_scoring.score_variants as predict
import variant_prioritization.prioritize_variants as prio
import variant_annotation.add_abb_score as add_abb
import variant_annotation.add_score_from_bigwig as add_score
import variant_annotation.add_segmentDuplication as add_segDup
import variant_annotation.add_simpleRepeats as add_repeats
import variant_annotation.annotate_with_vep as annotate
import yaml


if __name__=="__main__":
    parser = argparse.ArgumentParser(description = "AIdiva -- Augmented Intelligence Disease Variant Analysis")
    parser.add_argument("--config", type=str, dest="config", metavar="config.yaml", required=True, help="Config file specifying the parameters for AIdiva [required]")
    args = parser.parse_args()

    # parse configuration file
    config_file = open(args.config, "r")
    configuration = yaml.full_load(config_file)
    config_file.close()

    working_directory = configuration["Analysis-Input"]["work-dir"]

    if not working_directory.endswith("/"):
        working_directory = working_directory + "/"

    ref_path = configuration["Analysis-Input"]["ref-path"]

    # parse input files
    snp_vcf = configuration["Analysis-Input"]["vcf-snp"]
    indel_vcf = configuration["Analysis-Input"]["vcf-indel"]
    expanded_indel_vcf = configuration["Analysis-Input"]["vcf-expanded-indel"]

    # get machine learning models
    scoring_model_snp = configuration["Analysis-Input"]["scoring-model-snp"]
    scoring_model_indel = configuration["Analysis-Input"]["scoring-model-indel"]
    coding_region_file = configuration["Analysis-Input"]["coding-region"]

    # parse output files
    output_filename = configuration["Analysis-Output"]["out-filename"]

    # parse disease and inheritance information
    hpo_file = configuration["Analysis-Input"]["prioritization-information"]["hpo-list"]
    gene_exclusion_file = configuration["Analysis-Input"]["prioritization-information"]["gene-exclusion"]
    family_type = configuration["Analysis-Input"]["prioritization-information"]["family-type"]
    family_file = configuration["Analysis-Input"]["prioritization-information"]["family-file"]

    hpo_resources_folder = configuration["Internal-Parameters"]["hpo-resources"]

    ## TODO: Choose whether to delete the allele frequency list
    allele_frequency_list = configuration["Model-Features"]["allele-frequency-list"]
    feature_list = configuration["Model-Features"]["feature-list"]

    # convert splitted input data to vcf and annotate
    #snp_vcf_filename = ntpath.basename(snp_vcf)
    #indel_vcf_file = os.path.splitext(indel_vcf)[0]
    #indel_vcf_filename = os.path.basename(indel_vcf_file)
    #expanded_indel_vcf_file = os.path.splitext(expanded_indel_vcf)[0]
    #expanded_indel_vcf_filename = os.path.basename(expanded_indel_vcf_file)


    input_data_snp = convert_vcf.convert_vcf_to_pandas_dataframe(snp_vcf, False)
    input_data_indel = convert_vcf.convert_vcf_to_pandas_dataframe(indel_vcf, True)
    input_data_expanded_indel = convert_vcf.convert_vcf_to_pandas_dataframe(expanded_indel_vcf, True)

    input_data_combined_indel = combine_expanded_indels.combine_vcf_dataframes(input_data_indel, input_data_expanded_indel, feature_list)

    # predict pathogenicity score
    print("Score variants ...")
    coding_region = pd.read_csv(coding_region_file, sep="\t", names=["CHROM", "START", "END"], low_memory=False)
    predicted_data = predict.perform_pathogenicity_score_prediction(input_data_snp, input_data_combined_indel, scoring_model_snp, scoring_model_indel, allele_frequency_list, feature_list, coding_region)

    # prioritize and filter variants
    print("Filter variants and finalize score ...")
    prioritized_data = prio.prioritize_variants(predicted_data, family_file, family_type, hpo_resources_folder, hpo_file, gene_exclusion_file)

    write_result.write_result_vcf(prioritized_data, str(working_directory + output_filename + ".vcf"))
    prioritized_data.to_csv(str(working_directory + output_filename + ".csv"), sep="\t", index=False)
    print("Pipeline successfully finsished!")
