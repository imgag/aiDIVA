import argparse
import os
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
    #parser.add_argument("--annotate", dest="annotate", action="store_true", required=False, help="Flag indicating that the annotation with VEP needs to be done")
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
    input_vcf = configuration["Analysis-Input"]["vcf"] # should already be annotated with VEP
    scoring_model_snp = coding_region_file = os.path.dirname(__file__) + "/../data/rf_model_snp_scikit0-19-1.pkl"
    scoring_model_indel = coding_region_file = os.path.dirname(__file__) + "/../data/rf_model_inframe_indel_scikit0-19-1.pkl"
    coding_region_file = os.path.dirname(__file__) + "/../data/GRCh37_coding_sequences.bed"
    #scoring_model_snp = configuration["Analysis-Input"]["scoring-model-snp"]
    #scoring_model_indel = configuration["Analysis-Input"]["scoring-model-indel"]
    #coding_region_file = configuration["Analysis-Input"]["coding-region"]

    # parse output files
    output_filename = configuration["Analysis-Output"]["out-filename"]

    # parse disease and inheritance information
    hpo_file = configuration["Analysis-Input"]["prioritization-information"]["hpo-list"]
    gene_exclusion_file = configuration["Analysis-Input"]["prioritization-information"]["gene-exclusion"]
    family_type = configuration["Analysis-Input"]["prioritization-information"]["family-type"]
    family_file = configuration["Analysis-Input"]["prioritization-information"]["family-file"]

    vep_annotation_dict = configuration["VEP-Annotation"]
    prioritization_information_dict = configuration["Analysis-Input"]["prioritization-information"]
    internal_parameter_dict = configuration["Internal-Parameters"]

    allele_frequency_list = configuration["Model-Features"]["allele-frequency-list"]
    feature_list = configuration["Model-Features"]["feature-list"]

    # parse additional annotation files
    abb_score_file = configuration["Additional-Annotation"]["abb-score"]
    simple_repeat_file = configuration["Additional-Annotation"]["simple-repeats"]
    segment_duplication_file = configuration["Additional-Annotation"]["segment-duplication"]
    additional_bigwig_files = configuration["Additional-Annotation"]["additional-bigwig-files"]

    # convert splitted input data to vcf and annotate
    input_file = os.path.splitext(input_vcf)[0]
    input_filename = os.path.basename(input_file)

    # convert input vcf to pandas dataframe
    convert_vcf.split_vcf_file_in_indel_and_snps_set(input_vcf, str(working_directory + input_filename + "_snp.vcf"), str(working_directory + input_filename + "_indel.vcf"))
    expand_indels_and_create_vcf.convert_csv_to_vcf(str(working_directory + input_filename + "_indel.vcf"), str(working_directory + input_filename + "_indel_expanded.vcf"), ref_path)

    #input_data_snps = input_data[(input_data["Ref"].apply(len) == 1) & (input_data["Alt"].apply(len) == 1)]
    #input_data_indel = input_data[(input_data["Ref"].apply(len) > 1) | (input_data["Alt"].apply(len) > 1)]

    # annotate with VEP
    #if vep_annotation_dict["perform-vep-annotation"]:
    ## TODO: change to match the two input files
    print("Starting VEP annotation ...")
    annotate.call_vep_and_annotate_vcf(str(working_directory + input_filename + "_snp.vcf"), str(working_directory + input_filename + "_snp_annotated.vcf"), vep_annotation_dict, False)
    annotate.call_vep_and_annotate_vcf(str(working_directory + input_filename + "_indel.vcf"), str(working_directory + input_filename + "_indel_annotated.vcf"), vep_annotation_dict, True)
    annotate.call_vep_and_annotate_vcf(str(working_directory + input_filename + "_indel_expanded.vcf"), str(working_directory + input_filename + "_indel_expanded_annotated.vcf"), vep_annotation_dict, False)
    print("Finished VEP annotation!")

    # convert annotated vcfs back to pandas dataframes
    input_data_snp_annotated = convert_vcf.convert_vcf_to_pandas_dataframe(str(working_directory + input_filename + "_snp_annotated.vcf"), False)
    input_data_indel_annotated = convert_vcf.convert_vcf_to_pandas_dataframe(str(working_directory + input_filename + "_indel_annotated.vcf"), True)
    input_data_indel_expanded_annotated = convert_vcf.convert_vcf_to_pandas_dataframe(str(working_directory + input_filename + "_indel_expanded_annotated.vcf"), True)
    #elif configuration["Additional-Annotation"]["perform-additional-annotation"]:
    #    input_data_snps_annotated = convert_vcf.convert_vcf_to_pandas_dataframe(str(working_directory + input_filename + "_snp.vcf"), False)
    #    input_data_indel_annotated = convert_vcf.convert_vcf_to_pandas_dataframe(str(working_directory + input_filename + "_indel.vcf"), True)
    #    input_data_indel_expanded_annotated = convert_vcf.convert_vcf_to_pandas_dataframe(str(working_directory + input_filename + "_indel_expanded.vcf"), True)

    # add conservation scores from bigwig files
    # the following is only needed if not all custom annotations can be made with VEP (eg. if some plugins are missing)
    if configuration["Additional-Annotation"]["perform-additional-annotation"]:
        print("Starting conservation score annotation ...")
        for key, value in additional_bigwig_files.items():
            input_data_snp_annotated = add_score.group_and_process_data(value, input_data_snp_annotated, key)
            input_data_indel_expanded_annotated = add_score.group_and_process_data(value, input_data_indel_expanded_annotated, key)
        print("Finished conservation score annotation!")

        # add abb score to annotation
        if not abb_score_file is None:
            print("Starting ABB score annotation ...")
            input_data_snp_annotated = add_abb.group_and_process_data(abb_score_file, input_data_snp_annotated)
            input_data_indel_expanded_annotated = add_abb.group_and_process_data(abb_score_file, input_data_indel_expanded_annotated)
            print("Finished ABB score annotation!")

        # add simple repeats to annotation
        if not simple_repeat_file is None:
            print("Starting simpleRepeat annotation ...")
            input_data_snp_annotated = add_repeats.group_and_process_data(simple_repeat_file, input_data_snp_annotated)
            input_data_indel_expanded_annotated = add_repeats.group_and_process_data(simple_repeat_file, input_data_indel_expanded_annotated)
            print("Finished simpleRepeat annotation!")

        # add segment duplication to annotation
        if not segment_duplication_file is None:
            print("Starting segmentDuplication annotation ...")
            input_data_snp_annotated = add_segDup.group_and_process_data(segment_duplication_file, input_data_snp_annotated)
            input_data_indel_expanded_annotated = add_segDup.group_and_process_data(segment_duplication_file, input_data_indel_expanded_annotated)
            print("Finished segmentDuplication annotation!")

    # save annotated files
    #input_data_snp_annotated.to_csv(str(working_directory + input_filename + "_snp_annotated.csv"), sep="\t", index=False)
    #input_data_indel_annotated.to_csv(str(working_directory + input_filename + "_indel_annotated.csv"), sep="\t", index=False)
    #input_data_indel_expanded_annotated.to_csv(str(working_directory + input_filename + "_indel_expanded_annotated.csv"), sep="\t", index=False)

    # Workaround to make it work directly passing the tables to the combine method does not work (Why??? - maybe an error in pandas)
    # TODO fix this small bug
    #input_data_indel_annotated = pd.read_csv(str(working_directory + input_filename + "_indel_annotated.csv"), sep="\t", low_memory=False)
    #input_data_indel_expanded_annotated = pd.read_csv(str(working_directory + input_filename + "_indel_expanded_annotated.csv"), sep="\t", low_memory=False)

    ## TODO: get rid of multiple values in the feature columns

    # combine the two indel sets
    input_data_indel_combined_annotated = combine_expanded_indels.combine_vcf_dataframes(input_data_indel_annotated, input_data_indel_expanded_annotated, feature_list)
    #input_data_indel_combined_annotated.to_csv(str(working_directory + input_filename + "_indel_combined_annotated.csv"), sep="\t", index=False)

    # TODO decide how to handle allele ambiguity (especially if there are exactly two reported)
    #input_data["Alt"] = input_data["Alt"].map(lambda x: x.split(",")[0])
    # for now just use the variants that have only one allele
    #input_data[input_data["Alt"].apply(lambda x: len(x.split(","))) == 1]

    #input_data_indel_combined_annotated = convert_vcf.convert_vcf_to_pandas_dataframe(str(working_directory + input_filename + "_indel.vcf"), True)

    # predict pathogenicity score
    print("Score variants ...")
    coding_region = pd.read_csv(coding_region_file, sep="\t", names=None, low_memory=False)
    predicted_data = predict.perform_pathogenicity_score_prediction(input_data_snp_annotated, input_data_indel_combined_annotated, scoring_model_snp, scoring_model_indel, allele_frequency_list, feature_list, coding_region)
    #predicted_data.to_csv(str(working_directory + "test_with_sample" + "_complete_annotated_predicted.csv"), index=False, sep="\t")

    # prioritize and filter variants
    print("Filter variants and finalize score ...")
    prioritized_data = prio.prioritize_variants(predicted_data, family_file, family_type, "/mnt/users/ahboced1/AIdiva_project/HPO_resources/", hpo_file, gene_exclusion_file)

    write_result.write_result_vcf(prioritized_data, str(working_directory + output_filename + ".vcf"), bool(family_type == "SINGLE"))
    prioritized_data.to_csv(str(working_directory + output_filename + ".csv"), sep="\t", index=False)
    print("Pipeline successfully finsished!")
