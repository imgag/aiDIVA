import argparse
import os
import pandas as pd
import tempfile
import util.combine_expanded_indels_and_create_csv as combine_expanded_indels
import util.convert_indels_to_snps_and_create_vcf as expand_indels_and_create_vcf
import util.convert_vcf_to_csv as convert_vcf
import variant_scoring.predict as predict
import variant_prioritization.familySNP_gene_score as prio
import variant_annotation.add_abb_score as add_abb
import variant_annotation.add_score_from_bigwig as add_score
import variant_annotation.add_segmentDuplication as add_segDup
import variant_annotation.add_simpleRepeats as add_repeats
import variant_annotation.annotate_with_vep as annotate
import yaml


if __name__=='__main__':
    parser = argparse.ArgumentParser(description = 'AIdiva -- Augmented Intelligence Disease Variant Analysis')
    parser.add_argument('--config', type=str, dest='config', metavar='config.yaml', required=True, help='Config file specifying the parameters for AIdiva [required]')
    parser.add_argument('--annotate', dest='annotate', action='store_true', required=False, help='Flag indicating that the annotation with VEP needs to be done')
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
    scoring_model_snps = configuration["Analysis-Input"]["scoring-model-snps"]
    scoring_model_indels = configuration["Analysis-Input"]["scoring-model-indels"]
    
    # parse output files
    output_tsv = configuration["Analysis-Output"]["tsv"]
    output_filtered_tsv = configuration["Analysis-Output"]["tsv-filt"]
    
    # parse disease and inheritance information
    hpo_file = configuration["Analysis-Input"]["inheritance-information"]["hpo-file"]
    gene_exclusion_file = configuration["Analysis-Input"]["inheritance-information"]["gene-exclusion"]
    inheritance = configuration["Analysis-Input"]["inheritance-information"]["inheritance"]
    family_type = configuration["Analysis-Input"]["inheritance-information"]["family-type"]
    family_file = configuration["Analysis-Input"]["inheritance-information"]["family-file"]
    
    vep_annotation_dict = configuration["VEP-Annotation"]
    
    # parse additional annotation files
    abb_score_file = configuration["Additional-Annotation"]["abb-score"]
    simple_repeat_file = configuration["Additional-Annotation"]["simple-repeats"]
    segment_duplication_file = configuration["Additional-Annotation"]["segment-duplication"]
    additional_bigwig_files = configuration["Additional-Annotation"]["additional-bigwig-files"]
    
    # convert splitted input data to vcf and annotate 
    input_file = os.path.splitext(input_vcf)[0]
    input_filename = os.path.basename(input_file)
    # TODO keep sample information during the conversion
    #expand_indels_and_create_vcf.convert_csv_to_vcf(str(working_directory + input_filename + "_snps.vcf"), input_data_snps)
    #expand_indels_and_create_vcf.convert_csv_to_vcf(str(working_directory + input_filename + "_indel.vcf"), input_data_indel)
    
    # convert input vcf to pandas dataframe
    # TODO keep sample information during the conversion
    #split_vcf_file_in_indel_and_snps_set(filepath, filepath_snps, filepath_indel)
    #input_data = convert_vcf.convert_vcf_to_pandas_dataframe(input_vcf)
    convert_vcf.split_vcf_file_in_indel_and_snps_set(input_vcf, str(working_directory + input_filename + "_snps.vcf"), str(working_directory + input_filename + "_indel.vcf"))
    expand_indels_and_create_vcf.convert_csv_to_vcf(str(working_directory + input_filename + "_indel_expanded.vcf"), str(working_directory + input_filename + "_indel.vcf"), ref_path)
    
    #input_data_snps = input_data[(input_data["Ref"].apply(len) == 1) & (input_data["Alt"].apply(len) == 1)]
    #input_data_indel = input_data[(input_data["Ref"].apply(len) > 1) | (input_data["Alt"].apply(len) > 1)]
    
    # annotate with VEP
    print("Starting VEP annotation ...")
    annotate.call_vep_and_annotate_vcf(str(working_directory + input_filename + "_snps.vcf"), str(working_directory + input_filename + "_snps_annotated.vcf"), vep_annotation_dict, additional_bigwig_files, False)
    annotate.call_vep_and_annotate_vcf(str(working_directory + input_filename + "_indel.vcf"), str(working_directory + input_filename + "_indel_annotated.vcf"), vep_annotation_dict, additional_bigwig_files, True)
    annotate.call_vep_and_annotate_vcf(str(working_directory + input_filename + "_indel_expanded.vcf"), str(working_directory + input_filename + "_indel_expanded_annotated.vcf"), vep_annotation_dict, additional_bigwig_files, False)
    print("Finished VEP annotation!")
    
    # convert annotated vcfs back to pandas dataframes
    input_data_snps_annotated = convert_vcf.convert_vcf_to_pandas_dataframe(str(working_directory + input_filename + "_snps_annotated.vcf"), False)
    input_data_indel_annotated = convert_vcf.convert_vcf_to_pandas_dataframe(str(working_directory + input_filename + "_indel_annotated.vcf"), True)
    input_data_indel_expanded_annotated = convert_vcf.convert_vcf_to_pandas_dataframe(str(working_directory + input_filename + "_indel_expanded_annotated.vcf"), True)
    
    # add conservation scores from bigwig files
    print("Starting conservation score annotation ...")
    for key, value in additional_bigwig_files.items():
        input_data_snps_additional_annotation = add_score.group_and_process_data(value, input_data_snps_annotated, key)
        input_data_indel_expanded_additional_annotation = add_score.group_and_process_data(value, input_data_indel_expanded_annotated, key)
    print("Finished conservation score annotation!")
    
    # add abb score to annotation
    print("Starting ABB score annotation ...")
    input_data_snps_additional_annotation = add_abb.group_and_process_data(abb_score_file, input_data_snps_additional_annotation)
    input_data_indel_expanded_additional_annotation = add_abb.group_and_process_data(abb_score_file, input_data_indel_expanded_additional_annotation)
    print("Finished ABB score annotation!")
    
    # add simple repeats to annotation
    print("Starting ABB score annotation ...")
    input_data_snps_additional_annotation = add_repeats.group_and_process_data(simple_repeat_file, input_data_snps_additional_annotation)
    input_data_indel_expanded_additional_annotation = add_repeats.group_and_process_data(simple_repeat_file, input_data_indel_expanded_additional_annotation)
    print("Finished ABB score annotation!")
    
    # add segment duplication to annotation
    print("Starting ABB score annotation ...")
    input_data_snps_additional_annotation = add_segDup.group_and_process_data(segment_duplication_file, input_data_snps_additional_annotation)
    input_data_indel_expanded_additional_annotation = add_segDup.group_and_process_data(segment_duplication_file, input_data_indel_expanded_additional_annotation)
    print("Finished ABB score annotation!")
    
    # save annotated files
    input_data_snps_additional_annotation.to_csv(str(working_directory + input_filename + "_snps_annotated.csv"), sep="\t", index=False)
    input_data_indel_annotated.to_csv(str(working_directory + input_filename + "_indel_annotated.csv"), sep="\t", index=False)
    input_data_indel_expanded_additional_annotation.to_csv(str(working_directory + input_filename + "_indel_expanded_annotated.csv"), sep="\t", index=False)
    
    # Workaround to make it work directly passing the tables to the combine method does not work (Why??? - maybe an error in pandas)
    # TODO fix this small bug
    input_data_indel_annotated = pd.read_csv(str(working_directory + input_filename + "_indel_annotated.csv"), sep="\t", low_memory=False)
    input_data_indel_expanded_additional_annotation = pd.read_csv(str(working_directory + input_filename + "_indel_expanded_annotated.csv"), sep="\t", low_memory=False)
    
    # combine the two indel sets
    input_data_indel_combined_additional_annotation = combine_expanded_indels.combine_vcf_dataframes(input_data_indel_annotated, input_data_indel_expanded_additional_annotation)
    input_data_indel_combined_additional_annotation.to_csv(str(working_directory + input_filename + "_indel_combined_annotated.csv"), sep="\t", index=False)
    
    # TODO decide how to handle allele ambiguity (especially if there are exactly two reported)
    #input_data["Alt"] = input_data["Alt"].map(lambda x: x.split(",")[0])
    # for now just use the variants that have only one allele
    #input_data[input_data["Alt"].apply(lambda x: len(x.split(","))) == 1]
    
    # predict pathogenicity score
    predicted_data_snps, predicted_data_indel = predict.perform_pathogenicity_score_prediction(input_data_snps_additional_annotation, input_data_indel_combined_additional_annotation, scoring_model_snps, scoring_model_indels)
    
    # combine snps and indel data
    predicted_data_complete = pd.concat([predicted_data_snps, predicted_data_indel], sort=False)
    predicted_data_complete.sort_values(['Chr', 'Pos'], ascending=[True, True])

    # score pathogenicity of the variants
    print("Score variants ...")
    predicted_data_complete.to_csv(str(working_directory + input_filename + "_complete_annotated_predicted.csv"), index=False, sep="\t")
    
    # prioritize and filter variants
    print("Filter variants and finalize score ...")
    prio.main_program(str(working_directory + input_filename + "_complete_annotated_predicted.csv"), output_tsv, output_filtered_tsv, family_file, inheritance, family_type, hpo_file, gene_exclusion_file)
    
    print("Pipeline successfully finsished!")
