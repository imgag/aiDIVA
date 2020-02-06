import yaml
import argparse
import tempfile
import os
import pandas as pd
import util.add_segmentDuplication as add_segDup
import util.add_simpleRepeats as add_repeats
import util.add_abb_score as add_abb
import util.add_score_from_bigwig as add_score
import variant_scoring.predict as predict
import util.createTrainingVCF as vcf
import util.convert_vcf_to_csv as convert_vcf
import util.convert_indels_to_snps_and_create_vcf as expand_indels_and_create_vcf
import util.combine_expanded_indels_and_create_csv as combine_expanded_indels
import variant_prioritization.familySNP_gene_score as prio
import variant_annotation.annotate_with_vep as annotate


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
    annotate.call_vep_and_annotate_vcf(str(working_directory + input_filename + "_snps.vcf"), str(working_directory + input_filename + "_snps_annotated.vcf"), vep_annotation_dict, additional_bigwig_files, False)
    annotate.call_vep_and_annotate_vcf(str(working_directory + input_filename + "_indel.vcf"), str(working_directory + input_filename + "_indel_annotated.vcf"), vep_annotation_dict, additional_bigwig_files, True)
    annotate.call_vep_and_annotate_vcf(str(working_directory + input_filename + "_indel_expanded.vcf"), str(working_directory + input_filename + "_indel_expanded_annotated.vcf"), vep_annotation_dict, additional_bigwig_files, False)
    
    # convert annotated vcfs back to pandas dataframes
    input_data_snps_annotated = convert_vcf.convert_vcf_to_pandas_dataframe(str(working_directory + input_filename + "_snps_annotated.vcf"), False)
    input_data_indel_annotated = convert_vcf.convert_vcf_to_pandas_dataframe(str(working_directory + input_filename + "_indel_annotated.vcf"), True)
    input_data_indel_expanded_annotated = convert_vcf.convert_vcf_to_pandas_dataframe(str(working_directory + input_filename + "_indel_expanded_annotated.vcf"), True)
    
    input_data_snps_additional_annotation = add_score.group_and_process_data(additional_bigwig_files["phyloP46_mammal"], input_data_snps_annotated, "phyloP46_mammal")
    input_data_indel_expanded_additional_annotation = add_score.group_and_process_data(additional_bigwig_files["phyloP46_mammal"], input_data_indel_expanded_annotated, "phyloP46_mammal")
        
    input_data_snps_additional_annotation = add_score.group_and_process_data(additional_bigwig_files["phyloP46_primate"], input_data_snps_additional_annotation, "phyloP46_primate")
    input_data_indel_expanded_additional_annotation = add_score.group_and_process_data(additional_bigwig_files["phyloP46_primate"], input_data_indel_expanded_additional_annotation, "phyloP46_primate")
    
    input_data_snps_additional_annotation = add_score.group_and_process_data(additional_bigwig_files["phastCons46_mammal"], input_data_snps_additional_annotation, "phastCons46_mammal")
    input_data_indel_expanded_additional_annotation = add_score.group_and_process_data(additional_bigwig_files["phastCons46_mammal"], input_data_indel_expanded_additional_annotation, "phastCons46_mammal")
    
    input_data_snps_additional_annotation = add_score.group_and_process_data(additional_bigwig_files["phastCons46_primate"], input_data_snps_additional_annotation, "phastCons46_primate")
    input_data_indel_expanded_additional_annotation = add_score.group_and_process_data(additional_bigwig_files["phastCons46_primate"], input_data_indel_expanded_additional_annotation, "phastCons46_primate")
    
    print(input_data_snps_additional_annotation)
    print(input_data_indel_expanded_additional_annotation)
    
    # add abb score to annotation
    input_data_snps_additional_annotation = add_abb.group_and_process_data(abb_score_file, input_data_snps_additional_annotation)
    input_data_indel_expanded_additional_annotation = add_abb.group_and_process_data(abb_score_file, input_data_indel_expanded_additional_annotation)
    
    # add simple repeats to annotation
    input_data_snps_additional_annotation = add_repeats.group_and_process_data(simple_repeat_file, input_data_snps_additional_annotation)
    input_data_indel_expanded_additional_annotation = add_repeats.group_and_process_data(simple_repeat_file, input_data_indel_expanded_additional_annotation)
    
    # add segment duplication to annotation
    input_data_snps_additional_annotation = add_segDup.group_and_process_data(segment_duplication_file, input_data_snps_additional_annotation)
    input_data_indel_expanded_additional_annotation = add_segDup.group_and_process_data(segment_duplication_file, input_data_indel_expanded_additional_annotation)
    
    input_data_snps_additional_annotation.to_csv(str(working_directory + input_filename + "_snps_annotated.csv"), sep="\t", index=False)
    input_data_indel_annotated.to_csv(str(working_directory + input_filename + "_indel_annotated.csv"), sep="\t", index=False)
    input_data_indel_expanded_additional_annotation.to_csv(str(working_directory + input_filename + "_indel_expanded_annotated.csv"), sep="\t", index=False)
    
    # add annotation from addtional bigwig files
    # TODO wrap with a if condition to only use if no VEP annotation (is now part of the VEP annotation)
    #for key, value in additional_bigwig_files.items():
    #    input_data_snps_additional_annotation = add_score.group_and_process_data(value, input_data_snps_additional_annotation, key)
    #    input_data_indel_additional_annotation = add_score.group_and_process_data(value, input_data_indel_additional_annotation, key)
    
    # Workaround to make it work directly passing the tables to the combine method does not work (Why??? - maybe an error in pandas)
    input_data_indel_annotated = pd.read_csv(str(working_directory + input_filename + "_indel_annotated.csv"), sep="\t", low_memory=False)
    input_data_indel_expanded_additional_annotation = pd.read_csv(str(working_directory + input_filename + "_indel_expanded_annotated.csv"), sep="\t", low_memory=False)
    
    print(input_data_indel_annotated)
    print(input_data_indel_expanded_additional_annotation)
    
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

    #input_data_predicted.to_csv("tmp.predicted.tsv", index=False, sep="\t")
    
    #tmp = tempfile.NamedTemporaryFile(mode="w")
    #predicted_data_combined.to_csv(tmp.name, index=False, sep="\t")
    predicted_data_complete.to_csv(str(working_directory + input_filename + "_complete_annotated_predicted.csv"), index=False, sep="\t")
    
    #prio.main_program("tmp.predicted.tsv", output_tsv, output_filtered_tsv, family_file, inheritance, family_type, hpo_file, gene_exclusion_file)
    
    prio.main_program(str(working_directory + input_filename + "_complete_annotated_predicted.csv"), output_tsv, output_filtered_tsv, family_file, inheritance, family_type, hpo_file, gene_exclusion_file)
    #tmp.close()

