import argparse
import os
import pandas as pd
import tempfile
import helper_modules.combine_expanded_indels_and_create_csv as combine_expanded_indels
import helper_modules.create_result_vcf as write_result
import helper_modules.convert_indels_to_snps_and_create_vcf as expand_indels_and_create_vcf
import helper_modules.convert_vcf_to_csv as convert_vcf
import helper_modules.filter_vcf as filt_vcf
import helper_modules.split_vcf_in_indel_and_snp_set as split_vcf
import variant_scoring.score_variants as predict
import variant_prioritization.prioritize_variants as prio
#import variant_annotation.add_abb_score as add_abb
#import variant_annotation.add_score_from_bigwig as add_score
#import variant_annotation.add_segmentDuplication as add_segDup
#import variant_annotation.add_simpleRepeats as add_repeats
import variant_annotation.annotate_with_vep as annotate
import yaml


if __name__=="__main__":
    parser = argparse.ArgumentParser(description = "AIdiva -- Augmented Intelligence Disease Variant Analysis")
    parser.add_argument("--vcf", type=str, dest="vcf", metavar="input.vcf", required=True, help="VCF file with the variants to analyze [required]")
    #parser.add_argument("--out_prefix", type=str, dest="out_prefix", metavar="result", required=True, help="Prefix that is used to save the results [required]")
    parser.add_argument("--workdir", type=str, dest="workdir", metavar="workdir/", required=True, help="Path to the working directory (here all results are saved) [required]")
    parser.add_argument("--config", type=str, dest="config", metavar="config.yaml", required=True, help="Config file specifying the parameters for AIdiva [required]")
    parser.add_argument("--hpo_list", type=str, dest="hpo_list", metavar="hpo.txt", required=False, help="TXT file containing the HPO terms reported for the current patient")
    parser.add_argument("--gene_exclusion", type=str, dest="gene_exclusion", metavar="gene_exclusion.txt", required=False, help="TXT file containing the genes to exclude in the analysis")
    parser.add_argument("--family_file", type=str, dest="family_file", metavar="family.txt", required=False, help="TXT file showing the family relation of the current patient")
    parser.add_argument("--family_type", type=str, dest="family_type", metavar="SINGLE", required=False, help="String indicating the present family type [SINGLE, TRIO]")
    #parser.add_argument("--threads", type=int, dest="threads", metavar="1", required=False, help="Number of threads to use. (default: 1)")
    args = parser.parse_args()

    # parse configuration file
    config_file = open(args.config, "r")
    configuration = yaml.load(config_file, Loader=yaml.SafeLoader)
    config_file.close()

    working_directory = args.workdir

    if not working_directory.endswith("/"):
        working_directory = working_directory + "/"

    data_path = "/mnt/data/"
    # data_path = os.path.dirname(os.path.abspath(__file__)) + "/../data/"
    ref_path = configuration["Analysis-Input"]["ref-path"]

    # parse input files
    input_vcf = args.vcf
    scoring_model_snp = os.path.dirname(os.path.abspath(__file__)) + "/../data/" + configuration["Analysis-Input"]["scoring-model-snp"]
    scoring_model_indel = os.path.dirname(os.path.abspath(__file__)) + "/../data/" + configuration["Analysis-Input"]["scoring-model-indel"]
    
    # obtain number of threads to use during computation
    num_cores = int(configuration["VEP-Annotation"]["num-threads"])

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

    vep_annotation_dict = configuration["VEP-Annotation"]
    prioritization_information_dict = configuration["Analysis-Input"]["prioritization-information"]
    internal_parameter_dict = configuration["Internal-Parameters"]

    allele_frequency_list = configuration["Model-Features"]["allele-frequency-list"]
    feature_list = configuration["Model-Features"]["feature-list"]

    # convert splitted input data to vcf and annotate
    input_file = os.path.splitext(input_vcf)[0]
    input_filename = os.path.basename(input_file)
    input_filename = input_filename.split(".")[0]

    hpo_resources_folder = os.path.dirname(os.path.abspath(__file__)) + "/../data/hpo_resources/"

    print("Starting VCF preparation...")
    # sorting and filtering step to remove unsupported variants
    annotate.sort_vcf(input_vcf, str(working_directory + input_filename + "_sorted.vcf"), vep_annotation_dict)
    annotate.annotate_consequence_information(str(working_directory + input_filename + "_sorted.vcf"), str(working_directory + input_filename + "_consequence.vcf"), vep_annotation_dict, num_cores)
    filt_vcf.filter_coding_variants(str(working_directory + input_filename + "_consequence.vcf"), str(working_directory + input_filename + "_filtered.vcf"), "CONS")

    # convert input vcf to pandas dataframe
    split_vcf.split_vcf_file_in_indel_and_snps_set(str(working_directory + input_filename + "_filtered.vcf"), str(working_directory + input_filename + "_snp.vcf"), str(working_directory + input_filename + "_indel.vcf"))
    expand_indels_and_create_vcf.convert_csv_to_vcf(str(working_directory + input_filename + "_indel.vcf"), str(working_directory + input_filename + "_indel_expanded.vcf"), ref_path)

    # Annotation with VEP
    print("Starting VEP annotation...")
    annotate.call_vep_and_annotate_vcf(str(working_directory + input_filename + "_snp.vcf"), str(working_directory + input_filename + "_snp_vep.vcf"), vep_annotation_dict, False, False, num_cores)
    annotate.call_vep_and_annotate_vcf(str(working_directory + input_filename + "_indel.vcf"), str(working_directory + input_filename + "_indel_vep.vcf"), vep_annotation_dict, True, False, num_cores)
    annotate.call_vep_and_annotate_vcf(str(working_directory + input_filename + "_indel_expanded.vcf"), str(working_directory + input_filename + "_indel_expanded_vep.vcf"), vep_annotation_dict, False, True, num_cores)

    # Additional annotation with AnnotateFromVCF (a ngs-bits tool)
    # If VCF is used as output format VEP won't annotate from custom VCF files
    print("Starting AnnotateFromVCF annotation...")
    annotate.annotate_from_vcf(str(working_directory + input_filename + "_snp_vep.vcf"), str(working_directory + input_filename + "_snp_vep_annotated.vcf"), vep_annotation_dict, num_cores)
    annotate.annotate_from_vcf(str(working_directory + input_filename + "_indel_expanded_vep.vcf"), str(working_directory + input_filename + "_indel_expanded_vep_annotated.vcf"), vep_annotation_dict, num_cores)

    # convert annotated vcfs back to pandas dataframes
    input_data_snp_annotated = convert_vcf.convert_vcf_to_pandas_dataframe(str(working_directory + input_filename + "_snp_vep_annotated.vcf"), False, num_cores)
    input_data_indel_annotated = convert_vcf.convert_vcf_to_pandas_dataframe(str(working_directory + input_filename + "_indel_vep.vcf"), True, num_cores)
    input_data_indel_expanded_annotated = convert_vcf.convert_vcf_to_pandas_dataframe(str(working_directory + input_filename + "_indel_expanded_vep_annotated.vcf"), True, num_cores)

    # combine the two indel sets
    input_data_indel_combined_annotated = combine_expanded_indels.parallelized_indel_combination(input_data_indel_annotated, input_data_indel_expanded_annotated, feature_list, num_cores)

    # predict pathogenicity score
    print("Score variants...")
    #predicted_data = predict.perform_pathogenicity_score_prediction(input_data_snp_annotated, input_data_indel_combined_annotated, scoring_model_snp, scoring_model_indel, allele_frequency_list, feature_list, num_cores)
    predicted_data_snp = predict.perform_pathogenicity_score_prediction(scoring_model_snp, input_data_snp_annotated, allele_frequency_list, feature_list, num_cores)
    predicted_data_indel = predict.perform_pathogenicity_score_prediction(scoring_model_indel, input_data_indel_combined_annotated, allele_frequency_list, feature_list, num_cores)

    predicted_data = pd.concat([predicted_data_snp, predicted_data_indel])
    predicted_data.sort_values(["CHROM", "POS"], ascending=[True, True], inplace=True)
    predicted_data.reset_index(inplace=True, drop=True)
    predicted_data = predicted_data[predicted_data_snp.columns]

    # prioritize and filter variants
    print("Filter variants and finalize score...")
    prioritized_data = prio.prioritize_variants(predicted_data, hpo_resources_folder, num_cores, family_file, family_type, hpo_file, gene_exclusion_file)

    ## TODO: create additional output files according to the inheritance information (only filtered data)
    write_result.write_result_vcf(prioritized_data, str(working_directory + input_filename + "_aidiva_result.vcf"), bool(family_type == "SINGLE"))
    write_result.write_result_vcf(prioritized_data[prioritized_data["FILTER_PASSED"] == 1], str(working_directory + input_filename + "_aidiva_result_filtered.vcf"), bool(family_type == "SINGLE"))
    prioritized_data.to_csv(str(working_directory + input_filename + "_aidiva_result.tsv"), sep="\t", index=False)
    prioritized_data[prioritized_data["FILTER_PASSED"] == 1].to_csv(str(working_directory + input_filename + "_aidiva_result_filtered.tsv"), sep="\t", index=False)
    print("Pipeline successfully finsished!")
