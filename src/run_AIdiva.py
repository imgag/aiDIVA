import yaml
import argparse
import tempfile
import pandas as pd
import util.add_segmentDuplication as add_segDup
import util.add_simpleRepeats as add_repeats
import util.add_abb_score as add_abb
import util.add_score_from_bigwig as add_score
import variant_scoring.predict as predict
import util.createTrainingVCF as vcf
import util.convert_vcf_to_csv as convert_vcf
import variant_prioritization.familySNP_gene_score as prio



if __name__=='__main__':
    parser = argparse.ArgumentParser(description = 'AIdiva -- Augmented Intelligence Disease Variant Analysis')
    parser.add_argument('--config', type=str, dest='config', metavar='config.yaml', required=True, help='Config file specifying the parameters for AIdiva [required]')
    args = parser.parse_args()
    
    # parse configuration file
    config_file = open(args.config, "r")
    configuration = yaml.full_load(config_file)
    config_file.close()
    
    # parse input files
    input_vcf = configuration["Analysis-Input"]["vcf"] # should already be annotated with VEP
    scoring_model_snps = configuration["Analysis-Input"]["scoring-model-snps"]
    
    # parse output files
    output_tsv = configuration["Analysis-Output"]["tsv"]
    output_filtered_tsv = configuration["Analysis-Output"]["tsv-filt"]
    
    # parse disease and inheritance information
    hpo_file = configuration["Analysis-Input"]["inheritance-information"]["hpo-file"]
    gene_exclusion_file = configuration["Analysis-Input"]["inheritance-information"]["gene-exclusion"]
    inheritance = configuration["Analysis-Input"]["inheritance-information"]["inheritance"]
    family_type = configuration["Analysis-Input"]["inheritance-information"]["family-type"]
    family_file = configuration["Analysis-Input"]["inheritance-information"]["family-file"]
    
    # parse additional annotation files
    abb_score_file = configuration["Additional-Annotation"]["abb-score"]
    simple_repeat_file = configuration["Additional-Annotation"]["simple-repeats"]
    segment_duplication_file = configuration["Additional-Annotation"]["segment-duplication"]
    additional_bigwig_files = configuration["Additional-Annotation"]["additional-bigwig-files"]
    
    # convert input vcf to pandas dataframe
    input_data = convert_vcf.convert_vcf_to_pandas_dataframe(input_vcf)
    
    # add abb score to annotation
    input_data_additional_annotation = add_abb.group_and_process_data(abb_score_file, input_data)
    
    # add simple repeats to annotation
    input_data_additional_annotation = add_repeats.group_and_process_data(simple_repeat_file, input_data_additional_annotation)
    
    # add segment duplication to annotation
    input_data_additional_annotation = add_segDup.group_and_process_data(segment_duplication_file, input_data_additional_annotation)
    
    # add annotation from addtional bigwig files
    for key, value in additional_bigwig_files.items():
        input_data_additional_annotation = add_score.group_and_process_data(value, input_data_additional_annotation, key)
    
    input_data_additional_annotation.to_csv("temp2_additional_annotation_test_all_full.csv", index=False, sep="\t")
    
    #input_data_additional_annotation = pd.read_csv("temp2_additional_annotation_test_all_full.csv", low_memory=False, sep="\t")
    
    input_data_predicted = predict.perform_pathogenicity_score_prediction(input_data_additional_annotation, scoring_model_snps)
    
    #input_data_predicted.to_csv("tmp.predicted.tsv", index=False, sep="\t")
    
    tmp = tempfile.NamedTemporaryFile(mode="w")
    input_data_predicted.to_csv(tmp.name, index=False, sep="\t")
    
    #prio.main_program("tmp.predicted.tsv", output_tsv, output_filtered_tsv, family_file, inheritance, family_type, hpo_file, gene_exclusion_file)
    
    prio.main_program(tmp.name, output_tsv, output_filtered_tsv, family_file, inheritance, family_type, hpo_file, gene_exclusion_file)
    tmp.close()
    
    
    
    
    
