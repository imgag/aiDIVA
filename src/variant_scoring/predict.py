from sklearn.ensemble import RandomForestClassifier
from sklearn.datasets import make_classification
from sklearn.model_selection import RandomizedSearchCV
from sklearn.model_selection import GridSearchCV
import pandas as pd
import numpy as np
from pprint import pprint
import argparse
import pickle


random_seed = 14038


def import_model(model_file):
    model_to_import = open(model_file, "rb")
    rf_model = pickle.load(model_to_import)
    model_to_import.close()
    
    return rf_model


def read_input_data(input_file):
    input_data = pd.read_csv(input_file, sep=',', low_memory=False)
    
    return input_data


# make sure to use the same features, that were used for the training of the model
# TODO change features to final feature set (or do it dynamically through config file)
def prepare_input_data(input_data):    
    # fill SegDup missing values with -> 0
    # fill ABB_SCORE missing values with -> 0
    input_data['SegDupMax'].fillna(0, inplace=True)
    input_data['ABB_SCORE'].fillna(0, inplace=True)
    
    # fill Allele Frequence missing values with -> 0
    input_data['gnomAD_AF'].fillna(0, inplace=True)
    input_data['gnomAD_AFR_AF'].fillna(0, inplace=True)
    input_data['gnomAD_AMR_AF'].fillna(0, inplace=True)
    input_data['gnomAD_ASJ_AF'].fillna(0, inplace=True)
    input_data['gnomAD_EAS_AF'].fillna(0, inplace=True)
    input_data['gnomAD_FIN_AF'].fillna(0, inplace=True)
    input_data['gnomAD_NFE_AF'].fillna(0, inplace=True)
    input_data['gnomAD_OTH_AF'].fillna(0, inplace=True)
    input_data['gnomAD_SAS_AF'].fillna(0, inplace=True)
    input_data['AF'].fillna(0, inplace=True)
    input_data['AFR_AF'].fillna(0, inplace=True)
    input_data['AMR_AF'].fillna(0, inplace=True)
    input_data['EAS_AF'].fillna(0, inplace=True)
    input_data['EUR_AF'].fillna(0, inplace=True)
    input_data['SAS_AF'].fillna(0, inplace=True)
    input_data['AA_AF'].fillna(0, inplace=True)
    input_data['EA_AF'].fillna(0, inplace=True)
    
    # if multiple allele frequencies are reported only take the maximum
    input_data["gnomAD_AF"] = input_data.apply(lambda row: pd.Series(max([float(i) for i in str(row["gnomAD_AF"]).split("&")])), axis=1)
    input_data["gnomAD_AFR_AF"] = input_data.apply(lambda row: pd.Series(max([float(i) for i in str(row["gnomAD_AFR_AF"]).split("&")])), axis=1)
    input_data["gnomAD_AMR_AF"] = input_data.apply(lambda row: pd.Series(max([float(i) for i in str(row["gnomAD_AMR_AF"]).split("&")])), axis=1)
    input_data["gnomAD_ASJ_AF"] = input_data.apply(lambda row: pd.Series(max([float(i) for i in str(row["gnomAD_ASJ_AF"]).split("&")])), axis=1)
    input_data["gnomAD_EAS_AF"] = input_data.apply(lambda row: pd.Series(max([float(i) for i in str(row["gnomAD_EAS_AF"]).split("&")])), axis=1)
    input_data["gnomAD_FIN_AF"] = input_data.apply(lambda row: pd.Series(max([float(i) for i in str(row["gnomAD_FIN_AF"]).split("&")])), axis=1)
    input_data["gnomAD_NFE_AF"] = input_data.apply(lambda row: pd.Series(max([float(i) for i in str(row["gnomAD_NFE_AF"]).split("&")])), axis=1)
    input_data["gnomAD_OTH_AF"] = input_data.apply(lambda row: pd.Series(max([float(i) for i in str(row["gnomAD_OTH_AF"]).split("&")])), axis=1)
    input_data["gnomAD_SAS_AF"] = input_data.apply(lambda row: pd.Series(max([float(i) for i in str(row["gnomAD_SAS_AF"]).split("&")])), axis=1)
    input_data["AF"] = input_data.apply(lambda row: pd.Series(max([float(i) for i in str(row["AF"]).split("&")])), axis=1)
    input_data["AFR_AF"] = input_data.apply(lambda row: pd.Series(max([float(i) for i in str(row["AFR_AF"]).split("&")])), axis=1)
    input_data["AMR_AF"] = input_data.apply(lambda row: pd.Series(max([float(i) for i in str(row["AMR_AF"]).split("&")])), axis=1)
    input_data["EAS_AF"] = input_data.apply(lambda row: pd.Series(max([float(i) for i in str(row["EAS_AF"]).split("&")])), axis=1)
    input_data["EUR_AF"] = input_data.apply(lambda row: pd.Series(max([float(i) for i in str(row["EUR_AF"]).split("&")])), axis=1)
    input_data["SAS_AF"] = input_data.apply(lambda row: pd.Series(max([float(i) for i in str(row["SAS_AF"]).split("&")])), axis=1)
    input_data["AA_AF"] = input_data.apply(lambda row: pd.Series(max([float(i) for i in str(row["AA_AF"]).split("&")])), axis=1)
    input_data["EA_AF"] = input_data.apply(lambda row: pd.Series(max([float(i) for i in str(row["EA_AF"]).split("&")])), axis=1)

    # compute maximum Minor Allele Frequency (MAF)
    input_data[["MaxAF"]] = input_data.apply(lambda row: pd.Series(max([float(row["AFR_AF"]), float(row["AMR_AF"]), float(row["EAS_AF"]), float(row["EUR_AF"]), float(row["SAS_AF"]), float(row["AA_AF"]), float(row["EA_AF"]), float(row["gnomAD_AFR_AF"]), float(row["gnomAD_AMR_AF"]), float(row["gnomAD_ASJ_AF"]), float(row["gnomAD_EAS_AF"]), float(row["gnomAD_FIN_AF"]), float(row["gnomAD_NFE_AF"]), float(row["gnomAD_OTH_AF"]), float(row["gnomAD_SAS_AF"])])), axis=1)

    # fill remaining missing values in remaining columns with the median of the respective column
    input_data['CADD_PHRED'].fillna(input_data['CADD_PHRED'].median(), inplace=True)
    input_data['Condel'].fillna(input_data['Condel'].median(), inplace=True)
    input_data['phyloP46_primate'].fillna(input_data['phyloP46_primate'].median(), inplace=True)
    input_data['phyloP46_mammal'].fillna(input_data['phyloP46_mammal'].median(), inplace=True)
    input_data['phastCons46_primate'].fillna(input_data['phastCons46_primate'].median(), inplace=True)
    input_data['phastCons46_mammal'].fillna(input_data['phastCons46_mammal'].median(), inplace=True)
    input_data['Eigen-phred'].fillna(input_data['Eigen-phred'].median(), inplace=True)
    input_data['MutationAssessor_score'].fillna(input_data['MutationAssessor_score'].median(), inplace=True)

    input_features = np.asarray(input_data[['CADD_PHRED','Condel','SegDupMax','phyloP46_primate','phyloP46_mammal', 'phastCons46_primate', 'phastCons46_mammal', 'Eigen-phred','MaxAF','MutationAssessor_score','ABB_SCORE']])

    return input_data, input_features


def predict_rank(rf_model_snps, rf_model_indel, input_data_snps, input_features_snps, input_data_indel, input_features_indel):
    class_prediction_snps = rf_model_snps.predict(input_features_snps)
    score_prediction_snps = pd.DataFrame(rf_model.predict_proba(input_features_snps), columns=["Probability_Benign", "Probability_Pathogenic"])

    class_prediction_indel = rf_model_indel.predict(input_features_indel)
    score_prediction_indel = pd.DataFrame(rf_model_indel.predict_proba(input_features_indel), columns=["Probability_Benign", "Probability_Pathogenic"])

    input_data_snps["Rank"] = score_prediction_snps["Probability_Pathogenic"]
    input_data_indel["Rank"] = score_prediction_indel["Probability_Pathogenic"]
    
    return input_data_snps, input_data_indel


def perform_pathogenicity_score_prediction(input_data_snps, input_data_indel, rf_model_snps, rf_model_indel):
    prepared_input_data_snps, input_features_snps = prepare_input_data(input_data_snps)
    prepared_input_data_indel, input_features_indel = prepare_input_data(input_data_indel)
    rf_model_snps = import_model(rf_model_snps)
    rf_model_indel = import_model(rf_model_indel)
    predicted_data_snps, predicted_data_indel = predict_rank(rf_model_snps, rf_model_indel, prepared_input_data_snps, input_features_snps, prepared_input_data_indel, input_features_indel)
    
    return predicted_data_snps, predicted_data_indel


if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--in_data', type=str, dest='in_data', metavar='in.csv', required=True, help='CSV file containing the training data, used to train the random forest model\n')
    parser.add_argument('--out_data', type=str, dest='out_data', metavar='out.csv', required=True, help='CSV file containing the test data, used to compute the model statistics\n')
    parser.add_argument('--model_snps', type=str, dest='model_snps', metavar='model_snps.pkl', required=True, help='Specifies the name of the trained snps model to import\n')
    parser.add_argument('--model_indel', type=str, dest='model_indel', metavar='model_indel.pkl', required=True, help='Specifies the name of the trained indel model to import\n')
    args = parser.parse_args()
    
    rf_model_snps = import_model(args.model_snps)
    rf_model_indel = import_model(args.model_indel)
    input_data = read_input_data(args.in_data)
    
    # if multiple alleles are reported consider only the first one
    # TODO decide how to handle allele ambiguity
    input_data["Alt"] = input_data["Alt"].map(lambda x: x.split(",")[0])
    
    input_data_snps = input_data[(input_data["Ref"].apply(len) == 1) & (input_data["Alt"].apply(len) == 1)]
    input_data_indel = input_data[(input_data["Ref"].apply(len) > 1) | (input_data["Alt"].apply(len) > 1)]
    
    #TODO add indel handling call functions from other script to expand and then call vep and afterwards combine
    
    predicted_data_snps, predicted_data_indel = perform_pathogenicity_score_prediction(input_data_snps, input_data_indel, rf_model_snps, rf_model_indel)
    
    predicted_data_combined = pd.concat([predicted_data_snps, predicted_data_indel])
    predicted_data_combined.sort_values(['Chr', 'Pos'], ascending=[True, True])
    
    predicted_data_combined.to_csv(args.out_data, index=False, sep="\t", na_rep="NA")
