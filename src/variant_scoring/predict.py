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
    #input_data[["SegDupMax"]] = input_data.apply(lambda row: pd.Series(max([float(i) for i in str(row["SegmentDuplication"]).split("&")])), axis=1)
    #input_data["AF"] = input_data.apply(lambda row: pd.Series(max([float(i) for i in str(row["AF"]).split("&")])), axis=1)
    #input_data[["MaxAF"]] = input_data.apply(lambda row: pd.Series(np.max([row["AF"], row["gnomAD_AF"]])), axis=1)
    
    # fill SegDup missing values with -> 0
    # fill ABB_SCORE missing values with -> 0
    # fill Allele Frequence missing values with -> 0
    input_data['SegDupMax'].fillna(0, inplace=True)
    input_data['ABB_SCORE'].fillna(1, inplace=True)
    #input_data['MaxAF'].fillna(0, inplace=True) # remaining AF columns could also be filled with zeros
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

    input_data[["MaxAF"]] = input_data.apply(lambda row: pd.Series(max([float(row["AFR_AF"]), float(row["AMR_AF"]), float(row["EAS_AF"]), float(row["EUR_AF"]), float(row["SAS_AF"]), float(row["AA_AF"]), float(row["EA_AF"]), float(row["gnomAD_AFR_AF"]), float(row["gnomAD_AMR_AF"]), float(row["gnomAD_ASJ_AF"]), float(row["gnomAD_EAS_AF"]), float(row["gnomAD_FIN_AF"]), float(row["gnomAD_NFE_AF"]), float(row["gnomAD_OTH_AF"]), float(row["gnomAD_SAS_AF"])])), axis=1)

    #input_data[["Condel_val"]] = input_data["Condel"].str.extract(r'\((.*)\)')

    # fill remaining missing values in remaining columns with the median of the respective column        
    input_data['CADD_PHRED'].fillna(input_data['CADD_PHRED'].median(), inplace=True)
    #input_data['Condel_val'].fillna(input_data['Condel_val'].median(), inplace=True)
    input_data['Condel'].fillna(input_data['Condel'].median(), inplace=True)
    input_data['phyloP46_primate'].fillna(input_data['phyloP46_primate'].median(), inplace=True)
    input_data['phyloP46_mammal'].fillna(input_data['phyloP46_mammal'].median(), inplace=True)
    input_data['phastCons46_primate'].fillna(input_data['phastCons46_primate'].median(), inplace=True)
    input_data['phastCons46_mammal'].fillna(input_data['phastCons46_mammal'].median(), inplace=True)
    input_data['Eigen-phred'].fillna(input_data['Eigen-phred'].median(), inplace=True)
    input_data['MutationAssessor_score'].fillna(input_data['MutationAssessor_score'].median(), inplace=True)


    print(input_data[['CADD_PHRED','Condel','SegDupMax','phyloP46_primate','phyloP46_mammal', 'phastCons46_primate', 'phastCons46_mammal', 'Eigen-phred','MaxAF','MutationAssessor_score','ABB_SCORE']].isna().sum())
    
    input_features = np.asarray(input_data[['CADD_PHRED','Condel','SegDupMax','phyloP46_primate','phyloP46_mammal', 'phastCons46_primate', 'phastCons46_mammal', 'Eigen-phred','MaxAF','MutationAssessor_score','ABB_SCORE']])
    
    #input_data['Cadd2'].fillna(0, inplace=True)
    #input_data['Condel'].fillna(0, inplace=True)
    #input_data['SegMentDup'].fillna(0, inplace=True)
    #input_data['PrimatesPhyloP'].fillna(0, inplace=True)
    #input_data['PlacentalMammalPhyloP'].fillna(0, inplace=True)
    #input_data['PrimatesPhastCons'].fillna(0, inplace=True)
    #input_data['PlacentalMammalPhastCons'].fillna(0, inplace=True)
    #input_data['Eigen_Phred'].fillna(0, inplace=True)
    #input_data['ExAC_AF'].fillna(0, inplace=True)
    #input_data['MutAss'].fillna(input_data['MutAss'].median(), inplace=True)
    #input_data['ABB_score'].fillna(0, inplace=True)

    #input_features = np.asarray(input_data[['Cadd2','Condel','SegMentDup','PrimatesPhyloP','PlacentalMammalPhyloP', 'PrimatesPhastCons', 'PlacentalMammalPhastCons', 'Eigen_Phred','ExAC_AF','MutAss','ABB_score']])
    
    return input_data, input_features


def predict_rank(rf_model, input_data, input_features):
    class_prediction = rf_model.predict(input_features)
    score_prediction = pd.DataFrame(rf_model.predict_proba(input_features), columns=["Probability_Benign", "Probability_Pathogenic"])

    input_data["Rank"] = score_prediction["Probability_Pathogenic"]
    
    return input_data


def perform_pathogenicity_score_prediction(input_data, rf_model):
    prepared_input_data, input_features = prepare_input_data(input_data)
    rf_model_imported = import_model(rf_model)
    predicted_data = predict_rank(rf_model_imported, prepared_input_data, input_features)
    
    return predicted_data


if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--in_data', type=str, dest='in_data', metavar='in.csv', required=True, help='CSV file containing the training data, used to train the random forest model\n')
    parser.add_argument('--out_data', type=str, dest='out_data', metavar='out.csv', required=True, help='CSV file containing the test data, used to compute the model statistics\n')
    parser.add_argument('--model', type=str, dest='model', metavar='model.pkl', required=True, help='Specifies the name of the trained model to import\n')
    args = parser.parse_args()
    
    rf_model = import_model(args.model)
    input_data = read_input_data(args.in_data)
    prepared_input_data, input_features = prepare_input_data(input_data)
    predicted_data = predict_rank(rf_model, prepared_input_data, input_features)
    predicted_data.to_csv(args.out_data, index=False, sep="\t", na_rep="NA")
    
    #essential_data = input_data[['Chr', 'Pos', 'Ref', 'Alt', 'QUAL', 'FILTER', 'Consequence', 'SYMBOL', 'AF', 'gnomAD_AF', 'MaxAF', 'SegDupMax', 'phyloP46_mammal', 'phyloP46_primates', 'phastCons46_mammal', 'phastCons46_primates', 'MutationAssessor_score', 'Condel', 'CADD_PHRED', 'Eigen-phred', 'ABB_SCORE', 'SimpleTandemRepeatLength', 'SimpleTandemRepeatRegion', 'NA12878', 'DP.NA12878', 'REF.NA12878', 'ALT.NA12878', 'AF.NA12878', 'GQ.NA12878', 'NA12891', 'DP.NA12891', 'REF.NA12891', 'ALT.NA12891', 'AF.NA12891', 'GQ.NA12891', 'NA12892', 'DP.NA12892', 'REF.NA12892', 'ALT.NA12892', 'AF.NA12892', 'GQ.NA12892', 'Rank']]
