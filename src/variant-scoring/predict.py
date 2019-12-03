from sklearn.ensemble import RandomForestClassifier
from sklearn.datasets import make_classification
from sklearn.model_selection import RandomizedSearchCV
from sklearn.model_selection import GridSearchCV
#from sklearn.model_selection import input_test_split
import pandas as pd
import numpy as np
from pprint import pprint
import pickle
import sys

random_seed = 14038

model_to_import = open(sys.argv[1], "rb")
rf_model = pickle.load(model_to_import)
model_to_import.close()

input_data = pd.read_csv(sys.argv[2], sep='\t', low_memory=False)

input_data[["SegDupMax"]] = input_data.apply(lambda row: pd.Series(max([float(i) for i in str(row["SegmentDuplication"]).split("&")])), axis=1)
input_data["AF"] = input_data.apply(lambda row: pd.Series(max([float(i) for i in str(row["AF"]).split("&")])), axis=1)
input_data[["MaxAF"]] = input_data.apply(lambda row: pd.Series(np.max([row["AF"], row["gnomAD_AF"]])), axis=1)

# fill SegDup missing values with -> 0
# fill ABB_SCORE missing values with -> 1
# fill Allele Frequence missing values with -> 0
input_data['SegDupMax'].fillna(0, inplace=True)
input_data['ABB_SCORE'].fillna(1, inplace=True)
input_data['MaxAF'].fillna(0, inplace=True) # remaining AF columns could also be filled with zeros
input_data['AF'].fillna(0, inplace=True)
input_data['gnomAD_AF'].fillna(0, inplace=True)

#input_data[["Condel_val"]] = input_data["Condel"].str.extract(r'\((.*)\)')

# fill remaining missing values in remaining columns with the median of the respective column
input_data['CADD_PHRED'].fillna(input_data['CADD_PHRED'].median(), inplace=True)
#input_data['Condel_val'].fillna(input_data['Condel_val'].median(), inplace=True)
input_data['Condel'].fillna(input_data['Condel'].median(), inplace=True)
input_data['phyloP46_primates'].fillna(input_data['phyloP46_primates'].median(), inplace=True)
input_data['phyloP46_mammal'].fillna(input_data['phyloP46_mammal'].median(), inplace=True)
input_data['phastCons46_primates'].fillna(input_data['phastCons46_primates'].median(), inplace=True)
input_data['phastCons46_mammal'].fillna(input_data['phastCons46_mammal'].median(), inplace=True)
input_data['Eigen-phred'].fillna(input_data['Eigen-phred'].median(), inplace=True)
input_data['MutationAssessor_score'].fillna(input_data['MutationAssessor_score'].median(), inplace=True)


input_features = np.asarray(input_data[['CADD_PHRED','Condel','SegDupMax','phyloP46_primates','phyloP46_mammal', 'phastCons46_primates', 'phastCons46_mammal', 'Eigen-phred','MaxAF','MutationAssessor_score','ABB_SCORE']])

class_prediction = rf_model.predict(input_features)
score_prediction = pd.DataFrame(rf_model.predict_proba(input_features), columns=["Probability_Benign", "Probability_Pathogenic"])

input_data["Rank"] = score_prediction["Probability_Pathogenic"]

#input_data.to_csv(sys.argv[3], index=False, sep="\t")

essential_data = input_data[['Chr', 'Pos', 'Ref', 'Alt', 'QUAL', 'FILTER', 'Consequence', 'SYMBOL', 'AF', 'gnomAD_AF', 'MaxAF', 'SegDupMax', 'phyloP46_mammal', 'phyloP46_primates', 'phastCons46_mammal', 'phastCons46_primates', 'MutationAssessor_score', 'Condel', 'CADD_PHRED', 'Eigen-phred', 'ABB_SCORE', 'SimpleTandemRepeatLength', 'SimpleTandemRepeatRegion', 'NA12878', 'DP.NA12878', 'REF.NA12878', 'ALT.NA12878', 'AF.NA12878', 'GQ.NA12878', 'NA12891', 'DP.NA12891', 'REF.NA12891', 'ALT.NA12891', 'AF.NA12891', 'GQ.NA12891', 'NA12892', 'DP.NA12892', 'REF.NA12892', 'ALT.NA12892', 'AF.NA12892', 'GQ.NA12892', 'Rank']]

essential_data.to_csv(sys.argv[3], index=False, sep="\t", na_rep="NA")
