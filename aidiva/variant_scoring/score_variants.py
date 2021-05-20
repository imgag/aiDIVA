from sklearn.ensemble import RandomForestClassifier
from sklearn.datasets import make_classification
from sklearn.model_selection import RandomizedSearchCV
from sklearn.model_selection import GridSearchCV
import pandas as pd
import numpy as np
import multiprocessing as mp
import argparse
import pickle
import time


mean_dict = {"phastCons46mammal": 0.09691308336428194,
             "phastCons46primate": 0.12353343703613741,
             "phastCons46vertebrate": 0.1366339183101041,
             "phyloP46mammal": -0.0063575303590607925,
             "phyloP46primate": -0.012076641890840553,
             "phyloP46vertebrate": 0.06761867323083483,
             "phastCons100": 0.11273633387190414,
             "phyloP100": 0.052907788505469275,
             "MutationAssessor": 1.7961304794577417,
             "CONDEL": 0.49699016949707825,
             "EIGEN_PHRED": 4.342947928406315,
             "CADD_PHRED": 4.471745325,
             "FATHMM_XF": 0.35846023623584666,
             "SIFT": 0.35216996259535444,
             "REVEL": 0.28019263637740743,
             "PolyPhen": 0.5169017014355943,
             "oe_lof": 0.53667483333333332}

median_dict = {"MutationAssessor": 1.87,
               "CONDEL": 0.4805749233199981,
               "EIGEN_PHRED": 3.010301,
               "CADD_PHRED": 3.99,
               "FATHMM_XF": 0.209614,
               "SIFT": 0.153,
               "REVEL": 0.193,
               "PolyPhen": 0.547,
               "oe_lof": 0.48225}

random_seed = 14038

coding_region = None
rf_model_snp = None
rf_model_indel = None
allele_frequencies = []
features = []


def import_model(model_file):
    model_to_import = open(model_file, "rb")
    rf_model = pickle.load(model_to_import)
    model_to_import.close()

    return rf_model


def read_input_data(input_file):
    input_data = pd.read_csv(input_file, sep="\t", low_memory=False)

    return input_data


# make sure to use the same features, that were used for the training of the model
# fill SegDup missing values with -> 0
# fill ABB_SCORE missing values with -> 0
# fill Allele Frequency missing values with -> 0
# fill missing values from other features with -> median or mean
def prepare_input_data(input_data):
    # compute maximum Minor Allele Frequency (MAF) from population frequencies if not present
    if ("MaxAF" not in input_data.columns) and ("MAX_AF" not in input_data.columns):
        if allele_frequencies:
            for allele_frequency in allele_frequencies:
                input_data[allele_frequency] = input_data[allele_frequency].fillna(0)
                input_data[allele_frequency] = input_data.apply(lambda row: pd.Series(max([float(frequency) for frequency in str(row[allele_frequency]).split("&")], default=np.nan)), axis=1)
            input_data["MaxAF"] = input_data.apply(lambda row: pd.Series(max([float(frequency) for frequency in row[allele_frequency_list].tolist()])), axis=1)
        else:
            print("ERROR: Empty allele frequency list was given!")

    for feature in features:
        if (feature == "MaxAF") or (feature == "MAX_AF"):
            input_data[feature] = input_data[feature].fillna(0)
        elif feature == "homAF":
            input_data[feature] = input_data[feature].fillna(0)
        elif feature == "segmentDuplication":
            input_data[feature] = input_data.apply(lambda row: max([float(value) for value in str(row[feature]).split("&") if ((value != ".") and (value != "nan") and (value != ""))], default=np.nan), axis=1)
            input_data[feature] = input_data[feature].fillna(0)
        elif feature == "ABB_SCORE":
            input_data[feature] = input_data[feature].fillna(0)
        elif feature == "CAPICE":
            input_data[feature] = input_data[feature].fillna(0.5)
        elif "SIFT" == feature:
            input_data[feature] = input_data.apply(lambda row: min([float(value) for value in str(row[feature]).split("&") if ((value != ".") and (value != "nan") and (value != ""))], default=np.nan), axis=1)
            input_data[feature] = input_data[feature].fillna(median_dict["SIFT"])
        elif feature == "oe_lof":
            input_data[feature] = input_data.apply(lambda row: min([float(value) for value in str(row[feature]).split("&") if ((value != ".") and (value != "nan") and (value != "") and (":" not in value) and ("-" not in value))], default=np.nan), axis=1)
            input_data[feature] = input_data[feature].fillna(median_dict["oe_lof"])
        else:
            input_data[feature] = input_data.apply(lambda row: max([float(value) for value in str(row[feature]).split("&") if ((value != ".") and (value != "nan") and (value != ""))], default=np.nan), axis=1)
            if ("phastCons" in feature) or ("phyloP" in feature):
                input_data[feature] = input_data[feature].fillna(mean_dict[feature])
            else:
                input_data[feature] = input_data[feature].fillna(median_dict[feature])

    return input_data


def predict_pathogenicity(rf_model, input_data, input_features):
    score_prediction = pd.DataFrame(rf_model.predict_proba(input_features), columns=["Probability_Benign", "Probability_Pathogenic"])
    input_data["AIDIVA_SCORE"] = score_prediction["Probability_Pathogenic"]

    return input_data


def parallel_pathogenicity_prediction(rf_model, input_data, input_features, num_cores):
    splitted_input_data = np.array_split(input_data, num_cores)

    worker_pool = mp.Pool(num_cores)
    predicted_data = pd.concat(worker_pool.apply(rf_model.predict_proba(), (input_features)))

    input_data["AIDIVA_SCORE"] = predicted_data["Probability_Pathogenic"]

    worker_pool.close()
    worker_pool.join()

    return input_data


def parallelize_dataframe_processing(dataframe, function, num_cores):
    num_partitions = num_cores * 2

    if len(dataframe) <= num_partitions:
        dataframe_splitted = np.array_split(dataframe, 1)
    else:
        dataframe_splitted = np.array_split(dataframe, num_partitions)

    pool = mp.Pool(num_cores)
    dataframe = pd.concat(pool.map(function, dataframe_splitted))
    pool.close()
    pool.join()

    return dataframe


def perform_pathogenicity_score_prediction(rf_model, input_data, allele_frequency_list, feature_list, num_cores):
    global features
    features = feature_list

    global allele_frequencies
    allele_frequencies = allele_frequency_list

    rf_model = import_model(rf_model)
    prepared_input_data = parallelize_dataframe_processing(input_data, prepare_input_data, num_cores)
    input_features = np.asarray(prepared_input_data[features], dtype=np.float64)
    predicted_data = predict_pathogenicity(rf_model, prepared_input_data, input_features)

    # frameshift variants are not covered in the used model, set them to 1.0 if the MAX_AF is less or equal than 0.02
    predicted_data.loc[((abs(predicted_data["REF"].str.len() - predicted_data["ALT"].str.len()) % 3 != 0)), "AIDIVA_SCORE"] = np.nan # could also be set to 1.0

    # set splicing donor/acceptor variants to NaN or 1.0
    predicted_data.loc[(predicted_data["Consequence"].str.contains("splice_acceptor_variant") | predicted_data["Consequence"].str.contains("splice_donor_variant")), "AIDIVA_SCORE"] = np.nan

    # set synonymous variants to NaN (could also be set to 0.0)
    predicted_data.loc[(predicted_data["Consequence"].str.contains("synonymous")), "AIDIVA_SCORE"] = 0.0

    # exclude chromosomes MT
    ## TODO: can be removed (mitochondrial variants should already be filtered out)
    predicted_data.loc[((predicted_data["CHROM"] == "MT") | (predicted_data["CHROM"] == "chrMT") | (predicted_data["CHROM"] == "M") | (predicted_data["CHROM"] == "chrM")), "AIDIVA_SCORE"] = np.nan

    return predicted_data


if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--in_data", type=str, dest="in_data", metavar="in.csv", required=True, help="CSV file containing the training data, used to train the random forest model\n")
    parser.add_argument("--out_data", type=str, dest="out_data", metavar="out.csv", required=True, help="CSV file containing the test data, used to compute the model statistics\n")
    parser.add_argument("--model", type=str, dest="model", metavar="model.pkl", required=True, help="Specifies the name of the trained model to import\n")
    parser.add_argument("--feature_list", type=str, dest="feature_list", metavar="feature1,feature2,feature3", required=True, help="Comma separated list of the features used to train the model\n")
    parser.add_argument("--allele_frequency_list", type=str, dest="allele_frequency_list", metavar="frequency1,frequecy2,frequency3", required=False, help="Comma separated list of allele frequency sources that should be used as basis to get the maximum allele frequency\n")
    args = parser.parse_args()

    input_data = read_input_data(args.in_data)
    feature_list = args.feature_list.split(",")

    if args.allele_frequency_list:
        allele_frequency_list = args.allele_frequency_list.split(",")
    else:
        allele_frequency_list = []

    predicted_data = perform_pathogenicity_score_prediction(args.model, input_data, allele_frequency_list, feature_list)
    predicted_data.to_csv(args.out_data, index=False, sep="\t", na_rep="NA")
