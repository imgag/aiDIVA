import argparse
import gzip
import logging
import multiprocessing as mp
import numpy as np
import pandas as pd
import pickle

from functools import partial



MEAN_DICT = {"phastCons46mammal": 0.09691308336428194,
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

MEDIAN_DICT = {"MutationAssessor": 1.87,
               "CONDEL": 0.4805749233199981,
               "EIGEN_PHRED": 3.010301,
               "CADD_PHRED": 3.99,
               "FATHMM_XF": 0.209614,
               "SIFT": 0.153,
               "REVEL": 0.193,
               "PolyPhen": 0.547,
               "oe_lof": 0.48225}

SUPPORTED_CODING_VARIANTS = ["stop_gained",
                             "stop_lost",
                             "start_lost",
                             "frameshift_variant",
                             "inframe_insertion",
                             "inframe_deletion",
                             "missense_variant",
                             "protein_altering_variant",
                             "incomplete_terminal_codon_variant",
                             "start_retained_variant",
                             "stop_retained_variant",
                             "coding_sequence_variant"]

random_seed = 14038
logger = logging.getLogger(__name__)


def import_model(model_file):
    if model_file.endswith(".gz"):
        rf_model = pickle.load(gzip.open(model_file, "rb"))
    else:
        rf_model = pickle.load(open(model_file, "rb"))

    return rf_model


def read_input_data(input_file):
    input_data = pd.read_csv(input_file, sep="\t", low_memory=False)

    return input_data


# make sure to use the same features, that were used for the training of the model
# fill SegDup missing values with -> 0
# fill ABB_SCORE missing values with -> 0
# fill Allele Frequency missing values with -> 0
# fill CAPICE missing values with -> 0.5
# fill missing values from other features with -> median or mean
def prepare_input_data(feature_list, allele_frequency_list, input_data):
    # compute maximum Minor Allele Frequency (MAF) from population frequencies if MAX_AF not present
    if ("MaxAF" not in input_data.columns) and ("MAX_AF" not in input_data.columns):
        if allele_frequency_list:
            for allele_frequency in allele_frequency_list:
                input_data[allele_frequency] = input_data[allele_frequency].fillna(0)
                input_data[allele_frequency] = input_data.apply(lambda row: pd.Series(max([float(frequency) for frequency in str(row[allele_frequency]).split("&")], default=np.nan)), axis=1)
            input_data["MaxAF"] = input_data.apply(lambda row: pd.Series(max([float(frequency) for frequency in row[allele_frequency_list].tolist()])), axis=1)
        else:
            logger.error("Empty allele frequency list was given!")

    for feature in feature_list:
        if (feature == "MaxAF") or (feature == "MAX_AF"):
            input_data[feature] = input_data[feature].fillna(0)
        
        elif feature == "homAF":
            input_data[feature] = input_data[feature].fillna(0)
        
        elif feature == "CAPICE":
            input_data[feature] = input_data[feature].fillna(0.5)
        
        elif "SIFT" == feature:
            input_data[feature] = input_data.apply(lambda row: min([float(value) for value in str(row[feature]).split("&") if ((value != ".") and (value != "nan") and (value != ""))], default=np.nan), axis=1)
            input_data[feature] = input_data[feature].fillna(MEDIAN_DICT["SIFT"])
        
        elif feature == "oe_lof":
            input_data[feature] = input_data.apply(lambda row: min([float(value) for value in str(row[feature]).split("&") if ((value != ".") and (value != "nan") and (value != "") and (":" not in value) and ("-" not in value))], default=np.nan), axis=1)
            input_data[feature] = input_data[feature].fillna(MEDIAN_DICT["oe_lof"])
        
        else:
            input_data[feature] = input_data.apply(lambda row: max([float(value) for value in str(row[feature]).split("&") if ((value != ".") and (value != "nan") and (value != ""))], default=np.nan), axis=1)
            if ("phastCons" in feature) or ("phyloP" in feature):
                input_data[feature] = input_data[feature].fillna(MEAN_DICT[feature])
            else:
                input_data[feature] = input_data[feature].fillna(MEAN_DICT[feature])

    return input_data


def predict_pathogenicity(rf_model, input_data, input_features):
    score_prediction = pd.DataFrame(rf_model.predict_proba(input_features), columns=["Probability_Benign", "Probability_Pathogenic"])
    input_data["AIDIVA_SCORE"] = score_prediction["Probability_Pathogenic"]

    return input_data


def parallel_pathogenicity_prediction(rf_model, input_data, input_features, num_cores):
    try:
        worker_pool = mp.Pool(num_cores)
        predicted_data = pd.concat(worker_pool.apply(rf_model.predict_proba(), (input_features)))

        input_data["AIDIVA_SCORE"] = predicted_data["Probability_Pathogenic"]
    finally:
        worker_pool.close()
        worker_pool.join()

    return input_data


def parallelize_dataframe_processing(dataframe, function, num_cores):
    num_partitions = num_cores * 2

    if len(dataframe) <= num_partitions:
        dataframe_splitted = np.array_split(dataframe, 1)
    else:
        dataframe_splitted = np.array_split(dataframe, num_partitions)

    try:
        pool = mp.Pool(num_cores)
        dataframe = pd.concat(pool.map(function, dataframe_splitted))
    finally:
        pool.close()
        pool.join()

    return dataframe


def perform_pathogenicity_score_prediction(rf_model, input_data, allele_frequency_list, feature_list, num_cores):
    rf_model = import_model(rf_model)
    prepared_input_data = parallelize_dataframe_processing(input_data, partial(prepare_input_data, feature_list, allele_frequency_list), num_cores)
    input_features = np.asarray(prepared_input_data[feature_list], dtype=np.float64)
    predicted_data = predict_pathogenicity(rf_model, input_data, input_features)

    # frameshift variants are not covered in the used model, set them to 0.9 (1.0 is too high)
    predicted_data.loc[((abs(predicted_data["REF"].str.len() - predicted_data["ALT"].str.len()) % 3 != 0)), "AIDIVA_SCORE"] = 0.9

    # set splicing donor/acceptor variants to NaN if not additionally a supported consequence is reported for the variant 
    # add filter for splice_region variants
    # later in the pipeline the score from dbscSNV for splicing variants will be used
    predicted_data.loc[(~(any(term for term in SUPPORTED_CODING_VARIANTS if term in predicted_data["Consequence"])) & (~(predicted_data["rf_score"].isna()) | ~(predicted_data["ada_score"].isna()))), "AIDIVA_SCORE"] = np.nan

    # set synonymous variants to 0.0 if they are not at a splicing site
    predicted_data.loc[(predicted_data["Consequence"].str.contains("synonymous") & ~(predicted_data["Consequence"].str.contains("splice"))), "AIDIVA_SCORE"] = 0.0

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
