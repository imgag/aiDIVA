import argparse
import gzip
import logging
import multiprocessing as mp
import numpy as np
import pandas as pd
import pickle

from functools import partial


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


# extract the most severe consequence if overlapping consequences were found
def get_most_severe_consequence(variant, VARIANT_CONSEQUENCES):
    consequences = str(variant["Consequence"])
    found_consequences = consequences.split("&")

    # use most severe consequence for filtering if overlapping consequences are present
    if len(found_consequences) > 1:
        most_severe_consequence = found_consequences[0]

        for consequence in found_consequences:
            if VARIANT_CONSEQUENCES[consequence] < VARIANT_CONSEQUENCES[most_severe_consequence]:
                most_severe_consequence = consequence

    else:
        most_severe_consequence = found_consequences[0]

    return most_severe_consequence


# make sure to use the same features, that were used for the training of the model
# fill Allele Frequency missing values with -> 0
# fill homAF missing values with -> 0
# fill HIGH_IMPACT missing values with -> 0
# fill missing values from other features with -> median or mean
def prepare_input_data(feature_list, allele_frequency_list, MEAN_DICT, MEDIAN_DICT, SPLICE_VARIANTS, input_data):
    for feature in feature_list:
        if feature == "MAX_AF":
            input_data[feature] = input_data[feature].fillna(0)

        elif feature == "homAF":
            input_data[feature] = input_data[feature].fillna(0)
            
        elif feature == "HIGH_IMPACT":
            input_data[feature] = input_data[feature].fillna(0)

        # TODO if there occur problems we have to impute missing values (probably with 0)
        elif (feature == "IS_INDEL"):
            continue

        elif "SIFT" == feature:
            input_data[feature] = input_data.apply(lambda row: min([float(value) for value in str(row[feature]).split("&") if ((value != ".") and (value != "nan") and (value != ""))], default=np.nan), axis=1)
            input_data[feature] = input_data[feature].fillna(MEDIAN_DICT["SIFT"])

        elif feature == "oe_lof" or feature == "oe_mis" or feature == "oe_syn":
            input_data[feature] = input_data.apply(lambda row: min([float(value) for value in str(row[feature]).split("&") if ((value != ".") and (value != "nan") and (value != "") and (":" not in value) and ("-" not in value))], default=np.nan), axis=1)
            input_data[feature] = input_data[feature].fillna(MEDIAN_DICT[feature])

        # REVEL and ALPHA_MISSENSE score only missense variants so for missing high impact variants we fill with 1 instead of median/mean value
        elif feature == "REVEL" or feature == "ALPHA_MISSENSE_SCORE":
            input_data[feature] = input_data.apply(lambda row: max([float(value) for value in str(row[feature]).split("&") if ((value != ".") & (value != "nan") & (value != "") & (not ":" in value) & (not "-" in value))], default=np.nan), axis=1)

            # check if information about the variant impact is present in the input table, otherwise fallback to fill missing values with median/mean
            if "IMPACT" in input_data.columns:
                input_data.loc[(~(input_data["MOST_SEVERE_CONSEQUENCE"].str.contains("|".join(SPLICE_VARIANTS))) & (input_data["IMPACT"] == "HIGH") & (input_data[feature].isna())), feature] = 1.0
                input_data[feature] = input_data[feature].fillna(MEDIAN_DICT[feature])

            elif "HIGH_IMPACT" in input_data.columns:
                input_data.loc[(~(input_data["MOST_SEVERE_CONSEQUENCE"].str.contains("|".join(SPLICE_VARIANTS))) & (input_data["HIGH_IMPACT"] == 1) & (input_data[feature].isna())), feature] = 1.0
                input_data[feature] = input_data[feature].fillna(MEDIAN_DICT[feature])

            else:
                input_data[feature] = input_data[feature].fillna(MEDIAN_DICT[feature])

        else:
            input_data[feature] = input_data.apply(lambda row: max([float(value) for value in str(row[feature]).split("&") if ((value != ".") and (value != "nan") and (value != ""))], default=np.nan), axis=1)
            if ("phastCons" in feature) or ("phyloP" in feature):
                input_data[feature] = input_data[feature].fillna(MEAN_DICT[feature])

            else:
                input_data[feature] = input_data[feature].fillna(MEDIAN_DICT[feature])

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
        # do not split dataframe
        dataframe_splitted = [dataframe]

    else:
        ## TODO: replace np.array_split() with iloc to prevent problems in future pandas versions
        # usage of floor division (//) makes sure that we get an absolute number as result
        #dataframe_splitted = np.array_split(dataframe, num_partitions) # -> uses deprecated functionality that will behave differently in future pandas versions
        chunk_size = dataframe.shape[0] // num_partitions
        dataframe_splitted = [dataframe[i:i+chunk_size].copy() for i in range(0, dataframe.shape[0], chunk_size)]

    try:
        pool = mp.Pool(num_cores)
        dataframe = pd.concat(pool.map(function, dataframe_splitted))

    finally:
        pool.close()
        pool.join()

    return dataframe


def perform_pathogenicity_score_prediction(rf_model, input_data, allele_frequency_list, feature_list, CONSTANT_DICTIONARY, num_cores=1):
    # get constants
    VARIANT_CONSEQUENCES = CONSTANT_DICTIONARY["VARIANT_CONSEQUENCES"]
    CODING_VARIANTS = CONSTANT_DICTIONARY["CODING_VARIANTS"]
    SPLICE_VARIANTS = ["SPLICE_VARIANTS"]
    SYNONYMOUS_VARIANTS = ["SYNONYMOUS_VARIANTS"]
    MEAN_DICT = CONSTANT_DICTIONARY["MEAN_DICT"]
    MEDIAN_DICT = CONSTANT_DICTIONARY["MEDIAN_DICT"]

    # import trained RF model
    rf_model = import_model(rf_model)

    # compute maximum Minor Allele Frequency (MAF) from population frequencies if MAX_AF not present
    if "MAX_AF" not in input_data.columns:
        if allele_frequency_list:
            for allele_frequency in allele_frequency_list:
                input_data[allele_frequency] = input_data[allele_frequency].fillna(0)
                input_data[allele_frequency] = input_data.apply(lambda row: pd.Series(max([float(frequency) for frequency in str(row[allele_frequency]).split("&")], default=np.nan)), axis=1)

            input_data["MAX_AF"] = input_data.apply(lambda row: pd.Series(max([float(frequency) for frequency in row[allele_frequency_list].tolist()], default=np.nan)), axis=1)

        else:
            logger.error("Empty allele frequency list was given!")

    if "MOST_SEVERE_CONSEQUENCE" not in input_data.columns:
        if "Consequence" in input_data.columns:
            input_data["MOST_SEVERE_CONSEQUENCE"] = input_data.apply(lambda row: pd.Series(get_most_severe_consequence(row, VARIANT_CONSEQUENCES)), axis=1)

        else:
            logger.warning("Could not determine MOST_SEVERE_CONSEQUENCE for the current variant use term 'unknown'!")
            input_data["MOST_SEVERE_CONSEQUENCE"] = "unknown"

    prepared_input_data = parallelize_dataframe_processing(input_data, partial(prepare_input_data, feature_list, allele_frequency_list, MEAN_DICT, MEDIAN_DICT, SPLICE_VARIANTS), num_cores)
    input_features = np.asarray(prepared_input_data[feature_list], dtype=np.float64)
    predicted_data = predict_pathogenicity(rf_model, input_data, input_features)

    # frameshift variants are not covered in the used model, set them to 0.9 (1.0 is too high)
    predicted_data.loc[((abs(predicted_data["REF"].str.len() - predicted_data["ALT"].str.len()) % 3 != 0)), "AIDIVA_SCORE"] = 0.9

    # set splicing donor/acceptor variants to NaN if not additionally a supported consequence is reported for the variant 
    # add filter for splice_region variants
    # later in the pipeline the score from SpliceAI for splicing variants will be used
    predicted_data.loc[(~(predicted_data["MOST_SEVERE_CONSEQUENCE"].str.contains("|".join(CODING_VARIANTS))) & (predicted_data["MOST_SEVERE_CONSEQUENCE"].str.contains("|".join(SPLICE_VARIANTS)))), "AIDIVA_SCORE"] = np.nan

    # set synonymous variants to 0.0 if they are not at a splicing site
    # TODO check how to handle overlapping consequences
    predicted_data.loc[(predicted_data["MOST_SEVERE_CONSEQUENCE"].str.contains("|".join(SYNONYMOUS_VARIANTS)) & ~(predicted_data["MOST_SEVERE_CONSEQUENCE"].str.contains("|".join(CODING_VARIANTS))) & ~(predicted_data["MOST_SEVERE_CONSEQUENCE"].str.contains("|".join(SPLICE_VARIANTS)))), "AIDIVA_SCORE"] = 0.0

    return predicted_data


if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--in_data", type=str, dest="in_data", metavar="in.csv", required=True, help="CSV file containing the training data, used to train the random forest model\n")
    parser.add_argument("--out_data", type=str, dest="out_data", metavar="out.csv", required=True, help="CSV file containing the test data, used to compute the model statistics\n")
    parser.add_argument("--model", type=str, dest="model", metavar="model.pkl", required=True, help="Specifies the name of the trained model to import\n")
    parser.add_argument("--feature_list", type=str, dest="feature_list", metavar="feature1,feature2,feature3", required=True, help="Comma separated list of the features used to train the model\n")
    parser.add_argument("--allele_frequency_list", type=str, dest="allele_frequency_list", metavar="frequency1,frequecy2,frequency3", required=False, help="Comma separated list of allele frequency sources that should be used as basis to get the maximum allele frequency\n")
    parser.add_argument("--threads", type=str, dest="threads", metavar="1", required=False, help="Number of threads to use.\n")
    args = parser.parse_args()

    input_data = read_input_data(args.in_data)
    feature_list = args.feature_list.split(",")

    if args.allele_frequency_list:
        allele_frequency_list = args.allele_frequency_list.split(",")

    else:
        allele_frequency_list = []

    if args.threads:
        num_threads = int(args.threads)

    else:
        num_threads = 1

    predicted_data = perform_pathogenicity_score_prediction(args.model, input_data, allele_frequency_list, feature_list, num_threads)
    predicted_data.to_csv(args.out_data, index=False, sep="\t", na_rep="NA")
