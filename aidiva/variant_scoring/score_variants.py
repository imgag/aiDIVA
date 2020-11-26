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

coding_variants = ["splice_acceptor_variant",
                   "splice_donor_variant",
                   "stop_gained",
                   "frameshift_variant",
                   "stop_lost",
                   "start_lost",
                   "inframe_insertion",
                   "inframe_deletion",
                   "missense_variant",
                   "protein_altering_variant",
                   "splice_region_variant",
                   "incomplete_terminal_codon_variant",
                   "start_retained_variant",
                   "stop_retained_variant",
                   "synonymous_variant",
                   "coding_sequence_variant",
                   "5_prime_UTR_variant",
                   "3_prime_UTR_variant"]

random_seed = 14038
num_partitions = 10
num_cores = 5

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
def prepare_input_data(input_data):
    # fill SegDup missing values with -> 0
    # fill ABB_SCORE missing values with -> 0
    # fill Allele Frequence missing values with -> 0
    # fill missing values from other features with -> median or mean

    #for allele_frequency in allele_frequencies:
    #    input_data[allele_frequency] = input_data[allele_frequency].fillna(0)
    #    input_data[allele_frequency] = input_data.apply(lambda row: pd.Series(max([float(frequency) for frequency in str(row[allele_frequency]).split("&")], default=np.nan)), axis=1)

    # compute maximum Minor Allele Frequency (MAF)
    #if (not "MAX_AF" in input_data.columns) & (not "MaxAF" in input_data.columns):
    #    input_data["MaxAF"] = input_data.apply(lambda row: pd.Series(max([float(frequency) for frequency in row[allele_frequency_list].tolist()])), axis=1)

    global features

    for feature in features:
        if feature == "MaxAF" or feature == "MAX_AF":
            input_data[feature] = input_data[feature].fillna(0)
        elif feature == "homAF":
            input_data[feature] = input_data[feature].fillna(0)
        elif feature == "segmentDuplication":
            input_data[feature] = input_data.apply(lambda row: max([float(value) for value in str(row[feature]).split("&") if ((value != ".") & (value != "nan") & (value != ""))], default=np.nan), axis=1)
            input_data[feature] = input_data[feature].fillna(0)
        elif feature == "ABB_SCORE":
            input_data[feature] = input_data[feature].fillna(0)
        elif "SIFT" == feature:
            input_data[feature] = input_data.apply(lambda row: min([float(value) for value in str(row[feature]).split("&") if ((value != ".") & (value != "nan") & (value != ""))], default=np.nan), axis=1)
            input_data[feature] = input_data[feature].fillna(median_dict["SIFT"])
        elif feature == "oe_lof":
            input_data[feature] = input_data.apply(lambda row: min([float(value) for value in str(row[feature]).split("&") if ((value != ".") & (value != "nan") & (value != "") & (not ":" in value) & (not "-" in value))], default=np.nan), axis=1)
            input_data[feature] = input_data[feature].fillna(median_dict["oe_lof"])
        else:
            input_data[feature] = input_data.apply(lambda row: max([float(value) for value in str(row[feature]).split("&") if ((value != ".") & (value != "nan") & (value != ""))], default=np.nan), axis=1)
            if ("phastCons" in feature) | ("phyloP" in feature):
                input_data[feature] = input_data[feature].fillna(mean_dict[feature])
            else:
                input_data[feature] = input_data[feature].fillna(median_dict[feature])

    # TODO add workaround to handle the rare case that for one of the features only NaNs are present and therefor the median leads also to a nan
    # in that case the input features contains NaNs and lead to an error

    return input_data


def predict_pathogenicity(rf_model_snp, rf_model_indel, input_data_snp, input_features_snp, input_data_indel, input_features_indel):
    score_prediction_snp = pd.DataFrame(rf_model_snp.predict_proba(input_features_snp), columns=["Probability_Benign", "Probability_Pathogenic"])
    score_prediction_indel = pd.DataFrame(rf_model_indel.predict_proba(input_features_indel), columns=["Probability_Benign", "Probability_Pathogenic"])

    input_data_snp["AIDIVA_SCORE"] = score_prediction_snp["Probability_Pathogenic"]
    input_data_indel["AIDIVA_SCORE"] = score_prediction_indel["Probability_Pathogenic"]

    return input_data_snp, input_data_indel


def parallel_pathogenicity_prediction(rf_model, input_data, input_features, n_cores):
    if n_cores is None:
        num_cores = 1
    else:
        num_cores = n_cores
    splitted_input_data = np.array_split(input_data, num_cores)

    worker_pool = mp.Pool(num_cores)
    predicted_data = pd.concat(worker_pool.apply(rf_model.predict_proba(), (input_features)))

    input_data["AIDIVA_SCORE"] = predicted_data["Probability_Pathogenic"]

    worker_pool.close()
    worker_pool.join()

    return input_data


def check_coding(variant_to_check):
    if any(term for term in coding_variants if term in variant_to_check["Consequence"]):
        return 1
    else:
        return 0


def parallelized_coding_check(input_data):
    input_data["CODING"] = input_data.apply(lambda row: check_coding(row), axis=1)

    return input_data


def parallelize_dataframe_processing(dataframe, function):
    global num_partitions
    global num_cores
    if num_cores is None:
        num_cores = 1
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


## TODO: handle the case if one of the dataframes is empty
def perform_pathogenicity_score_prediction(input_data_snp, input_data_indel, rf_model_snp, rf_model_indel, allele_frequency_list, feature_list, n_cores=1):
    t = time.time()

    global features
    features = feature_list

    global allele_frequencies
    allele_frequencies = allele_frequency_list

    global num_cores
    num_cores = n_cores

    rf_model_snp = import_model(rf_model_snp)
    rf_model_indel = import_model(rf_model_indel)
    print("Finished model import in: %.2f seconds" % (time.time() - t))

    # filter variants, only coding variants should be scored

    t = time.time()
    input_data_snp = parallelize_dataframe_processing(input_data_snp, parallelized_coding_check)
    input_data_indel = parallelize_dataframe_processing(input_data_indel, parallelized_coding_check)
    print("Finished parallelized coding check in: %.2f seconds" % (time.time() - t))

    t = time.time()
    prepared_input_data_snp = parallelize_dataframe_processing(input_data_snp, prepare_input_data)
    prepared_input_data_indel = parallelize_dataframe_processing(input_data_indel, prepare_input_data)
    input_features_snp = np.asarray(prepared_input_data_snp[features], dtype=np.float64)
    input_features_indel = np.asarray(prepared_input_data_indel[features], dtype=np.float64)
    print("Finished data preparation in: %.2f seconds" % (time.time() - t))

    t = time.time()
    print(prepared_input_data_indel.dtypes)
    predicted_data_snp, predicted_data_indel = predict_pathogenicity(rf_model_snp, rf_model_indel, prepared_input_data_snp, input_features_snp, prepared_input_data_indel, input_features_indel)
    print(predicted_data_indel["MAX_AF"].unique())
    print(predicted_data_indel.dtypes)
    print(sum(predicted_data_indel["MAX_AF"].isna()))

    predicted_data_snp.to_csv("/mnt/storage1/users/ahboced1/test/snp_predicted.csv", sep="\t", index=False)
    predicted_data_indel.to_csv("/mnt/storage1/users/ahboced1/test/indel_predicted.csv", sep="\t", index=False)
    print("Finished prediction in: %.2f seconds" % (time.time() - t))

    t = time.time()
    # frameshift variants are not covered in the used model, set them to 1.0 if the MAX_AF is less or equal than 0.02
    # the following line might produce an SettingWithCopyWarning this Warning should be a false positive in this case
    predicted_data_indel.loc[((abs(predicted_data_indel["REF"].str.len() - predicted_data_indel["ALT"].str.len()) % 3 != 0)), "AIDIVA_SCORE"] = 1.0 #np.nan # could also be set to 1.0
    predicted_data_indel.loc[np.greater_equal(pd.to_numeric(prepared_input_data_indel["MAX_AF"]), 0.01), "AIDIVA_SCORE"] = np.nan

    # set splicing donor/acceptor variants to NaN or 1.0
    predicted_data_snp.loc[(predicted_data_snp["Consequence"].str.contains("splice_acceptor_variant") | predicted_data_snp["Consequence"].str.contains("splice_donor_variant")), "AIDIVA_SCORE"] = 1.0 #np.nan
    predicted_data_indel.loc[(predicted_data_indel["Consequence"].str.contains("splice_acceptor_variant") | predicted_data_indel["Consequence"].str.contains("splice_donor_variant")), "AIDIVA_SCORE"] = 1.0 #np.nan

    # set synonymous variants to NaN (could also be set to 0.0)
    predicted_data_snp.loc[(predicted_data_snp["Consequence"].str.contains("synonymous")), "AIDIVA_SCORE"] = 0.0 #np.nan
    predicted_data_indel.loc[(predicted_data_indel["Consequence"].str.contains("synonymous")), "AIDIVA_SCORE"] = 0.0 #np.nan

    # set score for non-coding variants to NaN
    # the models are only for coding variants
    predicted_data_snp.loc[((predicted_data_snp.CODING == 0)), "AIDIVA_SCORE"] = np.nan
    predicted_data_indel.loc[((predicted_data_indel.CODING == 0)), "AIDIVA_SCORE"] = np.nan

    # exclude chromosomes Y and MT
    #predicted_data_snp.loc[((predicted_data_snp["CHROM"] == "Y") | (predicted_data_snp["CHROM"] == "chrY")), "AIDIVA_SCORE"] = np.nan
    #predicted_data_indel.loc[((predicted_data_indel["CHROM"] == "Y") | (predicted_data_indel["CHROM"] == "chrY")), "AIDIVA_SCORE"] = np.nan
    predicted_data_snp.loc[((predicted_data_snp["CHROM"] == "MT") | (predicted_data_snp["CHROM"] == "chrMT") | (predicted_data_snp["CHROM"] == "M") | (predicted_data_snp["CHROM"] == "chrM")), "AIDIVA_SCORE"] = np.nan
    predicted_data_indel.loc[((predicted_data_indel["CHROM"] == "MT") | (predicted_data_indel["CHROM"] == "chrMT") | (predicted_data_indel["CHROM"] == "M") | (predicted_data_indel["CHROM"] == "chrM")), "AIDIVA_SCORE"] = np.nan
    print("Finished score adjustments in: %.2f seconds" % (time.time() - t))

    # combine snp and indel data
    ## TODO: the sort attribute is not present in pandas v0.19 so if this version should be
    #predicted_data_complete = pd.concat([predicted_data_snp, predicted_data_indel], sort=False)
    predicted_data_complete = pd.concat([predicted_data_snp, predicted_data_indel])
    predicted_data_complete.sort_values(["CHROM", "POS"], ascending=[True, True], inplace=True)
    predicted_data_complete.reset_index(inplace=True, drop=True)
    predicted_data_complete = predicted_data_complete[predicted_data_snp.columns]

    return predicted_data_complete


if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--in_data_snp", type=str, dest="in_data_snp", metavar="in.csv", required=True, help="CSV file containing the training data, used to train the random forest model\n")
    parser.add_argument("--in_data_indel", type=str, dest="in_data_indel", metavar="in.csv", required=True, help="CSV file containing the training data, used to train the random forest model\n")
    parser.add_argument("--out_data", type=str, dest="out_data", metavar="out.csv", required=True, help="CSV file containing the test data, used to compute the model statistics\n")
    parser.add_argument("--model_snp", type=str, dest="model_snp", metavar="model_snp.pkl", required=True, help="Specifies the name of the trained snp model to import\n")
    parser.add_argument("--model_indel", type=str, dest="model_indel", metavar="model_indel.pkl", required=True, help="Specifies the name of the trained indel model to import\n")
    parser.add_argument("--allele_frequency_list", type=str, dest="allele_frequency_list", metavar="frequency1,frequecy2,frequency3", required=False, help="Comma separated list of allele frequency sources that should be used as basis to get the maximum allele frequency\n")
    parser.add_argument("--feature_list", type=str, dest="feature_list", metavar="feature1,feature2,feature3", required=True, help="Comma separated list of the features used to train the model\n")
    args = parser.parse_args()

    input_data_snp = read_input_data(args.in_data_snp)
    input_data_indel = read_input_data(args.in_data_indel)

    if args.allele_frequency_list:
        allele_frequency_list = args.allele_frequency_list.split(",")
    else:
        allele_frequency_list = []
    feature_list = args.feature_list.split(",")

    # if multiple alleles are reported consider only the first one
    # TODO decide how to handle allele ambiguity
    #input_data["ALT"] = input_data["ALT"].map(lambda x: x.split(",")[0])

    predicted_data = perform_pathogenicity_score_prediction(input_data_snp, input_data_indel, args.model_snp, args.model_indel, allele_frequency_list, feature_list)
    predicted_data.to_csv(args.out_data, index=False, sep="\t", na_rep="NA")
