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
             "PolyPhen": 0.5169017014355943}

median_dict = {"MutationAssessor": 1.87,
               "CONDEL": 0.4805749233199981,
               "EIGEN_PHRED": 3.010301,
               "CADD_PHRED": 3.99,
               "FATHMM_XF": 0.209614,
               "SIFT": 0.153,
               "REVEL": 0.193,
               "PolyPhen": 0.547}

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
#def prepare_input_data(input_data, allele_frequency_list, feature_list):
def prepare_input_data(input_data):
    # fill SegDup missing values with -> 0
    # fill ABB_SCORE missing values with -> 0
    # fill Allele Frequence missing values with -> 0
    # fill missing values from other features with -> median or mean

    print(input_data[["CONDEL", "FATHMM_XF", "EIGEN_PHRED", "MutationAssessor"]])

    for allele_frequency in allele_frequencies:
        input_data[allele_frequency] = input_data[allele_frequency].fillna(0)
        input_data[allele_frequency] = input_data.apply(lambda row: pd.Series(max([float(frequency) for frequency in str(row[allele_frequency]).split("&")], default=np.nan)), axis=1)

    # compute maximum Minor Allele Frequency (MAF)
    if (not "MAX_AF" in input_data.columns) & (not "MaxAF" in input_data.columns):
        input_data["MaxAF"] = input_data.apply(lambda row: pd.Series(max([float(frequency) for frequency in row[allele_frequency_list].tolist()])), axis=1)

    for feature in features:
        if feature == "MaxAF" or feature == "MAX_AF":
            input_data[feature] = input_data[feature].fillna(0)
        elif feature == "segmentDuplication":
            input_data[feature] = input_data.apply(lambda row: max([float(value) for value in str(row[feature]).split("&") if ((value != ".") & (value != "nan") & (value != ""))], default=np.nan), axis=1)
            input_data[feature] = input_data[feature].fillna(0)
        elif feature == "ABB_SCORE":
            input_data[feature] = input_data[feature].fillna(0)
        elif "SIFT" in feature:
            input_data[feature] = input_data.apply(lambda row: min([float(value) for value in str(row[feature]).split("&") if ((value != ".") & (value != "nan") & (value != ""))], default=np.nan), axis=1)
            input_data[feature] = input_data[feature].fillna(median_dict["SIFT"])
        else:
            #print(sum(input_data[feature].isna()))
            input_data[feature] = input_data.apply(lambda row: max([float(value) for value in str(row[feature]).split("&") if ((value != ".") & (value != "nan") & (value != ""))], default=np.nan), axis=1)
            #print(sum(input_data[feature].isna()))
            if ("phastCons" in feature) | ("phyloP" in feature):
                input_data[feature] = input_data[feature].fillna(mean_dict[feature])
            else:
                input_data[feature] = input_data[feature].fillna(median_dict[feature])

    #input_features = np.asarray(input_data[feature_list])

    # TODO add workaround to handle the rare case that for one of the features only NaNs are present and therefor the median leads also to a nan
    # in that case the input features contains NaNs and lead to an error

    #return input_data, input_features

    #print(input_data[["CONDEL", "FATHMM_XF", "EIGEN_PHRED", "MutationAssessor"]])
    #print(input_data[["CONDEL", "FATHMM_XF", "EIGEN_PHRED", "MutationAssessor"]].isna().sum())

    return input_data


def predict_pathogenicity(rf_model_snp, rf_model_indel, input_data_snp, input_features_snp, input_data_indel, input_features_indel):
    #class_prediction_snp = rf_model_snp.predict(input_features_snp)
    score_prediction_snp = pd.DataFrame(rf_model_snp.predict_proba(input_features_snp), columns=["Probability_Benign", "Probability_Pathogenic"])

    #class_prediction_indel = rf_model_indel.predict(input_features_indel)
    score_prediction_indel = pd.DataFrame(rf_model_indel.predict_proba(input_features_indel), columns=["Probability_Benign", "Probability_Pathogenic"])

    input_data_snp["AIDIVA_SCORE"] = score_prediction_snp["Probability_Pathogenic"]
    input_data_indel["AIDIVA_SCORE"] = score_prediction_indel["Probability_Pathogenic"]

    return input_data_snp, input_data_indel


def parallel_pathogenicity_prediction(rf_model, input_data, input_features, n_cores):
    splitted_input_data = np.array_split(input_data, n_cores)

    worker_pool = mp.Pool(n_cores)
    predicted_data = pd.concat(worker_pool.apply(rf_model.predict_proba(), (input_features)))

    input_data["AIDIVA_SCORE"] = predicted_data["Probability_Pathogenic"]

    worker_pool.close()
    worker_pool.join()

    return input_data


def check_coding(coding_region, variant_to_check):
    coding = 0

    #if "chr" in str(variant_to_check["CHROM"]):
    #    chrom_id = str(variant_to_check["CHROM"]).replace("chr", "")
    #else:
    #    chrom_id = str(variant_to_check["CHROM"])

    #coding_region = coding_region[coding_region[0] == chrom_id]

    ## TODO: Exchange with consequence check
    #if not coding_region[(coding_region[1]).le(variant_to_check["POS"]) & (coding_region[2].ge(variant_to_check["POS"]))].empty:
    #if not coding_regions[(coding_regions[0] == chrom_id) & (coding_regions[1].le(variant_to_check["POS"])) & (coding_regions[2].ge(variant_to_check["POS"]))].empty:
    #if coding_regions[(coding_regions.CHROM.str.match(str(variant_to_check["CHROM"]))) & (coding_regions.START.le(variant_to_check["POS"])) & (coding_regions.END.ge(variant_to_check["POS"]))]:
    if any(term for term in coding_variants if term in variant_to_check["Consequence"]):
        coding = 1

    return coding


def parallelized_coding_check(input_data):
    global coding_region
    input_data["CODING"] = input_data.apply(lambda row: check_coding(coding_region, row), axis=1)

    return input_data


def parallelize_dataframe_processing(dataframe, function):
    global num_partitions
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


def perform_pathogenicity_score_prediction(input_data_snp, input_data_indel, rf_model_snp, rf_model_indel, allele_frequency_list, feature_list, coding_regions, n_cores=1):
    t = time.time()

    global coding_region
    coding_region = coding_regions

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
    #with mp.Pool(n_cores) as pool:

    t = time.time()
    input_data_snp = parallelize_dataframe_processing(input_data_snp, parallelized_coding_check)
    input_data_indel = parallelize_dataframe_processing(input_data_indel, parallelized_coding_check)

    print("Finished parallelized coding check in: %.2f seconds" % (time.time() - t))

    #input_data_snp["CODING"] = input_data_snp.apply(lambda row: check_coding(coding_regions, row), axis=1)
        #input_data_snp["CODING"] = pool.apply(check_coding, (coding_regions, input_data_snp))
    #input_data_indel["CODING"] = input_data_indel.apply(lambda row: check_coding(coding_regions, row), axis=1)
        #input_data_indel["CODING"] = pool.apply(check_coding, (coding_regions, input_data_indel))

    #print("Finished coding check in: %.2f seconds" % (time.time() - t))

    t = time.time()
    input_data_snp_coding = input_data_snp[input_data_snp["CODING"] == 1]
    input_data_snp_noncoding = input_data_snp[input_data_snp["CODING"] == 0]
    input_data_indel_coding = input_data_indel[input_data_indel["CODING"] == 1]
    input_data_indel_noncoding = input_data_indel[input_data_indel["CODING"] == 0]

    prepared_input_data_snp_coding = parallelize_dataframe_processing(input_data_snp_coding, prepare_input_data)
    prepared_input_data_indel_coding = parallelize_dataframe_processing(input_data_indel_coding, prepare_input_data)

    if not input_data_snp[input_data_snp["CODING"] == 0].empty:
        prepared_input_data_snp_noncoding = parallelize_dataframe_processing(input_data_snp_noncoding, prepare_input_data)
    else:
        prepared_input_data_snp_noncoding = pd.DataFrame()
    if not input_data_indel[input_data_indel["CODING"] == 0].empty:
        prepared_input_data_indel_noncoding = parallelize_dataframe_processing(input_data_indel_noncoding, prepare_input_data)
    else:
        prepared_input_data_indel_noncoding = pd.DataFrame()

    input_features_snp_coding = np.asarray(prepared_input_data_snp_coding[features])
    input_features_indel_coding = np.asarray(prepared_input_data_indel_coding[features])
    #input_features_snp_noncoding = np.asarray(prepared_input_data_snp_noncoding[features])
    #input_features_indel_noncoding = np.asarray(prepared_input_data_indel_noncoding[features])

    #prepared_input_data_snp_coding, input_features_snp_coding = prepare_input_data(input_data_snp_coding, allele_frequency_list, feature_list)
    #prepared_input_data_indel_coding, input_features_indel_coding = prepare_input_data(input_data_indel_coding, allele_frequency_list, feature_list)
    #prepared_input_data_snp_noncoding, input_features_snp_noncoding = prepare_input_data(input_data_snp_noncoding, allele_frequency_list, feature_list)
    #prepared_input_data_indel_noncoding, input_features_indel_noncoding = prepare_input_data(input_data_indel_noncoding, allele_frequency_list, feature_list)

    print("Finished data preparation in: %.2f seconds" % (time.time() - t))

    ## TODO: Add parallelization
    t = time.time()
    predicted_data_snp_coding, predicted_data_indel_coding = predict_pathogenicity(rf_model_snp, rf_model_indel, prepared_input_data_snp_coding, input_features_snp_coding, prepared_input_data_indel_coding, input_features_indel_coding)
    #predicted_data_snp_coding = parallel_pathogenicity_prediction(rf_model_snp, prepared_input_data_snp_coding, input_features_snp_coding, n_cores)
    #predicted_data_indel_coding = parallel_pathogenicity_prediction(rf_model_indel, prepared_input_data_indel_coding, input_features_indel_coding, n_cores)

    print("Finished prediction in: %.2f seconds" % (time.time() - t))

    # set the score for frameshift variants always to 1.0
    # the following line might produce an SettingWithCopyWarning this Warning should be a false positive in this case
    t = time.time()
    predicted_data_indel_coding.loc[(abs(predicted_data_indel_coding.REF.str.len() - predicted_data_indel_coding.ALT.str.len()) % 3 != 0), "AIDIVA_SCORE"] = 1.0

    # set score for non-coding variants to -1 or NaN
    # the models are only for coding variants
    if not prepared_input_data_snp_noncoding.empty:
        prepared_input_data_snp_noncoding.insert(len(prepared_input_data_snp_noncoding.columns), "AIDIVA_SCORE", np.nan, allow_duplicates=False)
    if not prepared_input_data_indel_noncoding.empty:
        prepared_input_data_indel_noncoding.insert(len(prepared_input_data_indel_noncoding.columns), "AIDIVA_SCORE", np.nan, allow_duplicates=False)
    predicted_data_snp_coding.loc[((predicted_data_snp_coding.CHROM == "Y") | (predicted_data_snp_coding.CHROM == "chrY") | (predicted_data_snp_coding.CHROM == "MT") | (predicted_data_snp_coding.CHROM == "chrM")), "AIDIVA_SCORE"] = np.nan
    predicted_data_indel_coding.loc[((predicted_data_indel_coding.CHROM == "Y") | (predicted_data_indel_coding.CHROM == "chrY") | (predicted_data_indel_coding.CHROM == "MT") | (predicted_data_snp_coding.CHROM == "chrM")), "AIDIVA_SCORE"] = np.nan

    # set splicing donor/acceptor variants to 1.0
    predicted_data_snp_coding.loc[(predicted_data_snp_coding.Consequence.str.contains("splice_acceptor_variant") | predicted_data_snp_coding.Consequence.str.contains("splice_donor_variant")), "AIDIVA_SCORE"] = 1.0
    predicted_data_indel_coding.loc[(predicted_data_indel_coding.Consequence.str.contains("splice_acceptor_variant") | predicted_data_indel_coding.Consequence.str.contains("splice_donor_variant")), "AIDIVA_SCORE"] = 1.0

    # set synonymous variants to 0.0
    predicted_data_snp_coding.loc[(predicted_data_snp_coding.Consequence.str.contains("synonymous")), "AIDIVA_SCORE"] = 0.0
    predicted_data_indel_coding.loc[(predicted_data_indel_coding.Consequence.str.contains("synonymous")), "AIDIVA_SCORE"] = 0.0

    print("Finished score adjustments in: %.2f seconds" % (time.time() - t))

    # combine snp and indel data
    ## TODO: the sort attribute is not present in pandas v0.19 so if this version should be
    #predicted_data_complete = pd.concat([predicted_data_snp, predicted_data_indel], sort=False)
    predicted_data_complete = pd.concat([predicted_data_snp_coding, prepared_input_data_snp_noncoding, predicted_data_indel_coding, prepared_input_data_indel_noncoding])
    predicted_data_complete.sort_values(["CHROM", "POS"], ascending=[True, True], inplace=True)
    predicted_data_complete.reset_index(inplace=True, drop=True)

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
    parser.add_argument("--coding_region", type=str, dest="coding_region", metavar="coding_regions.bed", required=True, help="Bed file containing the coding region of the reference assembly\n")
    args = parser.parse_args()

    input_data_snp = read_input_data(args.in_data_snp)
    input_data_indel = read_input_data(args.in_data_indel)

    if args.allele_frequency_list:
        allele_frequency_list = args.allele_frequency_list.split(",")
    else:
        allele_frequency_list = []
    feature_list = args.feature_list.split(",")

    coding_region = pd.read_csv(args.coding_region, sep="\t", names=["CHROM", "START", "END"], low_memory=False)

    # if multiple alleles are reported consider only the first one
    # TODO decide how to handle allele ambiguity
    #input_data["ALT"] = input_data["ALT"].map(lambda x: x.split(",")[0])

    #input_data_snp = input_data[(input_data["REF"].apply(len) == 1) & (input_data["ALT"].apply(len) == 1)]
    #input_data_indel = input_data[(input_data["REF"].apply(len) > 1) | (input_data["ALT"].apply(len) > 1)]

    predicted_data = perform_pathogenicity_score_prediction(input_data_snp, input_data_indel, args.model_snp, args.model_indel, allele_frequency_list, feature_list, coding_region)

    predicted_data.to_csv(args.out_data, index=False, sep="\t", na_rep="NA")
