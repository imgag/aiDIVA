from sklearn.ensemble import RandomForestClassifier
from sklearn.datasets import make_classification
from sklearn.model_selection import RandomizedSearchCV
from sklearn.model_selection import GridSearchCV
import pandas as pd
import numpy as np
import argparse
import pickle


mean_dict = {"phastCons46mammal": 0.09691308336428194,
             "phastCons46primate": 0.12353343703613741,
             "phastCons46vertebrate": 0.1366339183101041,
             "phyloP46mammal": -0.0063575303590607925,
             "phyloP46primate": -0.012076641890840553,
             "phyloP46vertebrate": 0.06761867323083483,
             "phastCons100": 0.11273633387190414,
             "phyloP100": 0.052907788505469275,
             "custom_MutationAssessor": 1.7961304794577417,
             "fannsdb_CONDEL": 0.49699016949707825,
             "custom_EIGEN_PHRED": 4.342947928406315,
             "CADD_PHRED": 4.471745325,
             "custom_FATHMM_XF": 0.35846023623584666,
             "SIFT": 0.35216996259535444,
             "REVEL": 0.28019263637740743,
             "PolyPhen": 0.5169017014355943}


median_dict = {"custom_MutationAssessor": 1.87,
               "fannsdb_CONDEL": 0.4805749233199981,
               "custom_EIGEN_PHRED": 3.010301,
               "CADD_PHRED": 3.99,
               "custom_FATHMM_XF": 0.209614,
               "SIFT": 0.153,
               "REVEL": 0.193,
               "PolyPhen": 0.547}


random_seed = 14038


def import_model(model_file):
    model_to_import = open(model_file, "rb")
    rf_model = pickle.load(model_to_import)
    model_to_import.close()

    return rf_model


def read_input_data(input_file):
    input_data = pd.read_csv(input_file, sep="\t", low_memory=False)

    return input_data


# make sure to use the same features, that were used for the training of the model
def prepare_input_data(input_data, allele_frequency_list, feature_list):
    # fill SegDup missing values with -> 0
    # fill ABB_SCORE missing values with -> 0
    # fill Allele Frequence missing values with -> 0
    # fill missing values from other features with -> median or mean

    for allele_frequency in allele_frequency_list:
        input_data[allele_frequency] = input_data[allele_frequency].fillna(0)
        input_data[allele_frequency] = input_data[allele_frequency].apply(lambda row: pd.Series(max([float(frequency) for frequency in str(row).split("&")], default=np.nan)))

    # compute maximum Minor Allele Frequency (MAF)
    if (not "MAX_AF" in input_data.columns) & (not "MaxAF" in input_data.columns):
        input_data["MaxAF"] = input_data.apply(lambda row: pd.Series(max([float(frequency) for frequency in row[allele_frequency_list].tolist()])), axis=1)

    for feature in feature_list:
        if feature == "MaxAF" or feature == "MAX_AF":
            input_data[feature] = input_data[feature].fillna(0)
        elif feature == "segmentDuplication":
            input_data[feature] = input_data[feature].apply(lambda row: max([float(value) for value in str(row).split("&") if ((value != ".") & (value != "nan") & (value != ""))], default=np.nan))
            input_data[feature] = input_data[feature].fillna(0)
        elif feature == "ABB_SCORE":
            input_data[feature] = input_data[feature].fillna(0)
        elif "SIFT" in feature:
            input_data[feature] = input_data[feature].apply(lambda row: min([float(value) for value in str(row).split("&") if ((value != ".") & (value != "nan") & (value != ""))], default=np.nan))
            input_data[feature] = input_data[feature].fillna(median_dict["SIFT"])
        else:
            input_data[feature] = input_data[feature].apply(lambda row: max([float(value) for value in str(row).split("&") if ((value != ".") & (value != "nan") & (value != ""))], default=np.nan))
            if ("phastCons" in feature) | ("phyloP" in feature):
                input_data[feature] = input_data[feature].fillna(mean_dict[feature])
            else:
                input_data[feature] = input_data[feature].fillna(median_dict[feature])

    input_features = np.asarray(input_data[feature_list])

    # TODO add workaround to handle the rare case that for one of the features only NaNs are present and therefor the median leads also to a nan
    # in that case the input features contains NaNs and lead to an error

    return input_data, input_features


def predict_pathogenicity(rf_model_snps, rf_model_indel, input_data_snps, input_features_snps, input_data_indel, input_features_indel):
    #class_prediction_snps = rf_model_snps.predict(input_features_snps)
    score_prediction_snps = pd.DataFrame(rf_model_snps.predict_proba(input_features_snps), columns=["Probability_Benign", "Probability_Pathogenic"])

    #class_prediction_indel = rf_model_indel.predict(input_features_indel)
    score_prediction_indel = pd.DataFrame(rf_model_indel.predict_proba(input_features_indel), columns=["Probability_Benign", "Probability_Pathogenic"])

    input_data_snps["AIDIVA_SCORE"] = score_prediction_snps["Probability_Pathogenic"]
    input_data_indel["AIDIVA_SCORE"] = score_prediction_indel["Probability_Pathogenic"]

    return input_data_snps, input_data_indel


def check_coding(coding_regions, variant_to_check):
    coding = 0

    if "chr" in str(variant_to_check["CHROM"]):
        chrom_id = str(variant_to_check["CHROM"]).replace("chr", "")
    else:
        chrom_id = str(variant_to_check["CHROM"])

    if not coding_regions[(coding_regions[0] == chrom_id) & (coding_regions[1].le(variant_to_check["POS"])) & (coding_regions[2].ge(variant_to_check["POS"]))].empty:
    #if coding_regions[(coding_regions.CHROM.str.match(str(variant_to_check["CHROM"]))) & (coding_regions.START.le(variant_to_check["POS"])) & (coding_regions.END.ge(variant_to_check["POS"]))]:
        coding = 1

    return coding


def perform_pathogenicity_score_prediction(input_data_snps, input_data_indel, rf_model_snps, rf_model_indel, allele_frequency_list, feature_list, coding_regions):
    prepared_input_data_snps, input_features_snps = prepare_input_data(input_data_snps, allele_frequency_list, feature_list)
    prepared_input_data_indel, input_features_indel = prepare_input_data(input_data_indel, allele_frequency_list, feature_list)

    rf_model_snps = import_model(rf_model_snps)
    rf_model_indel = import_model(rf_model_indel)

    # filter variants, only coding variants should be scored
    input_data_snps["CODING"] = input_data_snps.apply(lambda row: check_coding(coding_regions, row), axis=1)
    input_data_indel["CODING"] = input_data_indel.apply(lambda row: check_coding(coding_regions, row), axis=1)

    predicted_data_snps, predicted_data_indel = predict_pathogenicity(rf_model_snps, rf_model_indel, prepared_input_data_snps, input_features_snps, prepared_input_data_indel, input_features_indel)

    # set the score for frameshift variants always to 1.0
    # the following line might produce an SettingWithCopyWarning this Warning should be a false positive in this case
    predicted_data_indel.loc[(abs(predicted_data_indel.REF.str.len() - predicted_data_indel.ALT.str.len()) % 3 != 0), "AIDIVA_SCORE"] = 1.0

    # set score for non-coding variants to -1 or NaN
    # the models are only for coding variants
    predicted_data_snps.loc[(predicted_data_snps.CODING == 0), "AIDIVA_SCORE"] = np.NaN
    predicted_data_indel.loc[(predicted_data_indel.CODING == 0), "AIDIVA_SCORE"] = np.NaN
    predicted_data_snps.loc[((predicted_data_snps.CHROM == "Y") | (predicted_data_snps.CHROM == "chrY") | (predicted_data_snps.CHROM == "MT") | (predicted_data_snps.CHROM == "chrM")), "AIDIVA_SCORE"] = np.NaN
    predicted_data_indel.loc[((predicted_data_indel.CHROM == "Y") | (predicted_data_indel.CHROM == "chrY") | (predicted_data_indel.CHROM == "MT") | (predicted_data_snps.CHROM == "chrM")), "AIDIVA_SCORE"] = np.NaN

    # set splicing donor/acceptor variants to 1.0
    predicted_data_snps.loc[(predicted_data_snps.Consequence.str.contains("splice_acceptor_variant") | predicted_data_snps.Consequence.str.contains("splice_donor_variant")), "AIDIVA_SCORE"] = 1.0
    predicted_data_indel.loc[(predicted_data_indel.Consequence.str.contains("splice_acceptor_variant") | predicted_data_indel.Consequence.str.contains("splice_donor_variant")), "AIDIVA_SCORE"] = 1.0

    # set synonymous variants to 0.0
    predicted_data_snps.loc[(predicted_data_snps.Consequence.str.contains("synonymous")), "AIDIVA_SCORE"] = 0.0
    predicted_data_indel.loc[(predicted_data_indel.Consequence.str.contains("synonymous")), "AIDIVA_SCORE"] = 0.0

    # combine snps and indel data
    ## TODO: the sort attribute is not present in pandas v0.19 so if this version should be
    #predicted_data_complete = pd.concat([predicted_data_snps, predicted_data_indel], sort=False)
    predicted_data_complete = pd.concat([predicted_data_snps, predicted_data_indel])
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

    #input_data_snps = input_data[(input_data["REF"].apply(len) == 1) & (input_data["ALT"].apply(len) == 1)]
    #input_data_indel = input_data[(input_data["REF"].apply(len) > 1) | (input_data["ALT"].apply(len) > 1)]

    predicted_data = perform_pathogenicity_score_prediction(input_data_snp, input_data_indel, args.model_snp, args.model_indel, allele_frequency_list, feature_list, coding_region)

    predicted_data.to_csv(args.out_data, index=False, sep="\t", na_rep="NA")
