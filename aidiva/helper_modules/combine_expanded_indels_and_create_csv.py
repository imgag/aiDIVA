import argparse
import logging
import multiprocessing as mp
import numpy as np
import pandas as pd

from functools import partial


logger = logging.getLogger(__name__)


def annotate_indels_with_combined_snps_information(row, grouped_expanded_vcf, feature):
    if grouped_expanded_vcf[feature].get_group(row["INDEL_ID"]).empty:
        logger.error(f"Could not combine expanded InDels, INDEL_ID {row['INDEL_ID']} missing in data!")
        return np.nan
    else:
        return grouped_expanded_vcf[feature].get_group(row["INDEL_ID"]).median()


def combine_vcf_dataframes(feature_list, grouped_expanded_vcf, vcf_as_dataframe):
    for feature in feature_list:
        if (feature == "MaxAF") or (feature == "MAX_AF"):
            continue
        elif (feature == "simpleRepeat"):
            continue
        elif (feature == "oe_lof"):
            continue
        else:
            vcf_as_dataframe[feature] = vcf_as_dataframe.apply(lambda row : pd.Series(annotate_indels_with_combined_snps_information(row, grouped_expanded_vcf, feature)), axis=1)

    vcf_as_dataframe["ada_score"] = vcf_as_dataframe.apply(lambda row : pd.Series(annotate_indels_with_combined_snps_information(row, grouped_expanded_vcf, "ada_score")), axis=1)
    vcf_as_dataframe["rf_score"] = vcf_as_dataframe.apply(lambda row : pd.Series(annotate_indels_with_combined_snps_information(row, grouped_expanded_vcf, "rf_score")), axis=1)

    return vcf_as_dataframe


def parallelized_indel_combination(vcf_as_dataframe, expanded_vcf_as_dataframe, features, num_cores):
    feature_list = features

    for feature in feature_list:
        if (feature == "MaxAF") or (feature == "MAX_AF"):
            continue
        elif (feature == "simpleRepeat"):
            continue
        elif (feature == "oe_lof"):
            continue
        elif (feature == "SIFT"):
            expanded_vcf_as_dataframe[feature] = expanded_vcf_as_dataframe[feature].apply(lambda row: min([float(value) for value in str(row).split("&") if ((value != ".") & (value != "nan"))], default=np.nan))
        else:
            expanded_vcf_as_dataframe[feature] = expanded_vcf_as_dataframe[feature].apply(lambda row: max([float(value) for value in str(row).split("&") if ((value != ".") & (value != "nan"))], default=np.nan))

    expanded_vcf_as_dataframe["ada_score"] = expanded_vcf_as_dataframe["ada_score"].apply(lambda row: max([float(value) for value in str(row).split("&") if ((value != ".") & (value != "nan"))], default=np.nan))
    expanded_vcf_as_dataframe["rf_score"] = expanded_vcf_as_dataframe["rf_score"].apply(lambda row: max([float(value) for value in str(row).split("&") if ((value != ".") & (value != "nan"))], default=np.nan))

    grouped_expanded_vcf = expanded_vcf_as_dataframe.groupby("INDEL_ID")

    num_partitions = num_cores * 2

    if len(vcf_as_dataframe) <= num_partitions:
        dataframe_splitted = np.array_split(vcf_as_dataframe, 1)
    else:
        dataframe_splitted = np.array_split(vcf_as_dataframe, num_partitions)

    try:
        function_to_parallelize = partial(combine_vcf_dataframes, feature_list, grouped_expanded_vcf)
        pool = mp.Pool(num_cores)
        vcf_as_dataframe = pd.concat(pool.map(function_to_parallelize, dataframe_splitted))
    finally:
        pool.close()
        pool.join()

    return vcf_as_dataframe


def write_vcf_to_csv(vcf_combined_as_dataframe, out_file):
    vcf_combined_as_dataframe.to_csv(out_file, sep="\t", encoding="utf-8", index=False)


if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--in_data", type=str, dest="in_data", metavar="input.csv", required=True, help="InDel CSV file to combine\n")
    parser.add_argument("--in_data_expanded", type=str, dest="in_data_expanded", metavar="input_expanded.csv", required=True, help="Expanded InDel CSV file\n")
    parser.add_argument("--out_data", type=str, dest="out_data", metavar="output.csv", required=True, help="CSV file containing the combined CSV files\n")
    parser.add_argument("--feature_list", type=str, dest="feature_list", metavar="feature1,feature2,feature3", required=True, help="Comma separated list with the names of the previously annotated features\n")
    parser.add_argument("--threads", type=str, dest="threads", metavar="1", required=False, help="Number of threads to use\n")
    args = parser.parse_args()

    if args.threads is not None:
        num_cores = int(args.threads)
    else:
        num_cores = 1

    vcf_as_dataframe = pd.read_csv(args.in_data, sep="\t", low_memory=False)
    expanded_vcf_as_dataframe = pd.read_csv(args.in_data_expanded, sep="\t", low_memory=False)

    feature_list = args.feature_list.split(",")
    vcf_combined_as_dataframe = parallelized_indel_combination(vcf_as_dataframe, expanded_vcf_as_dataframe, feature_list, num_cores)
    write_vcf_to_csv(vcf_combined_as_dataframe, args.out_data)
