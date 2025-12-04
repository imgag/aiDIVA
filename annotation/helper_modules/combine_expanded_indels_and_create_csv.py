import argparse
import logging
import multiprocessing as mp
import numpy as np
import pandas as pd
import warnings

from functools import partial


logger = logging.getLogger(__name__)


def annotate_indels_with_combined_snps_information(row, grouped_expanded_vcf, feature, SPLICE_VARIANTS, SYNONYMOUS_VARIANTS):
    with warnings.catch_warnings():
        # we expect to see RuntimeWarnings if the values for a certain variant are missing (the following filter will prevent them from bloating the log file)
        warnings.filterwarnings(action='ignore', message='Mean of empty slice')

        if grouped_expanded_vcf[feature].get_group(row["INDEL_ID"]).empty:
            logger.error(f"Could not combine expanded InDels, INDEL_ID {row['INDEL_ID']} missing in data!")
            return np.nan

        else:
            current_group = grouped_expanded_vcf.get_group(row["INDEL_ID"])

            # we only use non-splicing und non-synonymous variants with impact High or Moderate
            current_group = current_group[(~(current_group["MOST_SEVERE_CONSEQUENCE"].str.contains("|".join(SYNONYMOUS_VARIANTS))) & ~(current_group["MOST_SEVERE_CONSEQUENCE"].str.contains("|".join(SPLICE_VARIANTS)))) & ((current_group["IMPACT"] == "HIGH") | (current_group["IMPACT"] == "MODERATE"))]

            # to prevent underscoring of High impact variants
            current_group.loc[(current_group["IMPACT"] == "HIGH") & (current_group[feature].isna()), feature] = 1.0

            return current_group[feature].mean()


def combine_vcf_dataframes(feature_list, SPLICE_VARIANTS, SYNONYMOUS_VARIANTS, grouped_expanded_vcf, vcf_as_dataframe):
    for feature in feature_list:
        if (feature == "MaxAF") or (feature == "MAX_AF"):
            continue

        if feature == "homAF":
            continue

        elif (feature == "simpleRepeat"):
            continue

        elif (feature == "oe_lof"):
            continue

        elif (feature == "HIGH_IMPACT"):
            continue

        elif (feature == "IS_INDEL"):
            continue

        else:
            vcf_as_dataframe[feature] = vcf_as_dataframe.apply(lambda row : pd.Series(annotate_indels_with_combined_snps_information(row, grouped_expanded_vcf, feature, SPLICE_VARIANTS, SYNONYMOUS_VARIANTS)), axis=1)

    return vcf_as_dataframe


def parallelized_indel_combination(vcf_as_dataframe, expanded_vcf_as_dataframe, features, num_cores, CONSTANT_DICTIONARY):
    # get constants
    SPLICE_VARIANTS = CONSTANT_DICTIONARY["SPLICE_VARIANTS"]
    SYNONYMOUS_VARIANTS = CONSTANT_DICTIONARY["SYNONYMOUS_VARIANTS"]

    feature_list = features

    for feature in feature_list:
        if (feature == "MaxAF") or (feature == "MAX_AF"):
            continue

        if feature == "homAF":
            continue

        elif (feature == "simpleRepeat"):
            continue

        elif (feature == "oe_lof"):
            continue

        elif (feature == "HIGH_IMPACT"):
            continue

        elif (feature == "IS_INDEL"):
            continue

        elif (feature == "SIFT"):
            expanded_vcf_as_dataframe[feature] = expanded_vcf_as_dataframe[feature].apply(lambda row: min([float(value) for value in str(row).split("&") if ((value != ".") & (value != "nan") & (value != "NA"))], default=np.nan))

        else:
            expanded_vcf_as_dataframe[feature] = expanded_vcf_as_dataframe[feature].apply(lambda row: max([float(value) for value in str(row).split("&") if ((value != ".") & (value != "nan")  & (value != "NA"))], default=np.nan))

    grouped_expanded_vcf = expanded_vcf_as_dataframe.groupby("INDEL_ID")

    num_partitions = num_cores * 2

    if len(vcf_as_dataframe) <= num_partitions:
        # do not split dataframe
        dataframe_splitted = [vcf_as_dataframe]

    else:
        # usage of floor division (//) makes sure that we get an absolute number as result
        chunk_size = vcf_as_dataframe.shape[0] // num_partitions
        dataframe_splitted = [vcf_as_dataframe[i:i+chunk_size].copy() for i in range(0, vcf_as_dataframe.shape[0], chunk_size)]

    try:
        function_to_parallelize = partial(combine_vcf_dataframes, feature_list, SPLICE_VARIANTS, SYNONYMOUS_VARIANTS, grouped_expanded_vcf)
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
