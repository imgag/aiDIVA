import pandas as pd
import numpy as np
import argparse


def extract_seg_dups(row, current_seg_dup_table):
    filtered_seg_dup_table = current_seg_dup_table.loc[(current_seg_dup_table[2] <= int(row["POS"])) & (current_seg_dup_table[3] >= int(row["POS"]))]

    occurences = len(filtered_seg_dup_table)

    if occurences == 1:
        value = filtered_seg_dup_table[26]
        value_rounded = ("%.6f" % round(value, 6))

        return value_rounded, value_rounded
    elif occurences > 1:
        values = filtered_seg_dup_table[26].tolist()
        values_rounded = [("%.6f" % round(value, 6)) for value in values]

        return "&".join(values_rounded), max(values_rounded)
    else:
        return "", ""


def extract_data(data, current_seg_dup_table):
    data[["SegmentDuplication", "SegDupMax"]] = data.apply(lambda row: pd.Series(extract_seg_dups(row, current_seg_dup_table)), axis=1)
    data["SegDupMax"] = data["SegDupMax"].replace(r"^\s*$", np.nan, regex=True)

    return data


def group_and_process_data(seg_dup_table, data):
    seg_dup = pd.read_csv(seg_dup_table, sep="\t", low_memory=False, header=None)
    seg_dup.set_index([1], inplace=True)

    data_grouped = [group for key, group in data.groupby("CHROM")]

    for group in data_grouped:
        if str(group["CHROM"].iloc[0]) == "MT":
            continue
        else:
            current_seg_dup_table = seg_dup.loc["chr" + str(group["CHROM"].iloc[0])]

        group = extract_data(group, current_seg_dup_table)

    data_combined = pd.concat(data_grouped)

    return data_combined


if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--in_data", type=str, dest="in_data", metavar="data.csv", required=True, help="CSV file containing the data, you want to extend with the segment duplication information\n")
    parser.add_argument("--duplication_data", type=str, dest="duplication_data", metavar="table.csv", required=True, help="CSV file containing the segment duplication information\n")
    parser.add_argument("--out_data", type=str, dest="out_data", metavar="out.csv", required=True, help="Specifies the extended output file\n")
    args = parser.parse_args()

    input_data = pd.read_csv(args.in_data, sep="\t", low_memory=False)

    processed_data = group_and_process_data(args.duplication_data, input_data)

    processed_data.to_csv(args.out_data, sep="\t", index=False)
