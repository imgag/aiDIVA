import pandas as pd
import numpy as np
import argparse


def extract_seg_dups(row, current_seg_dup_table):
    filtered_seg_dup_table = current_seg_dup_table.loc[(current_seg_dup_table.index.get_level_values(2) <= int(row["Pos"])) & (current_seg_dup_table.index.get_level_values(3) >= int(row["Pos"]))]
    
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
    seg_dup.set_index([1, 2, 3], inplace=True)

    seg_dup_chr1 = seg_dup.loc["chr1"]
    seg_dup_chr2 = seg_dup.loc["chr2"]
    seg_dup_chr3 = seg_dup.loc["chr3"]
    seg_dup_chr4 = seg_dup.loc["chr4"]
    seg_dup_chr5 = seg_dup.loc["chr5"]
    seg_dup_chr6 = seg_dup.loc["chr6"]
    seg_dup_chr7 = seg_dup.loc["chr7"]
    seg_dup_chr8 = seg_dup.loc["chr8"]
    seg_dup_chr9 = seg_dup.loc["chr9"]
    seg_dup_chr10 = seg_dup.loc["chr10"]
    seg_dup_chr11 = seg_dup.loc["chr11"]
    seg_dup_chr12 = seg_dup.loc["chr12"]
    seg_dup_chr13 = seg_dup.loc["chr13"]
    seg_dup_chr14 = seg_dup.loc["chr14"]
    seg_dup_chr15 = seg_dup.loc["chr15"]
    seg_dup_chr16 = seg_dup.loc["chr16"]
    seg_dup_chr17 = seg_dup.loc["chr17"]
    seg_dup_chr18 = seg_dup.loc["chr18"]
    seg_dup_chr19 = seg_dup.loc["chr19"]
    seg_dup_chr20 = seg_dup.loc["chr20"]
    seg_dup_chr21 = seg_dup.loc["chr21"]
    seg_dup_chr22 = seg_dup.loc["chr22"]
    seg_dup_chrX = seg_dup.loc["chrX"]
    seg_dup_chrY = seg_dup.loc["chrY"]
    
    data_grouped = [group for key, group in data.groupby("Chr")]
    
    for group in data_grouped:
        if str(group["Chr"].iloc[0]) == "1":
            current_seg_dup_table = seg_dup_chr1
        elif str(group["Chr"].iloc[0]) == "2":
            current_seg_dup_table = seg_dup_chr2
        elif str(group["Chr"].iloc[0]) == "3":
            current_seg_dup_table = seg_dup_chr3
        elif str(group["Chr"].iloc[0]) == "4":
            current_seg_dup_table = seg_dup_chr4
        elif str(group["Chr"].iloc[0]) == "5":
            current_seg_dup_table = seg_dup_chr5
        elif str(group["Chr"].iloc[0]) == "6":
            current_seg_dup_table = seg_dup_chr6
        elif str(group["Chr"].iloc[0]) == "7":
            current_seg_dup_table = seg_dup_chr7
        elif str(group["Chr"].iloc[0]) == "8":
            current_seg_dup_table = seg_dup_chr8
        elif str(group["Chr"].iloc[0]) == "9":
            current_seg_dup_table = seg_dup_chr9
        elif str(group["Chr"].iloc[0]) == "10":
            current_seg_dup_table = seg_dup_chr10
        elif str(group["Chr"].iloc[0]) == "11":
            current_seg_dup_table = seg_dup_chr11
        elif str(group["Chr"].iloc[0]) == "12":
            current_seg_dup_table = seg_dup_chr12
        elif str(group["Chr"].iloc[0]) == "13":
            current_seg_dup_table = seg_dup_chr13
        elif str(group["Chr"].iloc[0]) == "14":
            current_seg_dup_table = seg_dup_chr14
        elif str(group["Chr"].iloc[0]) == "15":
            current_seg_dup_table = seg_dup_chr15
        elif str(group["Chr"].iloc[0]) == "16":
            current_seg_dup_table = seg_dup_chr16
        elif str(group["Chr"].iloc[0]) == "17":
            current_seg_dup_table = seg_dup_chr17
        elif str(group["Chr"].iloc[0]) == "18":
            current_seg_dup_table = seg_dup_chr18
        elif str(group["Chr"].iloc[0]) == "19":
            current_seg_dup_table = seg_dup_chr19
        elif str(group["Chr"].iloc[0]) == "20":
            current_seg_dup_table = seg_dup_chr20
        elif str(group["Chr"].iloc[0]) == "21":
            current_seg_dup_table = seg_dup_chr21
        elif str(group["Chr"].iloc[0]) == "22":
            current_seg_dup_table = seg_dup_chr22
        elif str(group["Chr"].iloc[0]) == "X":
            current_seg_dup_table = seg_dup_chrX
        elif str(group["Chr"].iloc[0]) == "Y":
            current_seg_dup_table = seg_dup_chrY
        else:
            print("ERROR: The file containing the duplication information seems to be broken!", file=sys.stderr)

        group = extract_data(group, current_seg_dup_table)
    
    data_combined = pd.concat(data_grouped)
    
    return data_combined


if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--in_data', type=str, dest='in_data', metavar='data.csv', required=True, help='CSV file containing the data, you want to extend with the segment duplication information\n')
    parser.add_argument('--duplication_data', type=str, dest='duplication_data', metavar='table.csv', required=True, help='CSV file containing the segment duplication information\n')
    parser.add_argument('--out_data', type=str, dest='out_data', metavar='out.csv', required=True, help='Specifies the extended output file\n')
    args = parser.parse_args()

    input_data = pd.read_csv(args.in_data, sep="\t", low_memory=False)
    
    processed_data = group_and_process_data(args.duplication_data, input_data)

    processed_data.to_csv(args.out_data, sep="\t", index=False)
