import pandas as pd
import numpy as np
import sys


data = pd.read_csv(sys.argv[1], sep="\t", low_memory=False)
seg_dup = pd.read_csv(sys.argv[2], sep="\t", low_memory=False, header=None)
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


def extract_seg_dups(row):
    if str(row["Chr"]) == "1":
        current_seg_dup_table = seg_dup_chr1
    elif str(row["Chr"]) == "2":
        current_seg_dup_table = seg_dup_chr2
    elif str(row["Chr"]) == "3":
        current_seg_dup_table = seg_dup_chr3
    elif str(row["Chr"]) == "4":
        current_seg_dup_table = seg_dup_chr4
    elif str(row["Chr"]) == "5":
        current_seg_dup_table = seg_dup_chr5
    elif str(row["Chr"]) == "6":
        current_seg_dup_table = seg_dup_chr6
    elif str(row["Chr"]) == "7":
        current_seg_dup_table = seg_dup_chr7
    elif str(row["Chr"]) == "8":
        current_seg_dup_table = seg_dup_chr8
    elif str(row["Chr"]) == "9":
        current_seg_dup_table = seg_dup_chr9
    elif str(row["Chr"]) == "10":
        current_seg_dup_table = seg_dup_chr10
    elif str(row["Chr"]) == "11":
        current_seg_dup_table = seg_dup_chr11
    elif str(row["Chr"]) == "12":
        current_seg_dup_table = seg_dup_chr12
    elif str(row["Chr"]) == "13":
        current_seg_dup_table = seg_dup_chr13
    elif str(row["Chr"]) == "14":
        current_seg_dup_table = seg_dup_chr14
    elif str(row["Chr"]) == "15":
        current_seg_dup_table = seg_dup_chr15
    elif str(row["Chr"]) == "16":
        current_seg_dup_table = seg_dup_chr16
    elif str(row["Chr"]) == "17":
        current_seg_dup_table = seg_dup_chr17
    elif str(row["Chr"]) == "18":
        current_seg_dup_table = seg_dup_chr18
    elif str(row["Chr"]) == "19":
        current_seg_dup_table = seg_dup_chr19
    elif str(row["Chr"]) == "20":
        current_seg_dup_table = seg_dup_chr20
    elif str(row["Chr"]) == "21":
        current_seg_dup_table = seg_dup_chr21
    elif str(row["Chr"]) == "22":
        current_seg_dup_table = seg_dup_chr22
    elif str(row["Chr"]) == "X":
        current_seg_dup_table = seg_dup_chrX
    elif str(row["Chr"]) == "Y":
        current_seg_dup_table = seg_dup_chrY
    
    
    #occurences = sum((seg_dup.index.get_level_values(1) == "chr" + str(row["Chr"])) & (seg_dup.index.get_level_values(2) <= int(row["Pos"])) & (seg_dup.index.get_level_values(3) >= int(row["Pos"])))
    
    filtered_seg_dup_table = current_seg_dup_table.loc[(current_seg_dup_table.index.get_level_values(2) <= int(row["Pos"])) & (current_seg_dup_table.index.get_level_values(3) >= int(row["Pos"]))]
    
    occurences = len(filtered_seg_dup_table)
    
    if occurences == 1:
        value = filtered_seg_dup_table[26]
        value_rounded = ("%.6f" % round(value, 6))
        
        return value_rounded
    elif occurences > 1:        
        values = filtered_seg_dup_table[26].tolist() 
        values_rounded = [("%.6f" % round(value, 6)) for value in values]
        
        return "&".join(values_rounded)


def extract_data(data):
    data[["SegmentDuplication"]] = data.apply(lambda x: pd.Series(extract_seg_dups(x)), axis=1)
    return data


for group in data_grouped:
    group = extract_data(group)
    
data_combined = pd.concat(data_grouped)
data_combined.to_csv(sys.argv[3], sep="\t", index=False)
