import pandas as pd
import numpy as np
import sys


data = pd.read_csv(sys.argv[1], sep="\t", low_memory=False)

#data.rename(columns={"#CHROM": "Chr"}, inplace=True)
#data.rename(columns={"POS": "Pos"}, inplace=True)

repeats = pd.read_csv(sys.argv[2], sep="\t", low_memory=False, header=None)
repeats.set_index([1], inplace=True)

repeats_chr1 = repeats.loc["chr1"]
repeats_chr2 = repeats.loc["chr2"]
repeats_chr3 = repeats.loc["chr3"]
repeats_chr4 = repeats.loc["chr4"]
repeats_chr5 = repeats.loc["chr5"]
repeats_chr6 = repeats.loc["chr6"]
repeats_chr7 = repeats.loc["chr7"]
repeats_chr8 = repeats.loc["chr8"]
repeats_chr9 = repeats.loc["chr9"]
repeats_chr10 = repeats.loc["chr10"]
repeats_chr11 = repeats.loc["chr11"]
repeats_chr12 = repeats.loc["chr12"]
repeats_chr13 = repeats.loc["chr13"]
repeats_chr14 = repeats.loc["chr14"]
repeats_chr15 = repeats.loc["chr15"]
repeats_chr16 = repeats.loc["chr16"]
repeats_chr17 = repeats.loc["chr17"]
repeats_chr18 = repeats.loc["chr18"]
repeats_chr19 = repeats.loc["chr19"]
repeats_chr20 = repeats.loc["chr20"]
repeats_chr21 = repeats.loc["chr21"]
repeats_chr22 = repeats.loc["chr22"]
repeats_chrX = repeats.loc["chrX"]
repeats_chrY = repeats.loc["chrY"]

data_grouped = [group for key, group in data.groupby("Chr")] 


def extract_simple_repeats(row):
    if row["Chr"] == "1":
        current_repeat_table = repeats_chr1
    elif row["Chr"] == "2":
        current_repeat_table = repeats_chr2
    elif row["Chr"] == "3":
        current_repeat_table = repeats_chr3
    elif row["Chr"] == "4":
        current_repeat_table = repeats_chr4
    elif row["Chr"] == "5":
        current_repeat_table = repeats_chr5
    elif row["Chr"] == "6":
        current_repeat_table = repeats_chr6
    elif row["Chr"] == "7":
        current_repeat_table = repeats_chr7
    elif row["Chr"] == "8":
        current_repeat_table = repeats_chr8
    elif row["Chr"] == "9":
        current_repeat_table = repeats_chr9
    elif row["Chr"] == "10":
        current_repeat_table = repeats_chr10
    elif row["Chr"] == "11":
        current_repeat_table = repeats_chr11
    elif row["Chr"] == "12":
        current_repeat_table = repeats_chr12
    elif row["Chr"] == "13":
        current_repeat_table = repeats_chr13
    elif row["Chr"] == "14":
        current_repeat_table = repeats_chr14
    elif row["Chr"] == "15":
        current_repeat_table = repeats_chr15
    elif row["Chr"] == "16":
        current_repeat_table = repeats_chr16
    elif row["Chr"] == "17":
        current_repeat_table = repeats_chr17
    elif row["Chr"] == "18":
        current_repeat_table = repeats_chr18
    elif row["Chr"] == "19":
        current_repeat_table = repeats_chr19
    elif row["Chr"] == "20":
        current_repeat_table = repeats_chr20
    elif row["Chr"] == "21":
        current_repeat_table = repeats_chr21
    elif row["Chr"] == "22":
        current_repeat_table = repeats_chr22
    elif row["Chr"] == "X":
        current_repeat_table = repeats_chrX
    elif row["Chr"] == "Y":
        current_repeat_table = repeats_chrY
        
    current_repeat_table_entries = current_repeat_table.loc[(current_repeat_table[2] <= int(row["Pos"])) & (current_repeat_table[3] >= int(row["Pos"]))]
            
    occurences = len(current_repeat_table_entries)
    
    if occurences == 1:
        repeat_length = str(current_repeat_table_entries[5][0])
        repeat_region_start = str(current_repeat_table_entries[2][0])
        repeat_region_end = str(current_repeat_table_entries[3][0])
        
        repeat_region = ">".join(["chr" + str(row["Chr"]), repeat_region_start, repeat_region_end])
            
        return [repeat_region, repeat_length]
        
    elif occurences > 1:
        repeat_length_entries = [str(value) for value in current_repeat_table_entries[5].tolist()]
        repeat_region_start = [str(value) for value in current_repeat_table_entries[2].tolist()]
        repeat_region_end = [str(value) for value in current_repeat_table_entries[3].tolist()]
        
        repeat_region_entries = [">".join(["chr" + str(row["Chr"]), repeat_region_start[i], repeat_region_end[i]]) for i in range(len(repeat_region_start))]
        
        repeat_length = "&".join(repeat_length_entries)
        repeat_region = "&".join(repeat_region_entries)
                
        return [repeat_region, repeat_length]
    

def extract_data(data):
    data[["SimpleTandemRepeatRegion","SimpleTandemRepeatLength"]] = data.apply(lambda x: pd.Series(extract_simple_repeats(x)), axis=1)
    return data

for group in data_grouped:
    group = extract_data(group)

data_combined = pd.concat(data_grouped)
data_combined.to_csv(sys.argv[3], sep="\t", index=False)
