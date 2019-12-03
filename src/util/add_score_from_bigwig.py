import pandas as pd
import pyBigWig as pw
import sys


data = pd.read_csv(sys.argv[1], sep="\t", low_memory=False)

bw_file = pw.open(sys.argv[2])
    
if not bw_file.isBigWig():
    print("The given file is not in BigWig format!!!")

data_grouped = [group for key, group in data.groupby("Chr")] 

feature_list = ["gerp_new", "phastCons100", "phyloP100", "phastCons46_mammal", "phastCons46_primates", "phastCons46_vertebrate", "phyloP46_mammal", "phyloP46_primates", "phyloP46_vertebrate"]
data = data.reindex(columns=[*data.columns.tolist(), *feature_list], fill_value="")


def extract_value(row):
    start = row["Pos"] - 1
    end = row["Pos"]
    #end = start + len(row["Alt"])
    
    value = str(bw_file.values("chr" + str(row["Chr"]), start, end)[0])

    return value

def extract_data(data):
    data[[sys.argv[3]]] = data.apply(lambda x: pd.Series(extract_value(x)), axis=1)
    return data


for group in data_grouped:
    group = extract_data(group)

data_combined = pd.concat(data_grouped)
data_combined.to_csv(sys.argv[4], sep="\t", index=False)
