import pandas as pd
import argparse


def extract_simple_repeats(row, current_repeat_table):        
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
    else:
        return ["", ""]
    

def extract_data(data, current_repeat_table):
    data[["SimpleTandemRepeatRegion","SimpleTandemRepeatLength"]] = data.apply(lambda x: pd.Series(extract_simple_repeats(x, current_repeat_table)), axis=1)
    
    return data


def group_and_process_data(repeat_data, data):
    repeats = pd.read_csv(repeat_data, sep="\t", low_memory=False, header=None)
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

    for group in data_grouped:
        if group["Chr"].iloc[0] == "1":
            current_repeat_table = repeats_chr1
        elif group["Chr"].iloc[0] == "2":
            current_repeat_table = repeats_chr2
        elif group["Chr"].iloc[0] == "3":
            current_repeat_table = repeats_chr3
        elif group["Chr"].iloc[0] == "4":
            current_repeat_table = repeats_chr4
        elif group["Chr"].iloc[0] == "5":
            current_repeat_table = repeats_chr5
        elif group["Chr"].iloc[0] == "6":
            current_repeat_table = repeats_chr6
        elif group["Chr"].iloc[0] == "7":
            current_repeat_table = repeats_chr7
        elif group["Chr"].iloc[0] == "8":
            current_repeat_table = repeats_chr8
        elif group["Chr"].iloc[0] == "9":
            current_repeat_table = repeats_chr9
        elif group["Chr"].iloc[0] == "10":
            current_repeat_table = repeats_chr10
        elif group["Chr"].iloc[0] == "11":
            current_repeat_table = repeats_chr11
        elif group["Chr"].iloc[0] == "12":
            current_repeat_table = repeats_chr12
        elif group["Chr"].iloc[0] == "13":
            current_repeat_table = repeats_chr13
        elif group["Chr"].iloc[0] == "14":
            current_repeat_table = repeats_chr14
        elif group["Chr"].iloc[0] == "15":
            current_repeat_table = repeats_chr15
        elif group["Chr"].iloc[0] == "16":
            current_repeat_table = repeats_chr16
        elif group["Chr"].iloc[0] == "17":
            current_repeat_table = repeats_chr17
        elif group["Chr"].iloc[0] == "18":
            current_repeat_table = repeats_chr18
        elif group["Chr"].iloc[0] == "19":
            current_repeat_table = repeats_chr19
        elif group["Chr"].iloc[0] == "20":
            current_repeat_table = repeats_chr20
        elif group["Chr"].iloc[0] == "21":
            current_repeat_table = repeats_chr21
        elif group["Chr"].iloc[0] == "22":
            current_repeat_table = repeats_chr22
        elif group["Chr"].iloc[0] == "X":
            current_repeat_table = repeats_chrX
        elif group["Chr"].iloc[0] == "Y":
            current_repeat_table = repeats_chrY
        else:
            print("ERROR: The file containing the repeat information seems to be broken!", file=sys.stderr)
        
        group = extract_data(group, current_repeat_table)

    data_combined = pd.concat(data_grouped)
    
    return data_combined


if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--in_data', type=str, dest='in_data', metavar='data.csv', required=True, help='CSV file containing the data, you want to extend with the simple repeat information\n')
    parser.add_argument('--repeat_data', type=str, dest='repeat_data', metavar='table.csv', required=True, help='CSV file containing the simple repeat information\n')
    parser.add_argument('--out_data', type=str, dest='out_data', metavar='out.csv', required=True, help='Specifies the extended output file\n')
    args = parser.parse_args()

    input_data = pd.read_csv(args.in_data, sep="\t", low_memory=False)
    
    processed_data = group_and_process_data(args.repeat_data, input_data)

    processed_data.to_csv(args.out_data, sep="\t", index=False)
