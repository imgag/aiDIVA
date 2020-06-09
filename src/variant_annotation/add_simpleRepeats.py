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
    data[["SimpleTandemRepeatRegion", "SimpleTandemRepeatLength"]] = data.apply(lambda x: pd.Series(extract_simple_repeats(x, current_repeat_table)), axis=1)
    
    return data


def group_and_process_data(repeat_data, data):
    repeats = pd.read_csv(repeat_data, sep="\t", low_memory=False, header=None)
    repeats.set_index([1], inplace=True)
    
    data_grouped = [group for key, group in data.groupby("Chr")]

    for group in data_grouped:
        if str(group["Chr"].iloc[0]) == "MT":
            continue
        else:
            current_repeat_table = repeats.loc["chr" + str(group["Chr"].iloc[0])]
        
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
