import pandas as pd

def get_abb_score(row, abb):
    if row.Chr == 'X':
        if (23, int(row.Pos)) in abb.index:
            return abb.loc[(23, int(row.Pos))].ABB
    elif row.Chr == 'Y':
        if (24, int(row.Pos)) in abb.index:
            return abb.loc[(24, int(row.Pos))].ABB
    else:
        if (int(row.Chr), int(row.Pos)) in abb.index:
            return abb.loc[(int(row.Chr), int(row.Pos))].ABB


def extract_data(data, abb):
    data["ABB_SCORE"] = data.apply(lambda row: get_abb_score(row, abb), axis=1)
    
    return data


def group_and_process_data(abb_data, data): 
    abb = pd.read_csv(abb_data, sep="\t", low_memory=False)
    abb.set_index(["CHROM", "POS"], inplace=True)
    abb.sort_index(axis=0, inplace=True, sort_remaining=True)
    abb.sort_index(axis=1, inplace=True, sort_remaining=True)
    
    data_grouped = [group for key, group in data.groupby("Chr")]
    
    for group in data_grouped:
        group = extract_data(group, abb)
    
    data_combined = pd.concat(data_grouped)
    
    return data_combined


if __name__=='__main__':
    
    #!/usr/bin/python3
    import pandas as pd
    import sys
    import argparse
    
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--in_data', type=str, dest='in_data', metavar='in.csv', required=True, help='CSV file containing the data, you want to extend with the ABB score information\n')
    parser.add_argument('--abb_data', type=str, dest='abb_data', metavar='abb.csv', required=True, help='CSV file containing the ABB score information\n')
    parser.add_argument('--out_data', type=str, dest='out_data', metavar='out.csv', required=True, help='Specifies the extended output file\n')
    args = parser.parse_args()

    input_data = pd.read_csv(args.in_data, sep="\t", low_memory=False)
    processed_data = group_and_process_data(args.abb_data, input_data)
    processed_data.to_csv(args.out_data, sep="\t", index=False)
