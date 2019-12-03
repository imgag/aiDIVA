import pandas as pd
import numpy as np
import sys


def get_abb_score(row):
    if row.Chr == 'X':
        if (23, int(row.Pos)) in abb.index:
            return abb.loc[(23, int(row.Pos))].ABB
    elif row.Chr == 'Y':
        if (24, int(row.Pos)) in abb.index:
            return abb.loc[(24, int(row.Pos))].ABB
    else:
        if (int(row.Chr), int(row.Pos)) in abb.index:
            return abb.loc[(int(row.Chr), int(row.Pos))].ABB

data = pd.read_csv(sys.argv[1], sep="\t", low_memory=False)

abb = pd.read_csv(sys.argv[2], sep="\t", low_memory=False)

abb.set_index(["CHROM", "POS"], inplace=True)

data["ABB_SCORE"] = data.apply(lambda row: get_abb_score(row), axis=1)

data.to_csv(sys.argv[3], sep="\t", index=False)
