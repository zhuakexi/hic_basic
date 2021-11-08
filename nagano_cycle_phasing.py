import pandas as pd
import os
import gzip
import numpy as np
def contact_describe(cell_name:str,c1=1,p1=2,c2=3,p2=4) -> pd.Series:
    # get cell's basic statistics, defined in Nagano2017
    contacts = pd.read_table(cell_name, header=None, comment="#",low_memory=False)
    new_columns = list(contacts.columns)
    new_columns[c1] = "chr1"
    new_columns[p1] = "pos1"
    new_columns[c2] = "chr2"
    new_columns[p2] = "pos2"
    contacts.columns = new_columns
    intra = contacts.query(' chr1 == chr2 ')
    distances = abs(intra["pos1"] - intra["pos2"])
    
    all_ = len(distances[23_000 < distances])
    short = len(distances[(23_000 < distances) & (distances < 2_000_000)])
    mitotic = len(distances[(2_000_000 < distances) & (distances < 12_000_000)])
    farAvg = distances[(4_500_000 < distances) & (distances < 225_000_000)]
    
    mitotic_r = mitotic/all_
    short_r = short/all_
    
    # assign to different stages on Peter's cirtera
    if mitotic_r >= 0.3 and short_r <= 0.5:
        group = "Post-M"
    elif short_r > 0.5 and short_r + 1.8*mitotic_r > 1.0:
        group = "Pre-M"
    elif short_r <= 0.63:
        group = "G1"
    elif 0.63 < short_r <= 0.785:
        group = "early/mid-S"
    elif short_r > 0.785:
        group = "mid-S/G2"
    else:
        group = "blank"
    
    return {"short%":short_r, "mitotic%":mitotic_r, "farAvg":farAvg.mean(),"group":group }
def Nagano_ordering(metrics:pd.DataFrame):
    # ordering within group using Nagano2017's metric
    # Input:
    #   metrics: dataframe, must contain 
    #            short% mitotic% farAvg repli_score groups
    #            
    #            groups: group markers, same length with metric list
    #            factors within Pre-M, mid-S/G2, early/mid-S,
    #            Post-M, G1, and blank
    # Output:
    #    dataframe with new float col named ordering_value
    metrics_e = metrics.assign(
        scaled_short = (metrics["short%"] - metrics["short%"].mean())/metrics["short%"].std(),
        scaled_farAvg = (metrics["farAvg"] - metrics["farAvg"].mean())/metrics["farAvg"].std(),
        shallow_short = metrics["short%"]/metrics["short%"].var(),
        shallow_rs = metrics["repli_score"]/metrics["repli_score"].var()
    )
    ordering_value = []
    for index, row in metrics_e.iterrows():
        if row["group"] == "Post-M":
            ordering_value.append(
                -row["mitotic%"]
            )
        elif row["group"] == "G1":
            ordering_value.append(
                row["scaled_short"] + row["scaled_farAvg"]
            )
        elif row["group"] == "early/mid-S":
            ordering_value.append(
                row["shallow_short"] + row["shallow_rs"]
            )
        elif row["group"] == "mid-S/G2":
            ordering_value.append(
                row["shallow_short"] - row["shallow_rs"]
            )
        elif row["group"] == "Pre-M":
            ordering_value.append(
                row["mitotic%"]
            )
        elif row["group"] == "blank":
            ordering_value.append(np.nan)
    return metrics.assign(ordering_value = ordering_value)
def order_sample_within_group(chunk,shift):
    # Input
    # Output: new chunk with ordering col
    new_chunk = chunk.sort_values("ordering_value").assign(
        index_order = range(0,chunk.shape[0]))
    new_chunk["index_order"] = new_chunk["index_order"] + shift
    return new_chunk
def order_sample(ordering_values):
    # order
    # Input:
    #    ordering_values: dataframe, must have
    #      sample_name, group, ordering_value, cycle_num
    # Output:
    #    dataframe with new int col named index_order
    group_order = ["Post-M","G1","early/mid-S","mid-S/G2","Pre-M"]
    group_order = [group for group in group_order if group in ordering_values["group"].unique()]
    shift_map = list(ordering_values["group"].value_counts().loc[group_order].cumsum().values[:-1])
    shift_map = dict(zip(group_order,[0] + shift_map))
    index_order = []
    output = []
    for group, chunk in ordering_values.groupby("group"):
        if group != "blank":
            output.append(
                order_sample_within_group(chunk,shift_map[group])
            )
        else:
            output.append(
                chunk.assign(index_order=np.nan)
            )
    return pd.concat(output)