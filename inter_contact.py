import pandas as pd
import numpy

from hires_utils.hires_io import parse_pairs

def expected_inter(pairs:pd.DataFrame)->pd.DataFrame:
    # calculate expected inter contacts(exp(chra,chrb) = chra_total_frag * chrb_total_frag)
    # input: pairs data
    # output: n*n matrix
    
    ## calculate total contacts for each chromosome
    leg1s, leg2s = pairs[["chr1","pos1"]], pairs[["chr2","pos2"]]
    leg1s.columns = leg2s.columns = "chr pos".split()
    legs = leg1s.append(leg2s)
    leg_count = legs.groupby("chr").size()
    
    expected = pd.DataFrame(leg_count).dot(pd.DataFrame(leg_count).T)
    return expected
def pick_triu(matrix:pd.DataFrame) -> pd.DataFrame:
    # pick upper triangle of matrix
    # input: n*n matrix, with full columns and indices
    # return: multi index 1-D DataFrame
    # ignore diagnal elements
    n = len(matrix)
    index = []
    data = []
    for i in range(0, n):
        matrix_stripe = matrix.iloc[i, i + 1 :]
        index.extend([(matrix_stripe.name, sub_index) for sub_index in matrix_stripe.index])
        data.extend(matrix_stripe.values)
    return pd.Series(index=index, data=data)    
def cell_sig(file_path:str)->pd.Series:
    # calculate cell inter-contact-signature 
    # assume phased dip pairs cell as input, upper triangle not needed.
    # return 1035*1 pd.Series with proper index
    pairs = parse_pairs(file_path) 
    ## generate 46*46 contacts count matrix
    gp = pairs.groupby(["chr1","chr2"])
    pairs_counts = gp.size() # split-count method
    for_matrix = pairs_counts.reset_index()
    for_matrix.columns = "chr1 chr2 counts".split()
    matrix = for_matrix.pivot(index="chr1", columns="chr2",values="counts") # trans chromosome pair with 0 contacts will got NA here
    full_matrix = matrix.fillna(0)
    
    full_sig = full_matrix.T + full_matrix # in case input isn't sorted
    normed_sig = full_sig/ expected_inter(pairs) # expected cannot have zero
    
    res = pick_triu(normed_sig)
    return res