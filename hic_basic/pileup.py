import numpy as np
import pandas as pd
import skimage
from tqdm import tqdm

from .plot.hic import cool2mat

def cool2mat_OE(coolp, chrom, expected, balance=False):
    """
    Fetch cooler matrix and calculate OE.
    Input:
        coolp: path to cooler
        chrom: chromosome
    """
    expected = expected.query('(region1 == @chrom) and (region2 == @chrom)')
    if balance:
        expected = expected[
            ["dist","balanced.avg"]
            ].set_index("dist")["balanced.avg"].to_dict()
    else:
        expected = expected[
            ["dist","count.avg"]
            ].set_index("dist")["count.avg"].to_dict()
    raw_mat = cool2mat(coolp, chrom, balance=balance)
    if isinstance(raw_mat.index, pd.MultiIndex):
        raw_mat.index = raw_mat.index.get_level_values(1)
    if isinstance(raw_mat.columns, pd.MultiIndex):
        raw_mat.columns = raw_mat.columns.get_level_values(1)
    exp_mat = np.zeros_like(raw_mat.values).astype(float)
    for k, v in expected.items():
        np.fill_diagonal(exp_mat[k:], v)
        np.fill_diagonal(exp_mat[:,k:], v)
    return (raw_mat / exp_mat).fillna(0)
def add_diag_law(raw_mat, binsize=20000, power=0.25):
    """
    Add artificial diagonal law to the matrix.
    """
    w_mat = np.zeros_like(raw_mat.values).astype(float)
    for k in range(0, raw_mat.shape[0]):
        weight = (k * binsize + binsize)**(-power)
        np.fill_diagonal(w_mat[k:], weight)
        np.fill_diagonal(w_mat[:,k:], weight)
    return raw_mat * w_mat
def block_pileup(coolp, refs, expected=None, power=0.25, give_snips=False, balance=False, debug=False):
    """
    Fetch regions from cooler, iterate chrom by chrom.
    Input:
        coolp: path to cooler file
        refs: list of (chrom, start, end)
    """
    if all((i in tads.columns) for i in ["chrom1","start1","start2"]):
        format = "bedpe"
        chrom_col, start_col, end_col = "chrom1", "start1", "start2"
    elif all((i in tads.columns) for i in ["chrom","start","end"]):
        format =  "bed"
        chrom_col, start_col, end_col = "chrom", "start", "end"
    else:
        raise ValueError("No chrom1/start1/start2 or chrom/start/end found")
    print(f"ref is treated as {format}")
    chroms = tads[chrom_col].unique()
    #print(chroms)
    all_snips = []
    for chrom, tad_chunk in tqdm(tads.groupby(chrom_col), desc="chrom", total=len(chroms)):
        if expected is not None:
            chrom_mat = cool2mat_OE(str(coolp), chrom, expected, balance=balance)
            chrom_mat = add_diag_law(chrom_mat, power=power)
        else:
            chrom_mat = cool2mat(str(coolp), chrom, balance=balance)
        # NOTE: this is a temporary fix for the bug in cooler
        chrom_mat = chrom_mat.loc[
            ~chrom_mat.index.duplicated(keep="first"),
            ~chrom_mat.columns.duplicated(keep="first")
        ]
        # if 2-layer index, only keep layer 1 (pos)
        if isinstance(chrom_mat.index, pd.MultiIndex):
            chrom_mat = chrom_mat.droplevel(0, axis=0)
        if isinstance(chrom_mat.columns, pd.MultiIndex):
            chrom_mat = chrom_mat.droplevel(0, axis=1)
        chrom_snips = []
        for i, row in tad_chunk.iterrows():
            start, end = row[start_col], row[end_col]
            length = end - start
            left = start-length if start-length > 0 else 0
            # for multiindex dataframe, use layer 1 (pos)
            chrom_max = chrom_mat.index.get_level_values(1).max() if isinstance(chrom_mat.index, pd.MultiIndex) else chrom_mat.index.max()
            if debug:
                print(end, length, chrom_max)
            right = end+length if end+length < chrom_max else chrom_max
            snip = chrom_mat.loc[left:right, left:right]
            chrom_snips.append(snip)
        all_snips.extend(chrom_snips)
    #return chrom_snips
    resize_snips = []
    for snip in all_snips:
        snip = snip.values
        #snip = (snip - np.nanmean(snip)) / np.nanstd(snip)
        #resize_snips.append(cv2.resize(snip, (100,100)))
        resize_snips.append(skimage.transform.resize(
            snip, (90,90), preserve_range=True
            ))
    mat = np.nanmean(np.array(resize_snips), axis=0)
    if give_snips:
        return mat, resize_snips, all_snips
    else:
        return mat
def asymmetric_pileup(coolp, refs, expand, expected=None, binsize=None, power=0.25, give_snips=False, balance=False):
    """
    Fetch any regions from cooler, iterate chrom by chrom.
    TODO: support inter-chrom regions
    Input:
        coolp: path to cooler file
        refs: list of (chrom1, start1, end1, chrom2, start2, end2)
        expand: treat input refs as point, expand by expand bp
        expected: path to expected file, will give OBS/EXP if provided
        binsize: if provided, will use this to give index and columns for output matrix
        power: power to add to the diagonal, use this to visualize near-diagonal features such as TADs
        give_snips: return the snips used to calculate the pileup
        balance: use balanced matrix
    Output:
        mat: pileup matrix
    """
    # Input must be bedpe-like
    if all((i in refs.columns) for i in ["chrom1","start1","start2","chrom2","end1","end2"]):
        pass
    else:
        raise ValueError("Input refs should have columns: chrom1, start1, end1, chrom2, start2, end2")
    chroms = refs["chrom1"].unique()
    all_snips = []
    skipped = 0
    for chrom, ref_chunk in tqdm(refs.groupby("chrom1"), desc="chrom", total=len(chroms)):
        if expected is not None:
            chrom_mat = cool2mat_OE(str(coolp), chrom, expected, balance=balance)
            if power is not None:
                chrom_mat = add_diag_law(chrom_mat, power=power)
        else:
            chrom_mat = cool2mat(str(coolp), chrom, balance=balance)
        # NOTE: this is a temporary fix for the bug in cooler
        chrom_mat = chrom_mat.loc[
            ~chrom_mat.index.duplicated(keep="first"),
            ~chrom_mat.columns.duplicated(keep="first")
        ]
        chrom_snips = []
        for i, row in ref_chunk.iterrows():
            start1, end1, start2, end2 = row[["start1","end1","start2","end2"]]
            left1, left2 = start1 - expand, start2 - expand
            right1, right2 = end1 + expand, end2 + expand
            if any([
                left1 < 0,
                left2 < 0,
                right1 > chrom_mat.index.max(),
                right2 > chrom_mat.columns.max()
            ]):
                skipped += 1
                continue
            snip = chrom_mat.loc[left1:right1, left2:right2]
            chrom_snips.append(snip)
        all_snips.extend(chrom_snips)
    all_snip_values = [i.values for i in all_snips]
    print(f"Skipped {skipped} regions")
    mat = np.nanmean(np.array(all_snip_values), axis=0)
    mat = pd.DataFrame(mat)
    if binsize is not None:
        mat.index = np.arange(-expand, expand+binsize+1, binsize)
        mat.columns = np.arange(-expand, expand+binsize+1, binsize)
    if give_snips:
        return mat, all_snips
    else:
        return mat

# --- strength for different features --- #

def TAD_pileup_strength(mat, strip=3)->float:
    """
    Calculate TAD strength from pileup matrix.
    To remove influence from diagonal, we calculate region around domain loop.
    The size of region is determined by strip.
    Input:
        mat: 2D numpy array
        strip: min distance to diagonal, count by pixels
    Output:
        strength: float
    """
    if isinstance(mat, pd.DataFrame):
        mat = mat.values
    e = strip
    inter_tad_1 = mat[45+e:60,15+e:30]
    inter_tad_2 = mat[60:75-e,30:45-e]
    tad = mat[45+e:60,30:45-e]
    strength = tad.sum() / (0.5 * (inter_tad_1.sum() + inter_tad_2.sum()))
    #print(tad.sum(), inter_tad_1.sum(), inter_tad_2.sum())
    return strength