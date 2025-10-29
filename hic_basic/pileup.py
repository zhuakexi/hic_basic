import numpy as np
import pandas as pd
import skimage
from cooler import Cooler
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
    if all((i in refs.columns) for i in ["chrom1","start1","start2"]):
        format = "bedpe"
        chrom_col, start_col, end_col = "chrom1", "start1", "start2"
    elif all((i in refs.columns) for i in ["chrom","start","end"]):
        format =  "bed"
        chrom_col, start_col, end_col = "chrom", "start", "end"
    else:
        raise ValueError("No chrom1/start1/start2 or chrom/start/end found")
    print(f"ref is treated as {format}")
    chroms = refs[chrom_col].unique()
    #print(chroms)
    all_snips = []
    for chrom, tad_chunk in tqdm(refs.groupby(chrom_col), desc="chrom", total=len(chroms)):
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
def asymmetric_pileup(coolp, refs, expand, expected=None, binsize=None, power=0.25, give_snips=False, balance=False, show_progress=True):
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
    snips_added = 0
    skipped = 0
    coolbinsize = Cooler(coolp).binsize

    for chrom, ref_chunk in tqdm(refs.groupby("chrom1"), desc="chrom", total=len(chroms), disable=not show_progress):
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
        # strip first layer if chrom_mat is multiindex
        if isinstance(chrom_mat.index, pd.MultiIndex):
            chrom_mat.index = chrom_mat.index.get_level_values(1)
        if isinstance(chrom_mat.columns, pd.MultiIndex):
            chrom_mat.columns = chrom_mat.columns.get_level_values(1)
        for i, row in ref_chunk.iterrows():
            start1, end1, start2, end2 = row[["start1","end1","start2","end2"]]
            left1, left2 = start1 - expand, start2 - expand
            right1, right2 = end1 + expand - coolbinsize, end2 + expand - coolbinsize # loc will give the right-most bin
            if any([
                left1 < 0,
                left2 < 0,
                right1 > chrom_mat.index.max(),
                right2 > chrom_mat.columns.max()
            ]):
                skipped += 1
                continue
            snip = chrom_mat.loc[left1:right1, left2:right2]
            if give_snips:
                # store the snip
                all_snips.append(snip.values)
            if snips_added == 0:
                # initialize the matrix with the first snip
                mat_sum = snip.values
            else:
                # accumulate for stream mean
                mat_sum += snip.values
            snips_added += 1
    print(f"Skipped {skipped} regions")
    mat_mean = mat_sum / snips_added
    mat_mean = pd.DataFrame(mat_mean)
    if binsize is not None:
        mat_mean.index = np.arange(-expand, expand+binsize, binsize)
        mat_mean.columns = np.arange(-expand, expand+binsize, binsize)
    if give_snips:
        return mat_mean, all_snips
    else:
        return mat_mean

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
def kernel_strength(kernel, matrix, retmat=False):
    """
    Calculate the signal strength of matrix according to the kernel.

    The function computes the signal strength by averaging the values of the kernel
    where the mask is True, and dividing by the mean of the values of the False mask positions.

    Input:
        kernel (np.ndarray): A 2D numpy array representing the kernel. 0 for signal positions, 1 for background positions.
        Size of kernel must be smaller than the matrix size. Only the region defined by the mask is considered.
        matrix (np.ndarray): A 2D numpy array representing the matrix from which the signal strength is calculated.
    retmat (bool): If True, returns the matrix values corresponding to the kernel mask region.

    Returns:
    float: The signal strength.
    """
    assert kernel.shape[0] <= matrix.shape[0] and kernel.shape[1] <= matrix.shape[1], "Kernel size must be smaller than the matrix size."
    # kernel and matrix size must be odd
    assert kernel.shape[0] % 2 == 1 and kernel.shape[1] % 2 == 1, "Kernel size must be odd."
    signal_mask = kernel == 0
    background_mask = kernel == 1
    # pad the kernel to the size of the matrix
    pad_height, pad_width = (matrix.shape[0] - kernel.shape[0]) // 2, (matrix.shape[1] - kernel.shape[1]) // 2
    padded_signal_mask = np.pad(signal_mask, ((pad_height, pad_height), (pad_width, pad_width)), mode='constant', constant_values=False)
    padded_background_mask = np.pad(background_mask, ((pad_height, pad_height), (pad_width, pad_width)), mode='constant', constant_values=False)
    # calculate the signal strength
    signal_strength = np.nanmean(matrix[padded_signal_mask]) / np.nanmean(matrix[padded_background_mask])
    # print(padded_signal_mask.shape)
    # print(matrix[padded_signal_mask].shape)
    if retmat:
        mat = matrix[padded_signal_mask | padded_background_mask].reshape(kernel.shape[0], kernel.shape[1])
        return signal_strength, mat
    return signal_strength