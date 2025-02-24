import numpy as np
import warnings
def find_peak_prominence(arr, max_dist=None):
    # this function is from cooltools/lib/peaks.py
    """Find the local maxima of an array and their prominence.
    The prominence of a peak is defined as the maximal difference between the
    height of the peak and the lowest point in the range until a higher peak.

    Parameters
    ----------
    arr : array_like
    max_dist : int
        If specified, the distance to the adjacent higher peaks is limited
        by `max_dist`.

    Returns
    -------
    loc_max_poss : numpy.array
        The positions of local maxima of a given array.

    proms : numpy.array
        The prominence of the detected maxima.
    """

    arr = np.asarray(arr)
    n = len(arr)
    max_dist = len(arr) if max_dist is None else int(max_dist)

    # Finding all local minima and maxima (i.e. points the are lower/higher than
    # both immediate non-nan neighbors).
    arr_nonans = arr[~np.isnan(arr)]
    idxs_nonans2idx = np.arange(arr.size)[~np.isnan(arr)]

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)

        is_min_left = np.r_[False, arr_nonans[:-1] > arr_nonans[1:]]
        is_min_right = np.r_[arr_nonans[:-1] < arr_nonans[1:], False]
        is_loc_min = is_min_left & is_min_right
        loc_min_poss = np.where(is_loc_min)[0]
        loc_min_poss = idxs_nonans2idx[loc_min_poss]

        is_max_left = np.r_[False, arr_nonans[:-1] < arr_nonans[1:]]
        is_max_right = np.r_[arr_nonans[:-1] > arr_nonans[1:], False]
        is_loc_max = is_max_left & is_max_right
        loc_max_poss = np.where(is_loc_max)[0]
        loc_max_poss = idxs_nonans2idx[loc_max_poss]

    # For each maximum, find the position of a higher peak on the left and
    # on the right. If there are no higher peaks within the `max_dist` range,
    # just use the position `max_dist` away.
    left_maxs = -1 * np.ones(len(loc_max_poss), dtype=np.int64)
    right_maxs = -1 * np.ones(len(loc_max_poss), dtype=np.int64)

    for i, pos in enumerate(loc_max_poss):
        for j in range(pos - 1, -1, -1):
            if (arr[j] > arr[pos]) or (pos - j > max_dist):
                left_maxs[i] = j
                break

        for j in range(pos + 1, n):
            if (arr[j] > arr[pos]) or (j - pos > max_dist):
                right_maxs[i] = j
                break

    # Find the prominence of each peak with respect of the lowest point
    # between the peak and the adjacent higher peaks, on the left and the right
    # separately.
    left_max_proms = np.array(
        [
            (
                arr[pos] - np.nanmin(arr[left_maxs[i] : pos])
                if (left_maxs[i] >= 0)
                else np.nan
            )
            for i, pos in enumerate(loc_max_poss)
        ]
    )

    right_max_proms = np.array(
        [
            (
                arr[pos] - np.nanmin(arr[pos : right_maxs[i]])
                if (right_maxs[i] >= 0)
                else np.nan
            )
            for i, pos in enumerate(loc_max_poss)
        ]
    )

    # In 1D, the topographic definition of the prominence of a peak reduces to
    # the minimum of the left-side and right-side prominence.

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)

        max_proms = np.nanmin(np.vstack([left_max_proms, right_max_proms]), axis=0)

    # The global maximum, by definition, does not have higher peaks around it and
    # thus its prominence is explicitly defined with respect to the lowest local
    # minimum. This issue arises only if max_dist was not specified, otherwise
    # the prominence of the global maximum is already calculated with respect
    # to the lowest point within the `max_dist` range.
    # If no local minima are within the `max_dist` range, just use the
    # lowest point.
    global_max_mask = (left_maxs == -1) & (right_maxs == -1)
    if (global_max_mask).sum() > 0:
        global_max_idx = np.where(global_max_mask)[0][0]
        global_max_pos = loc_max_poss[global_max_idx]
        neighbor_loc_mins = (loc_min_poss >= global_max_pos - max_dist) & (
            loc_min_poss < global_max_pos + max_dist
        )
        if np.any(neighbor_loc_mins):
            max_proms[global_max_idx] = arr[global_max_pos] - np.nanmin(
                arr[loc_min_poss[neighbor_loc_mins]]
            )
        else:
            max_proms[global_max_idx] = arr[global_max_pos] - np.nanmin(
                arr[max(global_max_pos - max_dist, 0) : global_max_pos + max_dist]
            )

    return loc_max_poss, max_proms
def _mat_IS(dm_v:np.ndarray, w:int) -> np.ndarray:
    """
    Calculate the insulation score for each bin in the matrix.
    Input:
        dm_v: distance matrix
        w: window size
    Output:
        IS: insulation score, a 1D array
        counts: a 2D array, each row contains the count of left triangle, right triangle, and diamond
            only count the valid values
    """
    N = dm_v.shape[0]
    IS = []
    # catch = []
    counts = []
    x, y = np.indices(dm_v.shape)
    for v in range(N):
        # v for viewport id
        left_triangle = (x >= v - w) & (x < v) & (x < y) \
            & (y >= v - w) & (y < v)
        right_triangle = (x > v) & (x <= v + w) & (x < y) \
            & (y > v) & (y <= v + w)
        diamond = (x > v) & (x <= v + w) & (y >= v - w) & (y < v)
        left_triangle_count = (left_triangle & (~np.isnan(dm_v))).sum()
        right_triangle_count = (right_triangle & (~np.isnan(dm_v))).sum()
        diamond_count = (diamond & (~np.isnan(dm_v))).sum()
        full_triangle_count = w * (w - 1) / 2
        full_diamond_count = w * w
        # debug info
        # IS.append(left_triangle_count)
        # if v == 50:
        #     catch.append(left_triangle)
        # continue
        counts.append((left_triangle_count, right_triangle_count, diamond_count))
        if (left_triangle_count + right_triangle_count < full_triangle_count * 2) \
            or (diamond_count < full_diamond_count):
            IS.append(np.nan)
        else:
            #print("got")
            diamond_value = dm_v[diamond].sum()
            left_triangle_value = dm_v[left_triangle].sum()
            right_triangle_value = dm_v[right_triangle].sum()
            IS.append(
                diamond_value / (left_triangle_value + right_triangle_value)
            )
    return np.array(IS), np.array(counts)
def IS2blocks(borders):
    # prepare 0-1, 2-3, 4-5, ...
    chunks = []
    for i, chunk in borders.groupby(lambda x : x % 2):
        chunks.append(chunk)
    chunks[1].index = chunks[1].index - 1
    tad_blocks0 = pd.concat(
        [
            chunks[0][["chrom","start"]],
            chunks[1][["chrom","start"]]
            ],
        axis=1
    )
    tad_blocks0.columns = ["chrom1","start1","chrom2","start2"]
    # prepare 1-2, 3-4, 5-6, ...
    chunks = []
    for i, chunk in borders.groupby(lambda x : x % 2):
        chunks.append(chunk)
    chunks[0].index = chunks[0].index - 1
    tad_blocks1 = pd.concat(
        [
            chunks[1][["chrom","start"]],
            chunks[0][["chrom","start"]]
            ],
        axis=1
    )
    tad_blocks1.columns = ["chrom1","start1","chrom2","start2"]
    tad_blocks = pd.concat([tad_blocks0, tad_blocks1], axis=0).sort_index()
    tad_blocks = tad_blocks.dropna(how="any")
    tad_blocks = tad_blocks.query("chrom1 == chrom2").copy()
    tad_blocks["start1"] = tad_blocks["start1"].astype(int)
    tad_blocks["start2"] = tad_blocks["start2"].astype(int)
    return tad_blocks