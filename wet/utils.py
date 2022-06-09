import os

def check_input(filesp, cols):
    """
    Check if the input files are valid.
    Input:
        filesp: dataframe storing file paths.
        cols: cols to check; list or string.
    """
    missing = {}
    valid_samples = {}
    valid_files = {}
    missing_nums = {}
    if isinstance(cols, str):
        cols = [cols]
    for col in cols:
        missing[col] = []
        valid_samples[col] = []
        valid_files[col] = []
        missing_nums[col] = 0
        for sample, file in zip(filesp.index, filesp[col]):
            if not os.path.isfile(file):
                missing[col].append(file)
                continue
            valid_samples[col].append(sample)
            valid_files[col].append(file)
            missing_nums[col] += 1
        if len(missing[col]):
            print(missing[col])
            print("Warning: {} missing {} input file .".format(col, missing_nums[col]))
            print(",".join(missing[col]))
    if len(cols) == 1:
        return valid_samples[col], valid_files[col]
def two_sets(ref, check, warning=False):
    """
    Check two sample list.
    Input:
        ref: reference list, iterable.
        check: list to be checked, iterable.
    Output:
        [missing, extra]; list of sets
    """
    missing_sample = set(ref) - set(check)
    extra = set(check) - set(ref)
    if warning:
        if len(missing_sample):
            print("[two_sets] Warning: can't find value for %d samples:" % len(missing_sample))
            print(",".join(list(missing_sample)[:20]))
            if len(extra) > 20:
                print("...")
        if len(extra):
            print("[two_sets] Warning: find extra value for:")
            print(",".join(list(extra)[:20]))
            if len(extra) > 20:
                print("...")
    return missing_sample, extra