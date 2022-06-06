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