import itertools
import json
from pathlib import Path
# --- input --- #
def safe_input(file_pat, with_key=False, **kwargs):
    """
    Expand file_pat to see if the input files exist, if not, skip.
    Input:
        file_pat: file pattern
        with_key: return outkeys or not, bool
        kwargs: key word arguments, to expand file_pat with wildcards
    Output:
        list of input files or (list of input files, list of outkeys)
        filtered by the existence of the input files
    """
    # transform all values to str
    file_pat = str(file_pat)
    kwargs = {
        key : [str(value) for value in values]
        for key, values in kwargs.items()
    }
    # get all combinations of values
    key_combines = [dict(zip(kwargs.keys(), combination)) for combination in itertools.product(*kwargs.values())]
    # check file existence and filter
    inputs, outkeys = [], []
    for key_combine in key_combines:
        file = file_pat.format(**key_combine)
        if Path(file).exists():
            inputs.append(file)
            outkeys.append(key_combine)
        else:
            print(f"Waring: {file} not exists, skip {key_combine}")
            pass
    if with_key:
        return inputs, outkeys
    else:
        return inputs
def check_input(file_pat, **kwargs):
    # transform all values to str
    file_pat = str(file_pat)
    kwargs = {
        key : [str(value) for value in values]
        for key, values in kwargs.items()
    }
    key_combines = [dict(zip(kwargs.keys(), combination)) for combination in itertools.product(*kwargs.values())]
    good = True
    for key_combine in key_combines:
        file = file_pat.format(**key_combine)
        if Path(file).exists():
            pass
        else:
            print(f"Waring: {file} not exists")
            good = False
    return good
def check_key(file_pat, key_to_check, **kwargs):
    """
    Check one key against other keys in the file pattern.
    For exmaple, if key_to_check is "sample", pick sample
    that has all files that are expanted from file_pat using
    other keys in kwargs.
    Input:
        file_pat: file pattern
        key_to_check: key to check, one of the keys in kwargs, str
        kwargs: other keys, dict
    Output:
        list of values for key_to_check
    """
    # transform all values to str
    file_pat = str(file_pat)
    sub_dict = {
        key : [str(value) for value in values]
        for key, values in kwargs.items()
        if key != key_to_check
    }

    input_values = kwargs[key_to_check]
    good_values = []
    for input_value in input_values:
        #print(input_value)
        sub_dict.update({key_to_check: [input_value]})
        if check_input(file_pat, **sub_dict):
            good_values.append(input_value)
        else:
            print(f"Warning: {input_value} don't have all files")
    return good_values
# --- helpers --- #
def json_reader(file):
    with open(file, "r") as f:
        data = json.load(f)
    return data