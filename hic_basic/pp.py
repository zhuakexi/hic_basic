import pandas as pd
from sklearn import preprocessing
def standard_scaler(df, axis=0, with_std=False):
    """
    Return scaled dataframe.
    Input:
        df: pandas dataframe
        axis: 0 or 1. 0 to scale each col(col as feature), 1 to scale each row(row as feature).
        with_std: whether doing std scaling. No consensus but 
            usually don't scale it in omic-biology.
    Output:
        scaled dataframe. 
    """
    scaler = preprocessing.StandardScaler(with_std=with_std)
    if axis == 0:
        scaled = scaler.fit_transform(df.values)
        scaled = pd.DataFrame(scaled, index=df.index, columns=df.columns)
    elif axis == 1:
        scaled = scaler.fit_transform(df.values.T)
        scaled = pd.DataFrame(scaled.T, index=df.index, columns=df.columns)
    else:
        raise ValueError("axis should be 0 or 1.")
    return scaled