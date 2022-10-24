# trajectory inference methods
import numpy as np
import pandas as pd
def _circle_arctan(cordinate:pd.Series, center=None)->float:
    """
    Recover angle from 2d plot.
        Input: 
            cordinate: 2d cordinate (accept one point, apply this func to dataframe)
            center: center of circle, default to (0,0)
        Output: 
            scalar, 0 to 2pi
        e.g.
            embedding.apply(_circle_arctan, center=embedding.mean(), axis=1)
    """
    if not center is None:
        #print("Tilting center.")
        cordinate = cordinate + np.array(center)
    try:
        if cordinate.iloc[0] < 0:
            return 1.5*np.pi - np.arctan(cordinate.iloc[1]/cordinate.iloc[0])
        else:
            return 0.5*np.pi - np.arctan(cordinate.iloc[1]/cordinate.iloc[0])
    except ZeroDivisionError:
        if cordinate.iloc[1] == 0:
            # (0,0)
            return pd.NA
        elif cordinate.iloc[1] > 0:
            # y axes upper
            return 0.0
        elif cordinate.iloc[1] < 0:
            # y axes lower
            return np.pi
def angle_pseudotime(root,embedding,direction=1,colname="pseudotime")->pd.DataFrame:
    """
    Transform 2D embedding to pseudotime by naive circular ordering.
    Input:
        root: set what sample/angle to 0.
            str: name of sample, must be in data.index; float: angle value
        df: 2d embedding of samples; dataframe
        direction: 1 for clockwise, -1 for anticlockwise
    Output:
        0 - 1 pseudotime, dataframe with new column
    """
    data = embedding.apply(_circle_arctan, center=embedding.mean(axis=0),axis=1)
    if isinstance(root, str):
        rootv = data[root]
    else:
        rootv = float(root)
    data = (data - rootv) % (2*np.pi)
    assert direction in (-1,1)
    if direction == -1:
        data = 2*np.pi - data
    data = data / (2 * np.pi)
    return embedding.assign(**{colname:data})