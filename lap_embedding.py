import sys
import pandas as pd
import numpy as np
from scipy.spatial.distance import euclidean, pdist, squareform
from sklearn.manifold import SpectralEmbedding

def get_distance_matrix(cell_cdps:pd.DataFrame)->pd.DataFrame:
    # generate distance_matrix of cells from their contacts decay profiles
    distance_matrix = pd.DataFrame(squareform(pdist(cell_cdps.T))) # input is 'column for samples', transpose to row for samples
    # add columns and index for matrix
    distance_matrix.columns = cell_cdps.columns
    distance_matrix = distance_matrix.set_index(cell_cdps.columns)
    return distance_matrix
def circle_arctan(cordinate:pd.Series)->float:
    # recover ordering from 2d plot
    # 0 to 2pi; circle; 0 or 2pi is about M phase
    try:
        if cordinate.iloc[0] < 0:
            return 1.5*np.pi - np.arctan(cordinate.iloc[1]/cordinate.iloc[0])
        else:
            return 0.5*np.pi - np.arctan(cordinate.iloc[1]/cordinate.iloc[0])
    except ZeroDivisionError:
        if cordinate.iloc[1] == 0:
            return pd.NA
        elif cordinate.iloc[1] > 0:
            return 0.0
        elif cordinate.iloc[1] < 0:
            return np.pi
# test circle_arctan
#pd.DataFrame([[1,1],[1,-1],[-1,1],[-1,-1],[0,1],[0,-1],[0,0]]).apply(circle_arctan,axis=1)/np.pi
def spectral_ordering(distance_matrix:pd.DataFrame)->pd.DataFrame:
    # do spectral embedding and recover ordering
    # row for sample
    embedding = SpectralEmbedding(n_components=2, affinity='rbf')
    cdp_embedding_res = embedding.fit_transform(distance_matrix)
    cpd_embedding_res = pd.DataFrame(cdp_embedding_res, columns=["x","y"])
    cpd_embedding_res.index = distance_matrix.columns
    
    lap_index = cpd_embedding_res.apply(circle_arctan, axis=1)/np.pi #is lap index in fact
    lap_index.index = distance_matrix.columns
    lap_index.name = "lap_index"
    
    indexed_cells = pd.concat([cpd_embedding_res,lap_index],axis=1)
    return indexed_cells.reset_index()
    
