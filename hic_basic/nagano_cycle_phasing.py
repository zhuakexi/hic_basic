import pandas as pd
import os
import gzip
import numpy as np
from sklearn import preprocessing
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
import sklearn.neighbors as NN
from scipy.stats import entropy

def contact_describe(cell_name:str,c1=1,p1=2,c2=3,p2=4) -> pd.Series:
    # get cell's basic statistics, defined in Nagano2017
    contacts = pd.read_table(cell_name, header=None, comment="#",low_memory=False)
    new_columns = list(contacts.columns)
    new_columns[c1] = "chr1"
    new_columns[p1] = "pos1"
    new_columns[c2] = "chr2"
    new_columns[p2] = "pos2"
    contacts.columns = new_columns
    intra = contacts.query(' chr1 == chr2 ')
    distances = abs(intra["pos1"] - intra["pos2"])
    
    all_ = len(distances[23_000 < distances])
    short = len(distances[(23_000 < distances) & (distances < 2_000_000)])
    mitotic = len(distances[(2_000_000 < distances) & (distances < 12_000_000)])
    farAvg = distances[(4_500_000 < distances) & (distances < 225_000_000)]
    
    mitotic_r = mitotic/all_
    short_r = short/all_
    
    # assign to different stages on Peter's cirtera
    if mitotic_r >= 0.3 and short_r <= 0.5:
        group = "Post-M"
    elif short_r > 0.5 and short_r + 1.8*mitotic_r > 1.0:
        group = "Pre-M"
    elif short_r <= 0.63:
        group = "G1"
    elif 0.63 < short_r <= 0.785:
        group = "early/mid-S"
    elif short_r > 0.785:
        group = "mid-S/G2"
    else:
        group = "blank"
    
    return {"short%":short_r, "mitotic%":mitotic_r, "farAvg":farAvg.mean(),"group":group }
"""
def smooth_group(cdps,raw_groups):
    # fine-tuned group assignment using k-means clustering voting
    # Input:
    #   cdps: contacts decay profile, row for a sample
    #         n*m array
    #   groups: original group marker, length same as cdps row number
    #           n array
    # Output:
    #    smooth_group: new group marker
    #                  n array
    #    reassignment_ratio: percent of changing-group sample
    pass
"""
def smooth_cluster(cdps,annote,cluster_col="named_cluster",batchsize=5,pca_n=6):
    """
    # fine-tuned group assignment using k-means clustering voting
    Input:
        cdps: contact decay profiles
        annote: df that records original group assignment, must using sample names as index
        cluster_col: colname in annote that records originmal group assignment
        batchsize: number of samples in each voting batch
        pca_n: number of components used in distance calculation
    Output:
        annote with 2 new col: km_label km_label_group
    """
    # preprocessing
    scaler = preprocessing.StandardScaler()
    scaled = scaler.fit_transform(cdps.loc[annote.index].values)
    pca = PCA(whiten=True)
    pca_res = pca.fit_transform(scaled)
    # kmean clustering
    k_means = KMeans(init="k-means++", n_clusters = annote.shape[0]//batchsize)
    k_means.fit(pca_res[:,:pca_n])
    annote_ = annote.assign(km_label = k_means.labels_)
    # vote for new group
    km_label_group = annote_.groupby("km_label").apply(lambda x: x[cluster_col].value_counts().idxmax())
    km_label_group.name = "km_label_group"
    annote_ = pd.merge(
        annote_, 
        km_label_group,
        left_on="km_label",
        right_index = True
    )
    return annote_
def Nagano_ordering(metrics:pd.DataFrame):
    # ordering within group using Nagano2017's metric
    # Input:
    #   metrics: dataframe, must contain 
    #            short% mitotic% farAvg repli_score groups
    #            
    #            groups: group markers, same length with metric list
    #            factors within Pre-M, mid-S/G2, early/mid-S,
    #            Post-M, G1, and blank
    # Output:
    #    dataframe with new float col named ordering_value
    metrics_e = metrics.assign(
        scaled_short = (metrics["short%"] - metrics["short%"].mean())/metrics["short%"].std(),
        scaled_farAvg = (metrics["farAvg"] - metrics["farAvg"].mean())/metrics["farAvg"].std(),
        shallow_short = metrics["short%"]/metrics["short%"].var(),
        shallow_rs = metrics["repli_score"]/metrics["repli_score"].var()
    )
    ordering_value = []
    for index, row in metrics_e.iterrows():
        if row["group"] == "Post-M":
            ordering_value.append(
                -row["mitotic%"]
            )
        elif row["group"] == "G1":
            ordering_value.append(
                row["scaled_short"] + row["scaled_farAvg"]
            )
        elif row["group"] == "early/mid-S":
            ordering_value.append(
                row["shallow_short"] + row["shallow_rs"]
            )
        elif row["group"] == "mid-S/G2":
            ordering_value.append(
                row["shallow_short"] - row["shallow_rs"]
            )
        elif row["group"] == "Pre-M":
            ordering_value.append(
                row["mitotic%"]
            )
        elif row["group"] == "blank":
            ordering_value.append(np.nan)
    return metrics.assign(ordering_value = ordering_value)
def order_sample_within_group(chunk,shift):
    # Input
    # Output: new chunk with ordering col
    new_chunk = chunk.sort_values("ordering_value").assign(
        index_order = range(0,chunk.shape[0]))
    new_chunk["index_order"] = new_chunk["index_order"] + shift
    return new_chunk
def order_sample(ordering_values):
    # order
    # Input:
    #    ordering_values: dataframe, must have
    #      sample_name, group, ordering_value, cycle_num
    # Output:
    #    dataframe with new int col named index_order
    group_order = ["Post-M","G1","early/mid-S","mid-S/G2","Pre-M"]
    group_order = [group for group in group_order if group in ordering_values["group"].unique()]
    shift_map = list(ordering_values["group"].value_counts().loc[group_order].cumsum().values[:-1])
    shift_map = dict(zip(group_order,[0] + shift_map))
    index_order = []
    output = []
    for group, chunk in ordering_values.groupby("group"):
        if group != "blank":
            output.append(
                order_sample_within_group(chunk,shift_map[group])
            )
        else:
            output.append(
                chunk.assign(index_order=np.nan)
            )
    return pd.concat(output)

def schic_spectral_embedding(df):
    """
    non-linear dimensionality reduction of cdps curve
    Input:
        df: sample x feature matrix
    Output:
        2nd 3rd smallest eigenvectors of the graph Laplacian
    """
    # extract features
    y = df.values.astype(np.float)
    y[y==0] = 1 # add a pseudocount
    y = y/y.sum(axis=1)[:,np.newaxis]

    # compute distances
    celln = y.shape[0]
    dists = np.zeros((celln, celln))
    for i in range(celln):
        for j in range(i+1, celln):
            dists[i,j] = entropy(y[i,:], qk=y[j,:], base=None)
            dists[j,i] = entropy(y[j,:], qk=y[i,:], base=None)

    symdist = dists+dists.T

    # construct NN graph

    k_nn = 7
    nn = NN.NearestNeighbors(n_neighbors=k_nn, metric='precomputed', n_jobs=-1)
    nn.fit(symdist)
    dist2neighs, neighs =  nn.kneighbors()

    adj = np.zeros((symdist.shape[0], symdist.shape[0]))
    for i in range(symdist.shape[0]):
        for ki in range(k_nn):
            adj[i,neighs[i,ki]] = 1
    adj = np.maximum(adj, adj.T)

    # compute spectral embedding

    lap = np.diag(adj.sum(axis=1))-adj
    vals, vecs = np.linalg.eig(lap)
    tmp = np.argsort(vals)
    vals = vals[tmp]
    vecs = vecs[:,tmp]
    
    return pd.DataFrame(vecs[:,1:3],columns=["ev1","ev2"])