# transform pairs table to contact counts between binned coordinates
# pos-pos -> bin-bin / point -> pixel
import pandas as pd, numpy as np
#import datashader as ds
#import datashader.transfer_functions as tf
import pickle as pkl
#import xarray as xr
from .data import chromosomes

class GenomeIdeograph:
    def __init__(self, ref):
        """
        ref: reference_abbrevations(file stored in module) 
            or lengths file_path(csv format:: chrom, lengths)
        """
        self.chromosomes = chromosomes(ref)
    def breaks(self, binsize:int):
        # Get binned reference(int version)
        ##  Return:
        ##    breaks of bins(dict of list)
        binsize = int(binsize)
        data = self.chromosomes.iloc[:,0].to_dict()
        all_breaks = {}
        for chrom in data:
            length = data[chrom]
            breaks = list(range(0, length, binsize))
            breaks.append(length) # don't forget the rightmost point
            all_breaks[chrom] = breaks
        return all_breaks
    def bins(self, binsize:int, bed=False):
        """
        Get binned reference(IntervalIdex version)
        Input:
            binsize: int
            bed: bool, if True, return bed format
        Output:
           intervals of each bin(
                dict of IntervalIndex)        
        """

        binsize = int(binsize)
        breaks = self.breaks(binsize)
        bins = {chrom : pd.IntervalIndex.from_breaks(
                breaks[chrom], closed="left",
                name=chrom,dtype='interval[int64]')
                for chrom in breaks}
        if bed:
            bins = [
                pd.DataFrame(
                    {
                        "chrom" : chrom,
                        "start" : bins[chrom].left,
                        "end" : bins[chrom].right
                    }
                )
                for chrom in bins
            ]
            bins = pd.concat(bins)
        return bins
# def symmetry(X):
#     # flip lower-triangle-part and add 
#     #  it to upper-triangle
#     #  set lower triangle to zeros
#     # Input:
#     #  X: 2darray, assume square
#     X = np.tril(X,-1).T + X
#     X[np.tril_indices(X.shape[0],-1)] = 0
#     return X
# def mat_coarsen(mat, coarseness):
#     """
#     # coarsen matrix to lower resolution
#     # Input:
#     #    mat: matrix
#     # Output:
#     #    matrix
#     """
#     shape = np.array(mat.shape, dtype=float)
#     new_shape = coarseness * np.ceil(shape/ coarseness).astype(int)
#     # zero-padded array
#     zp_mat = np.zeros(new_shape)
#     zp_mat[:mat.shape[0],:mat.shape[1]] = mat
#     temp = zp_mat.reshape(
#         (   zp_mat.shape[0] // coarseness,
#             coarseness,
#             zp_mat.shape[1] // coarseness,
#             coarseness
#         )
#     )
#     zp_mat_c = np.sum(temp,axis=(1,3))
#     return zp_mat_c
def bin_cut(dat:pd.DataFrame, breaks:dict, bins:dict):
    # Binnify contacts between a pair of chromosomes(chr_pair)
    # Input:
    ##  dat: pairs, assume intra-contacts or 
    ##    inter-contacts between two chromosome
    ##  chr_pair: set with 1(for intra) or 2(inter) elements
    ##  breaks: binned reference 
    ##    chromosome_name : boundary of each bin in that chromosome} 
    ##    (keys should contain all eles in chr_pair)
    ##  bins: binned reference(IntervalIndex version)
    ##    use as index
    ##    chromosome_name : intervals of each bin in that chromosome} 
    ##    (keys should contain all eles in chr_pair)
    # Output:
    ##  pd.DataFrame with full interval_index
    
    # using first row to infer which chr_pair this
    chr1, chr2 = dat.iloc[0,[1,3]]
    # binnify
    b_dat, xi, yi = np.histogram2d(x=dat["pos1"],y=dat["pos2"],
        bins=[breaks[chr1],breaks[chr2]])
    # store in sparse matrix
    if chr1 == chr2:
        # upper-triangle for intra_contacts
        b_dat = symmetry(b_dat)
    b_dat = pd.DataFrame(
        b_dat).astype(pd.SparseDtype(int,0))
    # using Interval version of bins as index
    b_dat.index, b_dat.columns = \
        bins[chr1], bins[chr2]
    return b_dat
# def tiled_bin_cut(pairs:pd.DataFrame,ref:GenomeIdeograph,binsize:int)->Hicmap:
#     pairs_b = {}
#     for indi, dat in pairs.groupby(["chr1","chr2"]):
#         pairs_b[indi] = \
#         bin_cut(dat,ref.breaks(binsize),ref.bins(binsize))
#     # in-case input pairs isn't upper-triangle
#     norm_pairs_b = {}
#     for key in pairs_b:
#         if frozenset(key) in norm_pairs_b:
#             norm_pairs_b[frozenset(key)] += \
#             pairs_b[key]
#         else:
#             norm_pairs_b[frozenset(key)] = pairs_b[key]
#     return norm_pairs_b
# def shader_matrix_plot(mat:pd.DataFrame, width:int=500, height:int=500, short_length:int=0):
#     # plot binnified pairs matrix, difference samples aligned in same coords
#     # short_legth set short edge of resulting image, and keep h_w ratio(according to mat.shape)
#     #    useful when mat isn't standard square
#     x_mat = xr.DataArray(mat.values, coords=[("pos1",mat.index),("pos2",mat.columns)]) # transform to xarray form
#     if short_length != 0:
#         expand_r = short_length // min(*mat.shape)
#         cvs = ds.Canvas(plot_width=mat.shape[1]*expand_r, plot_height=mat.shape[0]*expand_r)
#     else:
#         cvs = ds.Canvas(plot_width=width, plot_height=height)
#     x_mat['_file_obj'] = None # work around for ds bug
#     return tf.shade(cvs.raster(x_mat))
# def write_matrix(mat:scipy.sparse.csc_matrix,file_name:str):
#     with open(file_name,'wb') as f:
#         pkl.dump(mat, f)
# def read_matrix(file_name:str)->pd.DataFrame:
#     with open(file_name,'rb') as f:
#         mat = pkl.load(f)
#     return pd.DataFrame.sparse.from_spmatrix(mat)
# class Hicmap:
#     def __init__(self, df, chromosomes, ref, binsize):
#         # df: sparse_matrix
#         # chromosomes: list
#         # ref: word or path str
#         # binsize: int
#         self.df = df.copy(deep=True)
#         self.chromosomes = chromosomes
#         self.binsize = binsize
#         if ref in ["hg19","hg19.dip","hg38","GRCh38"]:
#             self.ref = load_chrom_length()[ref]
#         else:
#             self.ref = pd.read_csv(ref)
#         self.bins = self.get_bins()
#         self.chrom_shifts = get_bin_shifts(self.bins)
#         self.len_bp = self.get_base_length()
#     def get_bins(self):
#         chroms_intervals = {} # each chrom
#         for chrom in self.chromosomes:
#             length = self.ref.loc[chrom,"length"]
#             breaks = list(range(0, length, self.binsize))
#             breaks.append(length) # don't forget the rightmost points
#             chroms_intervals[chrom] = pd.IntervalIndex.from_breaks(
#                 breaks,closed="left",name=chrom,dtype='interval[int64]')
#         return chroms_intervals # better in ordered dict
#     def __len__(self):
#         # caculate total number of all chromosomes' bins
#         return len(self.df)
#     def __str__(self):
#         return self.df.__str__()
#     def get_base_length(self):
#         length = 0
#         for chrom in self.chromosomes:
#             length += self.ref.loc[chrom].values[0]
#         return length
#     def pos2index(self, chrom, pos):
#         b_pos = pd.cut([pos], self.bins[chrom])
#         b_int_pos = b_pos.rename_categories(range(0,len(b_pos.categories)))
#         return self.chrom_shifts[chrom] + b_int_pos[0]
#     def loc(self, chrom1, start1, end1, chrom2, start2, end2):
#         start1_i = self.pos2index(chrom1, start1)
#         end1_i = self.pos2index(chrom1, end1)
#         start2_i = self.pos2index(chrom2, start2)
#         end2_i = self.pos2index(chrom2, end2)
#         #print(start1_i, end1_i)
#         return self.df.iloc[start1_i: end1_i, start2_i:end2_i]