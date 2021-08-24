# transform pairs table to contact counts between binned coordinates
# pos-pos -> bin-bin / point -> pixel
import pandas as pd, numpy as np
import datashader as ds
import datashader.transfer_functions as tf
import scipy
import pickle as pkl
import xarray as xr
from pkgutil import get_data
from io import StringIO
from . import ref

ref_dat = get_data(ref.__name__, "hg19.len.csv")
ref_f = StringIO(ref_dat.decode())
FULL_CHROM_NAMES=pd.read_csv(ref_f,index_col=0).index
ref_dat = get_data(ref.__name__, "hg19.dip.len.csv")
ref_f = StringIO(ref_dat.decode())
FULL_DIP_CHROM_NAMES=pd.read_csv(ref_f,index_col=0).index
def get_bins(chrom_lengths, chromosomes:list=None, resolution:int=1000000)->dict:
    # generate Intervals for each chromosome in list
    chroms_intervals = {} # each chrom
    for chrom in chromosomes:
        length = chrom_lengths.loc[chrom,"length"]
        breaks = list(range(0, length, resolution))
        breaks.append(length) # don't forget the rightmost points
        chroms_intervals[chrom] = pd.IntervalIndex.from_breaks(breaks,closed="left",name=chrom,dtype='interval[int64]')
    return chroms_intervals # better in ordered dict
#print(get_bins()["chr1"])
def bin_cut(contacts:pd.DataFrame, chrom_pair:tuple, bin_dict:dict):
    # cut single chromosome-pair
    # return int-labeled categories
    # (chra, chra) and (chrb, chra) are regard as different pairs
    chra, chrb = chrom_pair[0], chrom_pair[1]
    b_pos1, b_pos2 = pd.cut(contacts["pos1"], bin_dict[chra]), pd.cut(contacts["pos2"], bin_dict[chrb])
    b_int_pos1, b_int_pos2 = b_pos1.cat.rename_categories(range(0,len(b_pos1.cat.categories))), b_pos2.cat.rename_categories(range(0,len(b_pos2.cat.categories)))
    binned_data = pd.concat([b_int_pos1, b_int_pos2],axis=1)
    binned_data.columns = ["pos1", "pos2"]
    return binned_data
#print(bin_cut(chr1_pos1, ("chr1","chr1"), get_bins()))
# functions to align multiple chromosomes
def get_bin_shifts(chroms_intervals:dict)->pd.Series:
    # shifts for all chfomosome used
    # shifts are used in aligning chromsomes in a single map after per chromosome binning
    interval_nums = {}
    for chrom in chroms_intervals:
        interval_nums[chrom] = len(chroms_intervals[chrom])
    interval_nums = pd.Series(interval_nums)
    bin_shifts = [0]
    bin_shifts.extend(interval_nums.cumsum().values)
    bin_shifts.pop()
    return pd.Series(data=bin_shifts,index=interval_nums.index)
# print(get_bin_shifts(get_bins()))
def get_bin_sum(chroms_intervals:dict)->int:
    # caculate total number of all chromosomes' bins
    interval_nums = {}
    for chrom in chroms_intervals:
        interval_nums[chrom] = len(chroms_intervals[chrom])
    interval_nums = pd.Series(interval_nums)
    return interval_nums.sum()
def shift_binned_chromosomes(pairs_chunk:pd.DataFrame,chrom_pair:tuple, bin_shifts:pd.Series):
    # accept one contacts data block for one chrom combination
    # move binned contacts to align blocks of data on a single map
    # shift according to chromosome used, by default all chromosomes
    chra, chrb = chrom_pair[0], chrom_pair[1]
    new_pos1 = pairs_chunk["pos1"].astype(int) + bin_shifts[chra]
    new_pos2 = pairs_chunk["pos2"].astype(int) + bin_shifts[chrb]
    return pd.concat([new_pos1, new_pos2],axis=1)
def tiled_bin_cut(pairs:pd.DataFrame, chromosomes:list=None, reference:str="hg19", resolution:int=1000000):
    # binnify and align contacts between chromosomes(in chromosome list)
    # by default binnify all chromosomes
    ## get bin breaks
    refs = {
        "hg19":"hg19.len.csv", 
        "hg19.dip":"hg19.dip.len.csv",
        "mm10":"mm10.len.csv",
        "mm10.dip":"mm10.dip.len.csv"
        }
    if reference in refs:
        ref_dat = get_data(ref.__name__, refs[reference])
        ref_f = StringIO(ref_dat.decode())
    else:
        raise ValueError("reference not supported yet")
    ref_content = pd.read_csv(ref_f,index_col=0)
    if chromosomes == None:
        # if chromosome list not given using all chromosomes
        # from the given ref
        chromosomes = ref_content.index
    bin_dict = get_bins(ref_content, chromosomes, resolution)
    bins_shifts = get_bin_shifts(bin_dict)
    bin_sum = get_bin_sum(bin_dict)
    res = pd.DataFrame()
    # a dask parallel groupby will be better
    for group, df in pairs.groupby(['chr1','chr2']): 
        #binnify only chromosomes mentioned in list
        if group[0] in chromosomes and group[1] in chromosomes:
            # bin a single chromosome pair
            b_df = bin_cut(df, group, bin_dict)
            # shift according to chromosomes used
            sb_df = shift_binned_chromosomes(b_df, group, bins_shifts)
            # can append at ease since result is shifted
            res = res.append(sb_df)
    pixels = pd.DataFrame(np.zeros((bin_sum, bin_sum))) # all possible bin-combinations
    l_pixels = pixels.stack().reset_index()[["level_0","level_1"]] # change to longform
    l_pixels.columns = ["pos1","pos2"]
    ex_res = res.append(l_pixels) # extend binnify to ensure all possible combination has at least one value
    ex_bin_counts = ex_res.groupby(["pos1","pos2"]).size() # count bin values
    dex_bin_counts = ex_bin_counts - 1 # remove artificial contacts
    ex_filled = dex_bin_counts.unstack(fill_value=0) # get wide-form/matrix
    sp_ex_filled = scipy.sparse.csc_matrix(ex_filled.values) # transform to column sparse matrix
    return pd.DataFrame.sparse.from_spmatrix(sp_ex_filled)
def shader_pairs_plot(pairs:pd.DataFrame, chromosomes:list=None,ref_file:str="hic_basic/ref/hg19.len.csv",width:int=500, height:int=500):
    # direct plot pairs file, fast but can't align different samples
    ## shift chromosome data tile
    #cvs = ds.Canvas(plot_width=width, plot_height=height)
    #points = cvs.points(pairs, x="pos1",y="pos2")
    #return tf.shade(points) # this will pile all chromosome-pair together, useless
    pass # tedious, not usefull
def shader_matrix_plot(mat:pd.DataFrame, width:int=500, height:int=500, short_length:int=0):
    # plot binnified pairs matrix, difference samples aligned in same coords
    # short_legth set short edge of resulting image, and keep h_w ratio(according to mat.shape)
    #    useful when mat isn't standard square
    x_mat = xr.DataArray(mat.values, coords=[("pos1",mat.index),("pos2",mat.columns)]) # transform to xarray form
    if short_length != 0:
        expand_r = short_length // min(*mat.shape)
        cvs = ds.Canvas(plot_width=mat.shape[1]*expand_r, plot_height=mat.shape[0]*expand_r)
    else:
        cvs = ds.Canvas(plot_width=width, plot_height=height)
    x_mat['_file_obj'] = None # work around for ds bug
    return tf.shade(cvs.raster(x_mat))
def write_matrix(mat:scipy.sparse.csc_matrix,file_name:str):
    with open(file_name,'wb') as f:
        pkl.dump(mat, f)
def read_matrix(file_name:str)->pd.DataFrame:
    with open(file_name,'rb') as f:
        mat = pkl.load(f)
    return pd.DataFrame.sparse.from_spmatrix(mat)