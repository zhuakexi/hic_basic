import re
from pathlib import Path

import numpy as np
import pandas as pd
from hires_utils import chrom_rm_suffix


### --- functions about 1D genome --- ###


def parse_ucsc_region(region_str):
    """
    Parses a UCSC-style region string and returns a tuple of (chromosome, start, end).
    
    Parameters:
        region_str (str): A UCSC-style region string like 'chr1:2,345-6,789' or 'chrX:2345-6789'.
    
    Returns:
        tuple: A tuple containing the chromosome (str), start position (int), and end position (int).
    
    Raises:
        ValueError: If the input string is not in a valid UCSC region format.
    """
    # Regex to match UCSC style region strings
    # chr1, chrX, chr1(mat), chr1_mat
    chrom = r"([\w,\(,\)]+)"
    # 1, 1,000, 1000
    pos = r"(\d{1,3}(?:,\d{3})*|\d+)"

    cas_pats = []
    # chr1:1000-2000, chr1
    cas_pats.append(f"^{chrom}(?::{pos}-{pos})?$")
    # chr1:-1000, all chr1 regions before 1000
    cas_pats.append(f"^{chrom}:()-{pos}$")
    # chr1:1000-, all chr1 regions after 1000
    cas_pats.append(f"^{chrom}:{pos}-()$")
    for pat in cas_pats:
        match = re.match(pat, region_str)
        if match is not None:
            # match earlier pattern
            break
    if match is None:
        raise ValueError("Invalid UCSC region format.")
    
    chrom, start, end = match.groups()
    
    # Remove commas if present and convert to integers
    if start is not None:
        start = int(start.replace(',', '')) if len(start) > 0 else None
    if end is not None:
        end = int(end.replace(',', '')) if len(end) > 0 else None
    
    return (chrom, start, end)
class Region:
    """
    A class to represent a genomic region, will store a standard inner representation:
        [(chrom1,pos1), (chrom2,pos2)] no matter whether chrom1 == chrom2
    Relief burden on region-related argument parsing.
    """
    def __init__(self, region_arg, genome=None, binsize=None):
        """
        Initialize a Region object.
        Input:
            region_arg: all usable region representations are supported:
                "chr1", "chr1:1000000-2000000" or "chr1:1000000-":
                    ucsc-style region representation (requires genome if positions incomplete)
                [(chrom,pos), (chrom,pos)]:
                    list of tuples, each tuple is a pair of chrom and pos
                Region: another Region instance (will copy its properties)
            genome: genome name (default None). If None, genome-dependent features are disabled.
            binsize: if not None, will generate bins indexes for the region (requires genome)
        """
        # Handle Region instance input
        if isinstance(region_arg, Region):
            # Copy properties from the existing Region instance
            self.region_arg = region_arg.region_arg
            self.genome = genome if genome is not None else region_arg.genome
            self.binsize = binsize if binsize is not None else region_arg.binsize
            self.r = region_arg.r
            self.ir = region_arg.ir
            
            # Recalculate genome-dependent properties if genome changed
            if genome is not None and genome != region_arg.genome:
                from .data import get_chrom_data
                chrom_df = get_chrom_data(genome)
                chrom_idict = dict(zip(chrom_df.index, range(len(chrom_df))))
                
                # Update index representation with new genome
                ichrom1 = chrom_idict[self.r[0][0]]
                ichrom2 = chrom_idict[self.r[1][0]]
                self.ir = [(ichrom1, self.r[0][1]), (ichrom2, self.r[1][1])]
        else:
            # Original initialization for non-Region inputs
            self.region_arg = region_arg
            self.genome = genome
            self.binsize = binsize
            
            # Parse region (handles genome=None case)
            self.r, self.ir = self._parse_region(region_arg, genome)
        
        # Set slice representations (None if genome unavailable)
        self.slice = slice(*self.r) if self.genome is not None else None
        self.islice = slice(*self.ir) if self.ir is not None else None
        
        # Set relevant chromosomes (None if genome unavailable)
        self.region_chroms = self._get_relevant_chromosomes() if self.genome is not None else None
        
        # Generate bins if requested (skipped if genome unavailable)
        self.bins = None
        if self.binsize is not None and self.genome is not None:
            self.bins = self._get_bin_index(self.binsize)

    def _parse_region(self, region_arg, genome):
        """
        Parse region argument to standard representations.
        Returns:
            r: [(chrom_str, start), (chrom_str, end)]
            ir: [(chrom_idx, start), (chrom_idx, end)] or None if genome unavailable
        """
        if genome is not None:
            # Genome available - use chromosome info
            from .data import get_chrom_data
            chrom_df = get_chrom_data(genome)
            chrom_idict = dict(zip(chrom_df.index, range(len(chrom_df))))
            
            if isinstance(region_arg, str):
                chrom, pos1, pos2 = parse_ucsc_region(region_arg)
                ichrom = chrom_idict[chrom]
                # Fill missing positions using chromosome info
                pos1 = pos1 if pos1 is not None else 0
                pos2 = pos2 if pos2 is not None else chrom_df.loc[chrom, "length"]
                return [(chrom, pos1), (chrom, pos2)], [(ichrom, pos1), (ichrom, pos2)]
                
            elif isinstance(region_arg, list):
                assert len(region_arg) == 2, 'Only "two positions" is supported.'
                chrom1, chrom2 = region_arg[0][0], region_arg[1][0]
                ichrom1, ichrom2 = chrom_idict[chrom1], chrom_idict[chrom2]
                # Fill missing positions using chromosome info
                pos1 = region_arg[0][1] if len(region_arg[0]) > 1 else 0
                pos2 = region_arg[1][1] if len(region_arg[1]) > 1 else chrom_df.loc[chrom2, "length"]
                return [(chrom1, pos1), (chrom2, pos2)], [(ichrom1, pos1), (ichrom2, pos2)]
                
        else:
            # Genome unavailable - require complete specification
            if isinstance(region_arg, str):
                chrom, pos1, pos2 = parse_ucsc_region(region_arg)
                if None in (pos1, pos2):
                    raise ValueError("Incomplete region specification requires genome")
                return [(chrom, pos1), (chrom, pos2)], None
                
            elif isinstance(region_arg, list):
                assert len(region_arg) == 2, 'Only "two positions" is supported.'
                # Require explicit positions in all tuples
                if any(len(tup) < 2 for tup in region_arg):
                    raise ValueError("Incomplete region specification requires genome")
                chrom1, pos1 = region_arg[0][0], region_arg[0][1]
                chrom2, pos2 = region_arg[1][0], region_arg[1][1]
                return [(chrom1, pos1), (chrom2, pos2)], None
                
        raise NotImplementedError("Only list of tuples or UCSC-style region string is supported.")

    def _get_relevant_chromosomes(self):
        """Get relevant chromosomes (requires genome)"""
        from .data import get_chrom_data
        chrom_df = get_chrom_data(self.genome)
        return chrom_df.loc[self.r[0][0]:self.r[1][0]].index.tolist()

    def _get_bin_index(self, binsize):
        """Get bin indexes (requires genome)"""
        all_bins = GenomeIdeograph(self.genome).bins(binsize, bed=True, order=True)
        all_bins = all_bins.set_index(["chrom", "start"]).sort_index()
        return all_bins.loc[self.r[0]:self.r[1]].index.tolist()

    def __eq__(self, other):
        """Check if two regions have the same genomic coordinates."""
        if not isinstance(other, Region):
            return NotImplemented
        return self.r == other.r

    def __contains__(self, other):
        """Check if a point or region is contained within this region."""
        if isinstance(other, tuple) and len(other) == 2:
            # Point containment: (chrom, pos)
            chrom, pos = other
            # Single chromosome region
            if self.r[0][0] == self.r[1][0]:
                return chrom == self.r[0][0] and self.r[0][1] <= pos <= self.r[1][1]
            # Multi-chromosome region
            if self.genome is None:
                raise ValueError("Multi-chromosome containment requires genome reference")
            from .data import get_chrom_data
            chrom_df = get_chrom_data(self.genome)
            chrom_order = chrom_df.index.tolist()
            start_idx = chrom_order.index(self.r[0][0])
            end_idx = chrom_order.index(self.r[1][0])
            target_idx = chrom_order.index(chrom)
            if not (start_idx <= target_idx <= end_idx):
                return False
            if target_idx == start_idx:
                return pos >= self.r[0][1]
            if target_idx == end_idx:
                return pos <= self.r[1][1]
            return True
        elif isinstance(other, Region):
            # Region containment
            # Check if both endpoints of 'other' are contained in self
            return other.r[0] in self and other.r[1] in self
        return NotImplemented


class GenomeIdeograph:
    def __init__(self, chromosomes):
        """Initialize GenomeIdeograph with chromosome data.
        Input:
            chromosomes: RefGenome object, Chromosomes object, DataFrame, 
                        genome name string, or file path to chromosome lengths
        """
        # Use get_chrom_data to handle all input types
        from .data import get_chrom_data
        chrom_df = get_chrom_data(chromosomes, order=False)
        
        # Ensure DataFrame has required format
        assert "length" in chrom_df.columns, "DataFrame must contain 'length' column"
        
        # Ensure index is CategoricalIndex for proper ordering
        if not isinstance(chrom_df.index, pd.CategoricalIndex):
            chrom_df.index = pd.CategoricalIndex(
                chrom_df.index, 
                categories=chrom_df.index, 
                ordered=True
            )
        
        self.chromosomes = chrom_df
    
    def breaks(self, binsize: int, flavor="hickit", within_chromosome=True):
        """
        Get binned reference(int version)
        Input:
            binsize: int
            flavor: str, "hickit", "bedtools", "cooler_compat", "default", or "nearest_boundary"
            within_chromosome: bool, if True, each bin is confined within a chromosome;
                if False, bins can span multiple chromosomes with boundary handling specified by flavor
        Return:
            If within_chromosome=True: breaks of bins(dict of list)
            If within_chromosome=False: list of regions, each as [(start_chrom, start_pos), (end_chrom, end_pos)]
        Note:
            For hickit-flavored binning, the size of the last bin 
            >= 0.5 * binsize and < 1.5 * binsize. For bedtools-flavored binning, 
            the size of the last bin > 0 and <= binsize.
            For cooler_compat-flavored binning, the size of the last bin
            > binsize is trimmed to binsize.
            When within_chromosome is False, flavor can be:
                - "default": no boundary attachment
                - "nearest_boundary": snap to nearest chromosome boundary
        """
        valid_flavors = ["hickit", "bedtools", "cooler_compat", "default", "nearest_boundary"]
        assert flavor in valid_flavors, f"flavor should be one of {valid_flavors}"
        
        if not within_chromosome and flavor not in ["default", "nearest_boundary"]:
            raise ValueError("When within_chromosome=False, flavor must be 'default' or 'nearest_boundary'")
        
        binsize = int(binsize)
        data = self.chromosomes.iloc[:, 0].to_dict()
        
        if within_chromosome:
            # Original within-chromosome binning logic
            all_breaks = {}
            for chrom in data:
                length = data[chrom]
                breaks = list(range(0, length, binsize))
                
                if flavor in ["hickit", "cooler_compat"]:
                    if breaks and length - breaks[-1] < 0.5 * binsize:
                        breaks.pop()
                    if len(breaks) == 0:  # in case length < 0.5 * binsize
                        breaks = [0]
                
                if flavor != "cooler_compat":
                    breaks.append(length)  # don't forget the rightmost point
                else:
                    if breaks and length - breaks[-1] > binsize:
                        breaks.append(breaks[-1] + binsize)
                    else:
                        breaks.append(length)
                
                all_breaks[chrom] = breaks
            return all_breaks
        else:
            # Cross-chromosome binning logic
            # First create global bins across all chromosomes
            total_length = sum(data.values())
            global_breaks = list(range(0, total_length, binsize))
            global_breaks.append(total_length)  # Ensure we cover the entire genome
            
            # Calculate chromosome offsets
            chrom_offsets = {}
            current_offset = 0
            for chrom in data:
                chrom_offsets[chrom] = current_offset
                current_offset += data[chrom]
            
            # Create a mapping from offset to chromosome
            offset_to_chrom = {}
            for chrom, offset in chrom_offsets.items():
                offset_to_chrom[offset] = chrom
            
            # Convert global breaks to chromosome coordinates
            regions = []
            for i in range(len(global_breaks) - 1):
                start_global = global_breaks[i]
                end_global = global_breaks[i + 1]
                
                # Find start chromosome and position
                start_chrom = None
                start_pos = None
                for chrom, offset in sorted(chrom_offsets.items(), key=lambda x: x[1]):
                    if offset <= start_global < offset + data[chrom]:
                        start_chrom = chrom
                        start_pos = start_global - offset
                        break
                
                # Find end chromosome and position
                end_chrom = None
                end_pos = None
                for chrom, offset in sorted(chrom_offsets.items(), key=lambda x: x[1]):
                    if offset <= end_global < offset + data[chrom]:
                        end_chrom = chrom
                        end_pos = end_global - offset
                        break
                
                # If end_global is exactly at the end of a chromosome
                if end_chrom is None:
                    for chrom, offset in sorted(chrom_offsets.items(), key=lambda x: x[1]):
                        if end_global == offset + data[chrom]:
                            end_chrom = chrom
                            end_pos = data[chrom]
                            break
                
                # Handle boundary snapping if requested
                if flavor == "nearest_boundary":
                    # For start position, check if it's close to a boundary
                    if start_pos > 0 and start_pos < binsize / 2:
                        start_pos = 0
                    
                    # For end position, check if it's close to a boundary
                    if end_pos < data[end_chrom] and (data[end_chrom] - end_pos) < binsize / 2:
                        end_pos = data[end_chrom]
                
                regions.append([(start_chrom, start_pos), (end_chrom, end_pos)])
            
            return regions
    
    def bins(self, binsize: int, bed=False, order=False, flavor="hickit", within_chromosome=True):
        """
        Get binned reference(IntervalIndex version)
        Input:
            binsize: int
            bed: bool, if True, return bed format
            order: bool, if True, sort by chromosome and start position
            flavor: str, binning flavor
            within_chromosome: bool, if True, bins are confined within chromosomes;
                if False, bins can span multiple chromosomes
        Output:
            If within_chromosome=True: intervals of each bin(dict of IntervalIndex)
            If within_chromosome=False: list of regions, each as [(start_chrom, start_pos), (end_chrom, end_pos)]
        """
        if order:
            assert bed, "order only works with bed"
        
        binsize = int(binsize)
        breaks = self.breaks(binsize, flavor=flavor, within_chromosome=within_chromosome)
        
        if within_chromosome:
            bins = {
                chrom: pd.IntervalIndex.from_breaks(
                    breaks[chrom], closed="left",
                    name=chrom, dtype='interval[int64]'
                )
                for chrom in breaks
            }
            
            if bed:
                bins_list = [
                    pd.DataFrame({
                        "chrom": chrom,
                        "start": bins[chrom].left,
                        "end": bins[chrom].right
                    })
                    for chrom in bins
                ]
                bins_df = pd.concat(bins_list)
                bins_df = bins_df.reset_index(drop=True)
                
                if order:
                    bins_df["chrom"] = pd.Categorical(
                        bins_df["chrom"],
                        categories=self.chromosomes.index,
                        ordered=True
                    )
                    bins_df = bins_df.sort_values(["chrom", "start"])
                    bins_df = bins_df.reset_index(drop=True)
                
                return bins_df
            
            return bins
        else:
            # For cross-chromosome bins, return the list of regions directly
            return breaks
    def append_bins(self, 
                    df: pd.DataFrame, 
                    binsize: int, 
                    chr1_col: str = "chr1",
                    pos1_col: str = "pos1", 
                    chr2_col: str = None,
                    pos2_col: str = None,
                    flavor: str = "hickit"):
        """
        Append bins to dataframe.
        
        Input:
            df: input dataframe
            binsize: binsize, see GenomeIdeograph.bins
            chr1_col: column name for first chromosome
            pos1_col: column name for first position
            chr2_col: column name for second chromosome (optional)
            pos2_col: column name for second position (optional)
            flavor: flavor of the binning, see GenomeIdeograph.bins
            
        Output:
            df_e: dataframe with new bin columns
                For pairs mode: adds "chrom1","start1","end1","chrom2","start2","end2"
                For single mode: adds "chrom1","start1","end1"
        """
        bins = self.bins(binsize, flavor=flavor)
        
        # Determine mode based on whether chr2 and pos2 columns are provided
        is_pairs_mode = chr2_col is not None and pos2_col is not None
        
        if is_pairs_mode:
            # Pairs mode - process both sets of coordinates
            chunk_e_list = []
            for (chrom1, chrom2), chunk in df.groupby([chr1_col, chr2_col], observed=True):
                cut1 = pd.cut(chunk[pos1_col], bins[chrom1])
                cut2 = pd.cut(chunk[pos2_col], bins[chrom2])
                
                chunk_e = chunk.assign(
                    chrom1=chrom1,
                    start1=pd.Index(cut1.astype(pd.IntervalDtype())).left,
                    end1=pd.Index(cut1.astype(pd.IntervalDtype())).right,
                    chrom2=chrom2,
                    start2=pd.Index(cut2.astype(pd.IntervalDtype())).left,
                    end2=pd.Index(cut2.astype(pd.IntervalDtype())).right,
                )
                chunk_e_list.append(chunk_e)
            
            df_e = pd.concat(chunk_e_list, axis=0)
            
        else:
            # Single mode - process only first set of coordinates
            chunk_e_list = []
            for chrom1, chunk in df.groupby(chr1_col, observed=True):
                if chrom1 not in bins:
                    print(f"Warning: chromosome {chrom1} not found in bins, skipping.")
                    continue  # Skip chromosomes not in bins
                cut1 = pd.cut(chunk[pos1_col], bins[chrom1])
                
                chunk_e = chunk.assign(
                    chrom1=chrom1,
                    start1=pd.Index(cut1.astype(pd.IntervalDtype())).left,
                    end1=pd.Index(cut1.astype(pd.IntervalDtype())).right,
                )
                chunk_e_list.append(chunk_e)
            
            df_e = pd.concat(chunk_e_list, axis=0)
        
        return df_e
    def coarsen_grouper(self, binsize1:int, binsize2:int, flavor:str="hickit")->pd.Series:
        """
        Get grouper to map from binsize1 to a bigger binsize2.
        Input:
            binsize1: int
            binsize2: int
            flavor: str, "hickit" or "bedtools"
        Output:
            pd.Series, 2-level index (chrom, start) from binsize1 as index,
                tuple (chrom, start) from binsize2 as value
        """
        def new_start(row, chrom_intervals=None):
            #print(row)
            middle = (row["start"] + row["end"]) // 2
            interval = chrom_intervals[row["chrom"]][middle]
            return (row["chrom"], interval.left)
        small_bed = self.bins(binsize1, bed=True, order=True, flavor=flavor)
        big_intervals = self.bins(binsize2, bed=False, order=False, flavor=flavor)
        big_intervals = {
            chrom : intervals.to_series()
            for chrom, intervals in big_intervals.items()
        }
        small_bed = small_bed.assign(
            new_pos = small_bed.apply(
                new_start,
                axis=1,
                chrom_intervals = big_intervals
            )
        )
        small_bed = small_bed.set_index(
            ["chrom","start"]
        ).drop(columns="end")
        return small_bed["new_pos"]
    def pixel_id(self, row, binsize:int, intra=True):
        """
        Get pixel_id from row
        Input:
            row: pd.Series, assume bedpe format
            binsize: int
        Output:
            int
        """
        chrom1, start1, end1, chrom2, start2, end2 = row[:6]
        if intra: # not triu, use full to avoid confusion
            offsets_l1 = self.chromosomes["length"] // binsize
            offsets_l2 = ((offsets_l1) ** 2).cumsum()
            offsets_l2 = pd.Series(
                np.insert(offsets_l2.values, 0, 0)[0:-1],
                index = offsets_l2.index
            )
            offsets_l1, offsets_l2 = offsets_l1.to_dict(), offsets_l2.to_dict()
            pixel_id = offsets_l2[chrom1] + (start1 // binsize) * offsets_l1[chrom1] + (start2 // binsize)
        else: # not implemented
            raise NotImplementedError("inter-chrom not implemented")
        return pixel_id
    def join_pixel_id(self, df, binsize, intra=True, filep=None):
        """
        Add a 'pixel_id' column to the input DataFrame in a vectorized manner.
        Input:
            df: DataFrame, bedpe-like with columns 'chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2'
            binsize: int, size of the bins used in pixelation
            intra: bool, if True, process intra-chromosomal pixels
        Output:
            DataFrame with an additional 'pixel_id' column
        """
        print("binnify: \n",type(df))
        print("binnify: \n", df)
        if filep is not None:
            print("binnify: \n", filep)
        if intra:
            # 计算每个染色体长度的累积和，用于确定每个染色体的像素偏移
            offsets_l1 = (self.chromosomes["length"] // binsize).to_dict()
            cumul_offsets_l2 = ((self.chromosomes["length"] // binsize) ** 2).cumsum()
            offsets_l2 = pd.Series(np.insert(cumul_offsets_l2.values, 0, 0)[:-1], index=cumul_offsets_l2.index).to_dict()
            # 计算 pixel_id
            df = df.assign(
                pixel_id = (df['chrom1'].map(offsets_l2).astype(int)
                            + (df['start1'] // binsize) * (df['chrom1'].map(offsets_l1).astype(int))
                            +(df['start2'] // binsize))
                )
        else:
            raise NotImplementedError("inter-chrom not implemented")

        return df
    def join_positions(self, df, binsize, intra=True):
        """
        Add columns 'chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2' to the input DataFrame in a vectorized manner.
        Input:
            df: DataFrame, bedpe-like with column 'pixel_id'
            binsize: int, size of the bins used in pixelation
        Output:
            DataFrame with additional columns 'chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2'
        """
        if intra:
            count_l1 = (self.chromosomes["length"] // binsize)
            cumul_l2 = ((self.chromosomes["length"] // binsize) ** 2).cumsum()
            offsets_l2 = pd.Series(np.insert(cumul_l2.values, 0, 0)[:-1], index=cumul_l2.index).to_dict()
            l2_boundaries = np.insert(cumul_l2.values, 0, 0)
            l2_labels = cumul_l2.index
            # get chr1
            chrom1 = pd.cut(df["pixel_id"], bins=l2_boundaries, labels=l2_labels, right=False)
            sub_chrom_offsets = df["pixel_id"] - chrom1.map(offsets_l2).astype(int)
            dim1_unit = chrom1.map(count_l1).astype(int) # cost how many dim2 for dim1 to gain 1
            start1 = (sub_chrom_offsets // dim1_unit) * binsize
            end1 = start1 + binsize
            start2 = (sub_chrom_offsets % dim1_unit) * binsize
            end2 = start2 + binsize
            df = df.assign(
                chrom1 = chrom1,
                start1 = start1,
                end1 = end1,
                chrom2 = chrom1,
                start2 = start2,
                end2 = end2
            )
        else:
            raise NotImplementedError("inter-chrom not implemented")
        return df
### --- functions about 2D genome --- ###



class RegionPair:
    """
    A class to represent a pair of genomic regions, will store a standard inner representation:
        [Region1, Region2]
    Relief burden on region-related argument parsing.
    """
    def __init__(self, region_pair_arg, genome=None, binsize=None):
        """
        Initialize a RegionPair object.
        Input:
            region_pair_arg: all usable region pair representations are supported:
                [[(chrom1,pos1),(chrom1,pos2)],[(chrom2,pos1),(chrom2,pos2)]]
                    list of two regions, each region is a list of tuples
                ["chr1:1000000-2000000", "chr2:1000000-2000000"]
                    list of two ucsc-style region representations
        Output:
            region_pair: standard inner representation
                [Region1, Region2]
        """
        self.region_pair_arg = region_pair_arg
        self.genome = genome
        self.binsize = binsize
        self.region_pair = self._parse_region_pair(region_pair_arg)
        self.r = [r.r for r in self.region_pair]
    def _parse_region_pair(self, region_pair_arg):
        """
        Parse region pair argument to standard inner representation.
        Input:
            region_pair_arg: see __init__
            genome: see __init__
        Output:
            region_pair: standard inner representation
                [Region1, Region2]
        """
        assert len(region_pair_arg) == 2, 'Only "two regions" is supported.'
        if all([isinstance(region_arg, Region) for region_arg in region_pair_arg]):
            return region_pair_arg
        return [Region(region_arg, genome=self.genome, binsize=self.binsize) for region_arg in region_pair_arg]


### --- functions for genome spreadsheet --- ###


def index_merge_hom(df):
    """
    Cancel ploidiness in index. You can then merge homologous chromosomes use functions like groupby.
    This will change index inplace.
    Input:
        df: DataFrame with MultiIndex ["chrom", "start"], chom is ploidified, like "chr1(mat)"
    Output:
        a dataframe with modified level 0 index value, like "chr1(mat)" -> "chr1"
    """
    chrom = chrom_rm_suffix(
        df.index.get_level_values(0)
        )
    start = df.index.get_level_values(1)
    new_index = pd.MultiIndex.from_arrays(
        [chrom, start],
        names = df.index.names
    )
    df.index = new_index
    return df
def sort_chrom(df, genome):
    """
    Sort dataframe by natural chromosome order.
    TODO:
        1. ask whether to rm chromint after sorting
        2. ask whether to reset index
    Input:
        df: must have a column named "chrom"
    Output:
        df: sorted by natural chromosome order
    """
    from .data import get_chrom_data
    genome = get_chrom_data(genome, order=True)
    genome = genome.reset_index()
    genome = genome.assign(
        chromint = genome.index
    )
    df = pd.merge(
        df,
        genome,
        on = "chrom",
        how = "left"
    )
    df = df.sort_values("chromint", ascending=True)
    return df