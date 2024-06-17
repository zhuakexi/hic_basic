import re

from .data import chromosomes
from .binnify import GenomeIdeograph
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
        [(chrom,pos), (chrom,pos)]
    Relief burden on region-related argument parsing.
    """
    def __init__(self, region_arg, genome="GRCh38", binsize=None):
        """
        Initialize a Region object.
        Input:
            region_arg: all usable region representations are supported:
                "chr1", "chr1:1000000-2000000" or "chr1:1000000-":
                    ucsc-style region representation
                [(chrom,pos), (chrom,pos)]:
                    list of tuples, each tuple is a pair of chrom and pos
                ...
                TODO: add more region representations
            genome: genome name, used to determine start and end of region
            binsize: if not None, will generate bins indexes for the region
                [(chrom, pos), ...]
        """
        self.region_arg = region_arg
        self.genome = genome
        self.r, self.ir = self._parse_region(region_arg, genome)
        self.slice = slice(*self.r)
        self.islice = slice(*self.ir)
        self.region_chroms = self._get_relevant_chromosomes()
        if binsize is not None:
            self.bins = self._get_bin_index(binsize)
    def _parse_region(self, region_arg, genome):
        """
        Parse region argument to standard inner representation.
        Input:
            region_arg: see __init__
            genome: see __init__
        Output:
            region: standard inner representation
                [(chrom,pos), (chrom,pos)]
        """
        if isinstance(region_arg, str):
            chrom, pos1, pos2 = parse_ucsc_region(region_arg)
            return [(chrom, pos1), (chrom, pos2)], [(chrom, pos1), (chrom, pos2)]
        elif isinstance(region_arg, list):
            assert len(region_arg) == 2
            chrom_df = chromosomes(genome)
            chrom_idict = dict(zip(chrom_df.index, range(len(chrom_df))))
            chrom1, chrom2 = region_arg[0][0], region_arg[1][0]
            ichrom1, ichrom2 = chrom_idict[chrom1], chrom_idict[chrom2]
            pos1 = region_arg[0][1] if len(region_arg[0]) > 1 else 0
            pos2 = region_arg[1][1] if len(region_arg[1]) > 1 else chrom_df.loc[chrom2, "length"]
            return [(chrom1, pos1), (chrom2, pos2)], [(ichrom1, pos1), (ichrom2, pos2)]
        else:
            raise NotImplementedError("Only list of tuples or UCSC-style region string is supported.")
    def _get_relevant_chromosomes(self):
        """
        Get relevant chromosomes for the region.
        Input:
            genome: genome name
        Output:
            chrom_df: relevant chromosomes
                DataFrame
        """
        chrom_df = chromosomes(self.genome)
        chroms = chrom_df.loc[self.r[0][0]:self.r[1][0]].index.tolist()
        return chroms
    def _get_bin_index(self, binsize):
        """
        Get bin indexes for the region.
        Input:
            binsize: binsize
        Output:
            bins: bin indexes for the region
                [(chrom, pos), ...]
        """
        all_bins = GenomeIdeograph(self.genome).bins(
            binsize,bed=True,order=True)
        all_bins = all_bins.set_index(["chrom", "start"]).sort_index()
        return all_bins.loc[self.r[0]:self.r[1]].index.tolist()
class RegionPair:
    """
    A class to represent a pair of genomic regions, will store a standard inner representation:
        [Region1, Region2]
    Relief burden on region-related argument parsing.
    """
    def __init__(self, region_pair_arg, genome="GRCh38", binsize=None):
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
        return [Region(region_arg, genome=self.genome, binsize=self.binsize) for region_arg in region_pair_arg]