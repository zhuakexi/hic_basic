import gzip
import json
import os
import pickle
import re
import tarfile

from io import StringIO
from pathlib import Path
from itertools import product
from typing import Union, List, Dict, Any

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import pysam
from hires_utils.hires_io import parse_3dg, parse_seg
from scipy.io import mmread
from scipy.sparse import coo_matrix, csr_matrix

from .genome import Region

class DevNull:
    def write(self, _):
        pass

    def flush(self):
        pass

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        pass

def divide_name(filename):
    #home-made os.path.splitext, for it can't handle "name.a.b.c" properly
    basename = os.path.basename(filename)
    parts = basename.split(".") #split return >= 1 length list
    if len(parts) == 1:
        return parts[0], ""
    else:
        return parts[0], "."+".".join(parts[1:])
import pandas as pd

def parse_bed(bed_path: str) -> pd.DataFrame:
    """
    Load and parse BED format files, automatically assigning BED12 column names.
    
    Inputs:
        bed_path (str): Path to the BED file.
        
    Returns:
        pd.DataFrame: Parsed BED data with columns named according to BED12 specification:
            - chrom (str): Chromosome name.
            - start (int): Start coordinate (0-based).
            - end (int): End coordinate (exclusive).
            - name (str, optional): Feature name.
            - score (float, optional): Score value (0-1000).
            - strand (str, optional): Strand ('+' or '-').
            - thickStart (int, optional): Thick start coordinate.
            - thickEnd (int, optional): Thick end coordinate.
            - itemRgb (str, optional): RGB color (e.g., '255,0,0').
            - blockCount (int, optional): Number of blocks.
            - blockSizes (str, optional): Comma-separated block sizes.
            - blockStarts (str, optional): Comma-separated block start positions.
            
    Notes:
        - Automatically detects the number of columns in the BED file
        - Assigns BED12 column names in order: chrom, start, end, name, score, 
          strand, thickStart, thickEnd, itemRgb, blockCount, blockSizes, blockStarts
        - If the file has fewer than 12 columns, only the first N columns are named
        - If the file has more than 12 columns, only the first 12 are read
        - Performs automatic type conversion for numeric columns
        - Coordinates are 0-based and end-exclusive by BED convention
    """
    # Define all BED12 column names in standard order
    bed12_cols = [
        'chrom', 'start', 'end', 'name', 'score', 'strand',
        'thickStart', 'thickEnd', 'itemRgb', 'blockCount',
        'blockSizes', 'blockStarts'
    ]
    
    # First, determine the number of columns in the BED file
    with open(bed_path, 'r') as f:
        first_line = f.readline().strip()
        if not first_line:  # Empty file
            n_cols = 0
        else:
            n_cols = len(first_line.split('\t'))
    
    # Read only up to 12 columns
    n_read = min(n_cols, 12)
    
    # Read the BED file with appropriate column names
    if n_read == 0:
        return pd.DataFrame(columns=bed12_cols)
    
    # Create column names for the actual number of columns
    col_names = bed12_cols[:n_read]
    
    # Read the file
    df = pd.read_csv(
        bed_path, 
        sep='\t', 
        header=None,
        comment='#',  # Skip comment lines
        usecols=range(n_read),  # Read only available columns
        names=col_names
    )
    
    # Define numeric columns and convert them
    numeric_cols = ['start', 'end', 'score', 'thickStart', 'thickEnd', 'blockCount']
    
    for col in numeric_cols:
        if col in df.columns:
            # Convert to numeric, coerce errors to NaN
            df[col] = pd.to_numeric(df[col], errors='coerce')
            # Convert to integer if possible (except score which can be float)
            if col != 'score':
                df[col] = df[col].astype(pd.Int64Dtype())  # Use nullable integer type
    
    return df
def parse_pairs(filename:str)->pd.DataFrame:
    '''
    read from 4DN's standard .pairs format
    compatible with all hickit originated pairs-like format 
    '''
    # read comments
    with gzip.open(filename,"rt") as f:
        comments = []
        chromosomes = []
        lengths = []
        for line in f.readlines():
            if line[0] != "#":
                break
            if line.startswith("#chromosome") or line.startswith("#chromsize"):
                chrom, length = line.split(":")[1].strip().split()
                chromosomes.append(chrom)
                lengths.append(int(length))
            if line.startswith("#columns:"):
                columns = line.split(":")[1].strip().split()
            ## comment lines are stored in dataframe.attrs["comment"]
            comments.append(line)
    dtype_array = {"readID":"category",
            "chr1":pd.CategoricalDtype(categories=chromosomes),
            "pos1":"int",
            "chr2":pd.CategoricalDtype(categories=chromosomes),
            "pos2":"int",
            "strand1":pd.CategoricalDtype(categories=["+","-"]),
            "strand2":pd.CategoricalDtype(categories=["+","-"]),
            "phase0":pd.CategoricalDtype(categories=["1","0","."]),
            "phase1":pd.CategoricalDtype(categories=["1","0","."]),
            "phase_prob00":"float",
            "phase_prob01":"float",
            "phase_prob10":"float",
            "phase_prob11":"float"}
    dtypes = {key:value for key, value in dtype_array.items() if key in columns}
    #read table format data
    pairs = pd.read_table(
        filename, 
        header=None, 
        comment="#",
        dtype=dtypes,
        names=columns
        )
    pairs.attrs["comments"] = comments
    pairs.attrs["name"], _ = divide_name(filename) # infer real sample name
    pairs.attrs["chromosomes"] = chromosomes
    pairs.attrs["lengths"] = lengths
    #assign column names
    #sys.stderr.write("pairs_parser: %s parsed \n" % filename)
    return pairs
def parse_pairs_like(filename:str)->pd.DataFrame:
    '''
    read from 4DN's standard .pairs format
    compatible with all hickit originated pairs-like format 
    '''
    #comment lines are stored in dataframe.attrs["comment"]
    name_array = "readID chr1 pos1 chr2 pos2 strand1 strand2 phase0 phase1 phase_prob00 phase_prob01 phase_prob10 phase_prob11".split()
    dtype_array = {"readID":"category",
            "chr1":"string",
            "pos1":"int",
            "chr2":"string",
            "pos2":"int",
            "strand1":"string",
            "strand2":"string",
            "phase0":"string",
            "phase1":"string",
            "phase_prob00":"float",
            "phase_prob01":"float",
            "phase_prob10":"float",
            "phase_prob11":"float"}
    #read comment line
    with gzip.open(filename,"rt") as f:
        comments = []
        for line in f.readlines():
            if line[0] != "#":
                break
            comments.append(line)
    #infer number of columns
    line_length = len(line.strip().split("\t"))
    #pick used eles from builtin arrays
    columns = name_array[0:line_length]
    dtypes = {key:value for key, value in dtype_array.items() if key in columns}
    #read table format data
    pairs = pd.read_table(
        filename, 
        header=None, 
        comment="#",
        dtype=dtypes,
        names=columns
        )
    pairs.attrs["comments"] = comments
    pairs.attrs["name"], _ = divide_name(filename) # infer real sample name
    #assign column names
    #sys.stderr.write("pairs_parser: %s parsed \n" % filename)
    return pairs
def parse_gtf(file,ID=True,name=True,mgi_id=True, **args):
    """
    Parsing gtf file. Read all in memory. Extract gene_id to df if set true.
    Input:
        file: path to .gtf file.
        ID: whether to parse gene_id from attribute string.
        name: whether to parse gene_name from attribute string.
        **args: arguments passed to pd.read_table
    Return:
        pd.DataFrame
    """
    header = ["seqname","source","feature","start","end","score","strand","frame","attributes"]
    dtypes = dict(zip(header, ["category","category","category","int","int","string","category","category","string"]))
    gtf = pd.read_table(file, comment ="#", header=None, names = header, dtype=dtypes,**args)
    # get gene_id from attributes
    if ID:
        ID = gtf["attributes"].str.extract(r"gene_id \"([\w,.]+)\";", expand=True)
        ID.columns = ["gene_id"]
        gtf = pd.concat([gtf, ID], axis=1, join="inner")
    if name:
        name = gtf["attributes"].str.extract(r"gene_name \"(\S+)\";", expand=True)
        name.columns = ["gene_name"]
        gtf = pd.concat([gtf, name], axis=1, join="inner")
    if mgi_id:
        mgi_id = gtf["attributes"].str.extract(r"mgi_id \"(MGI:\d+)\";", expand=True)
        mgi_id.columns = ["mgi_id"]
        gtf = pd.concat([gtf, mgi_id], axis=1, join="inner")
    return gtf
def parse_gff(file, ID=False, Name=False):
    """
    Parsing gff file. Read all in memory. Extract ID(gene), and Name if set true.
    """
    header = ["seqid","annotation_source","feature_type","start","end","score","strand","phase","attributes"]
    dtypes = dict(zip(header,["category","category","category","int","int","string","category","category","string"]))
    gff = pd.read_table(file,
                        comment="#",
                        header=None, 
                        names = header,
                        dtype = dtypes
                       )
    additions = []
    if ID:
        ID = gff["attributes"].str.extract(r"ID=gene:(\w+);",expand=True)
        ID.columns = ["ID"]
        additions.append(ID)
    if Name:
        Name = gff["attributes"].str.extract(r"Name=(\w+);",expand=True)
        Name.columns = ["Name"]
        additions.append(Name)
    gff = pd.concat([gff] + additions, axis=1, join="inner")
    return gff

def parse_vcf(vcf_path, region, return_alleles=False):
    """
    Extracts genotype data for all samples within a specified genomic region from a VCF file.
    Optionally returns genotypes in allele form (e.g., 'GA/GA' instead of '1/1').

    Parameters:
        vcf_path (str): Path to the VCF file.
        chromosome (str): Target chromosome (e.g., '1', 'chr1').
        start (int): Start position (1-based).
        end (int): End position (1-based).
        return_alleles (bool): If True, returns genotype in allele form. Default is False.

    Returns:
        pd.DataFrame: A DataFrame where each row represents a variant and each column represents a sample.
                      Genotypes are returned either as indices (e.g., '1/1') or as allele strings (e.g., 'GA/GA').
    """
    ((chrom, start), (_, end)) = Region(region).r
    # # --- check query chrom --- #
    # for i in range(10):
    #     chromnames = []
    #     try:
    #         chromnames.append(vcf.get_reference_name(i))
    #     except ValueError:
    #         continue
    # if len(chromnames) == 0:
    #     raise ValueError("No valid chromosome names found in VCF header.")
    # else:
    #     flavor = check_chromname_flavor(chromnames)
    with pysam.VariantFile(vcf_path) as vcf_file:

        samples = list(vcf_file.header.samples)
        genotype_data = []

        for record in vcf_file.fetch(chrom, start, end):

            variant_info = {
                'CHROM': record.chrom,
                'POS': record.pos,
                'ID': record.id or '.',
                'REF': record.ref,
                'ALT': ','.join(record.alts) if record.alts else '.',
                'QUAL': record.qual if record.qual is not None else '.',
                'FILTER': ';'.join(record.filter) if record.filter else 'PASS'
            }

            for sample in samples:
                if return_alleles:
                    alleles = record.samples.get(sample).alleles
                    alleles = tuple(i if i is not None else "N" for i in alleles)
                    genotype_str = '/'.join(alleles)
                else:
                    genotype = record.samples[sample].get('GT', './.')
                    if isinstance(genotype, tuple):
                        genotype_str = '/'.join(map(str, genotype))
                    else:
                        genotype_str = str(genotype)
                variant_info[sample] = genotype_str

            genotype_data.append(variant_info)

        df = pd.DataFrame(genotype_data)
        standard_cols = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER']
        df = df[standard_cols + samples]

    return df


### --- read DNA file --- ###
def parse_sam_line(sam_line: Union[str, List[str]]) -> pd.Series:
    """
    Parse a SAM format line into a pandas Series with standardized fields and optional tags.
    
    This function parses a SAM (Sequence Alignment/Map) format record, extracting both
    the 11 mandatory fields and any optional TAG:TYPE:VALUE fields according to the
    SAM specification v1.6+.
    
    Args:
        sam_line: A SAM format line as a string or pre-split list of strings.
                  If string, it will be split by tabs.
                  
    Returns:
        pd.Series: A pandas Series where:
                   - Index includes: 'QNAME', 'FLAG', 'RNAME', 'POS', 'MAPQ', 'CIGAR',
                     'RNEXT', 'PNEXT', 'TLEN', 'SEQ', 'QUAL' (mandatory fields)
                   - Additional indices for any optional tags (e.g., 'NM', 'MD', 'AS')
                   - Values are appropriately typed (int, float, str) based on SAM type codes
                   
    Raises:
        ValueError: If the input doesn't contain at least 11 fields (mandatory SAM fields)
        
    Notes:
        - POS, PNEXT: 1-based coordinates (0 indicates unmapped)
        - FLAG: Bitwise flag, returned as integer for further decoding
        - Optional tags follow pattern: TAG:TYPE:VALUE
        - TYPE codes: i (int), f (float), Z (string), A (char), H (hex), B (byte array)
        - For 'B' type (byte arrays), values are returned as list of appropriate type
        
    Examples:
        >>> sam_line = "read1\t16\tchr1\t100\t255\t30M\t=\t200\t50\tACGT\tIIII\tNM:i:1\tMD:Z:30"
        >>> series = parse_sam_line(sam_line)
        >>> series['QNAME']  # Returns 'read1'
        >>> series['FLAG']   # Returns 16
        >>> series['NM']     # Returns 1
        >>> series['MD']     # Returns '30'
    """
    # Split line if input is string
    if isinstance(sam_line, str):
        fields = sam_line.strip().split('\t')
    else:
        fields = sam_line
    
    # Validate minimum field count
    if len(fields) < 11:
        raise ValueError(
            f"SAM line must have at least 11 fields, got {len(fields)}. "
            f"Line: {' '.join(fields[:11])}"
        )
    
    # Initialize result dictionary
    result = {}
    
    # --- Parse mandatory fields ---
    # Field 1: QNAME (Query template NAME) - string
    result['QNAME'] = fields[0]
    
    # Field 2: FLAG (Bitwise FLAG) - integer
    result['FLAG'] = int(fields[1])
    
    # Field 3: RNAME (Reference sequence NAME) - string
    # '*' indicates no alignment
    result['RNAME'] = fields[2]
    
    # Field 4: POS (1-based leftmost mapping POSition) - integer
    # 0 indicates unmapped
    result['POS'] = int(fields[3])
    
    # Field 5: MAPQ (MAPping Quality) - integer
    result['MAPQ'] = int(fields[4])
    
    # Field 6: CIGAR (CIGAR string) - string
    result['CIGAR'] = fields[5]
    
    # Field 7: RNEXT (Reference name of the mate/next read) - string
    result['RNEXT'] = fields[6]
    
    # Field 8: PNEXT (Position of the mate/next read) - integer
    result['PNEXT'] = int(fields[7])
    
    # Field 9: TLEN (observed Template LENgth) - integer
    result['TLEN'] = int(fields[8])
    
    # Field 10: SEQ (segment SEQuence) - string
    result['SEQ'] = fields[9]
    
    # Field 11: QUAL (ASCII of Phred-scaled base QUALity+33) - string
    result['QUAL'] = fields[10]
    
    # --- Parse optional tags (fields 12+) ---
    for field in fields[11:]:
        # Skip empty fields
        if not field:
            continue
            
        # Split tag field (format: TAG:TYPE:VALUE)
        parts = field.split(':', 2)
        
        # Validate tag format
        if len(parts) != 3:
            # Some tools may produce malformed tags; store as-is with warning
            import warnings
            warnings.warn(
                f"Optional field '{field}' does not follow TAG:TYPE:VALUE format. "
                f"Storing as raw string."
            )
            result[field] = field
            continue
        
        tag, type_code, value_str = parts
        
        # Parse value based on type code
        try:
            if type_code == 'i':  # Integer
                result[tag] = int(value_str)
            elif type_code == 'f':  # Floating point
                result[tag] = float(value_str)
            elif type_code == 'Z':  # String
                result[tag] = value_str
            elif type_code == 'A':  # Single character
                if len(value_str) != 1:
                    raise ValueError(f"Type 'A' expects single character, got '{value_str}'")
                result[tag] = value_str
            elif type_code == 'H':  # Hex string
                # Store as string; conversion to bytes if needed can be done externally
                result[tag] = value_str
            elif type_code == 'B':  # Byte array (numeric array)
                # Format: B:[cCsSiIf]:[value1,value2,...]
                # First character after 'B:' indicates array element type
                if value_str[0] not in 'cCsSiIf':
                    raise ValueError(f"Invalid byte array type specifier: {value_str[0]}")
                
                # Remove the type specifier and parse comma-separated values
                array_values = value_str[2:] if value_str[1] == ':' else value_str[1:]
                
                # Handle empty array
                if not array_values:
                    result[tag] = []
                    continue
                    
                # Parse based on element type
                elem_type = value_str[0]
                if elem_type in 'cCsS':  # Signed/unsigned byte/short
                    result[tag] = [int(x) for x in array_values.split(',')]
                elif elem_type == 'i':  # Integer
                    result[tag] = [int(x) for x in array_values.split(',')]
                elif elem_type == 'f':  # Float
                    result[tag] = [float(x) for x in array_values.split(',')]
                else:
                    # Fallback: store as string
                    result[tag] = array_values
            else:
                # Unknown type code - store as string
                import warnings
                warnings.warn(
                    f"Unknown type code '{type_code}' for tag '{tag}'. "
                    f"Storing value as string."
                )
                result[tag] = value_str
                
        except (ValueError, IndexError) as e:
            # If parsing fails, store raw string with error indication
            import warnings
            warnings.warn(
                f"Failed to parse tag '{tag}:{type_code}:{value_str}': {str(e)}. "
                f"Storing as raw string."
            )
            result[tag] = field
    
    # Convert to pandas Series with appropriate dtype detection
    # We'll let pandas infer dtypes, but ensure certain fields maintain correct type
    series = pd.Series(result)
    
    # Ensure numeric fields remain numeric (pandas may convert to float if NaN present)
    int_fields = ['FLAG', 'POS', 'MAPQ', 'PNEXT', 'TLEN']
    for field in int_fields:
        if field in series.index and pd.notna(series[field]):
            # Convert to nullable Int64 if NaN possible, else regular int
            if series[field] == 0 and field in ['POS', 'PNEXT']:
                # 0 has special meaning in SAM (unmapped), keep as int
                series[field] = int(series[field])
    
    return series


# Optional utility function for decoding FLAG field
def decode_sam_flag(flag: int) -> Dict[str, bool]:
    """
    Decode SAM bitwise FLAG into human-readable components.
    
    Args:
        flag: Integer FLAG value from SAM field 2
        
    Returns:
        Dictionary with flag component names as keys and boolean values
        
    Reference: SAM v1.6 specification
    """
    flag_dict = {
        'READ_PAIRED': bool(flag & 0x1),
        'PROPER_PAIR': bool(flag & 0x2),
        'READ_UNMAPPED': bool(flag & 0x4),
        'MATE_UNMAPPED': bool(flag & 0x8),
        'READ_REVERSE_STRAND': bool(flag & 0x10),
        'MATE_REVERSE_STRAND': bool(flag & 0x20),
        'FIRST_IN_PAIR': bool(flag & 0x40),
        'SECOND_IN_PAIR': bool(flag & 0x80),
        'NOT_PRIMARY_ALIGNMENT': bool(flag & 0x100),
        'READ_FAILS_QC': bool(flag & 0x200),
        'DUPLICATE_READ': bool(flag & 0x400),
        'SUPPLEMENTARY_ALIGNMENT': bool(flag & 0x800)
    }
    return flag_dict



### --- read scRNA-seq file --- ###
def read_expr(path,sep=None)->pd.DataFrame:
    """
    Read expression matrix from file.
    Don't guarantee OV(obs*var) or VO(var*obs) output. It just read in the file.
    Input:
        path: path to matrix file
        sep: separator, if None will infer from file extension
    Output:
        expression matrix
    """
    mat = None

    # in case a dataframe is passed in
    if isinstance(path, pd.DataFrame):
        mat = path
        return mat

    # formats that don't require separator inference
    # like .parquet, .pkl, .pkl.gz
    if path.endswith("parquet"):
        mat = pd.read_parquet(path)
    elif path.endswith("pkl") or path.endswith("pkl.gz"):
        mat = pd.read_pickle(path)
    else:
        pass
    
    if mat is not None:
        return mat

    # format that requires separator inference
    if sep is None:
        if path.endswith("csv.gz") or path.endswith("csv"):
            mat = pd.read_csv(path,index_col=0)
        elif path.endswith("tsv.gz") or path.endswith("tsv"):
            mat = pd.read_table(path,index_col=0)
        else:
            print("Unknown file extension, will assume tab separated")
            mat = pd.read_table(path,index_col=0)
    else:
        mat = pd.read_table(path,sep=sep,index_col=0)
        if mat.shape ==(0,0):
            # purge index and column if no content
            mat = pd.DataFrame()
    mat.columns = mat.columns.astype("string")
    mat.index = mat.index.astype("string")
    return mat
def match_10x(file_list:list, retind=False)->dict:
    """
    Check if input file names contains a 10x result.
    Input:
        file_list: list of file path
        retind: return index of matched file instead of file path
    Output:
        {"matrix":matrix_file, "genes":gene_file, "barcodes":barcodes_file}
    """
    checker = {
        "matrix" : ["matrix"],
        "genes" : ["genes","features"],
        "barcodes" : ["barcodes"],
    }
    scores = {}
    result = {}
    file_list = [str(file) for file in file_list] # also accept iterable
    for filetype, keywords in checker.items():
        # score for each file
        scores = [
            sum([re.search(keyword, str(file).lower()) is not None for keyword in keywords])
            for file in file_list]
        # find file with highest score
        best_score = max(scores)
        if best_score > 0:
            if retind: # return index
                result[filetype] = scores.index(best_score)
            else: # return file path
                result[filetype] = file_list[scores.index(best_score)]
        else:
            result[filetype] = None
            print("Warning: No file found for %s" % filetype)
    return result
def combine_10x(annoted_file:dict)->pd.DataFrame:
    """
    Extract data from files and combine them into a expression matrix.
    Input:
        annoted_file: {"matrix":matrix_file, "genes":gene, "barcodes":barcodes}
    Output:
        a cell * gene matrix
    """
    pass
def read_10x(fp, gene_column=1, cell_column=0)->pd.DataFrame:
    """
    Read 10x-flavor matrix market file.
    Input:
        fp: file path
        gene_column: use this column from gene/feature matrix
        cell_column: use this column from cell matrix
    Output:
        a cell * gene matrix
    """
    fp = str(fp) # use end to decide input type
    if fp.endswith(".tar.gz"):
        with tarfile.open(fp, "r:gz") as tar:
            members = tar.getmembers()
            filenames = [member.name for member in members]
            annoted_members = {
                filetype : members[ind]
                for filetype, ind in match_10x(filenames, retind=True).items()
            }
            #print(annoted_members)
            for filetype, member in annoted_members.items():
                with tar.extractfile(member) as gf:
                    with gzip.open(gf,"rt") as f:
                        if filetype == "matrix":
                            data = mmread(f)
                        elif filetype == "genes":
                            genes = pd.read_table(
                                f,header=None)[gene_column].rename("var")
                        elif filetype == "barcodes":
                            barcodes = pd.read_table(
                                f,header=None)[cell_column].rename("obs")
                        else:
                            raise ValueError("Unknown file type: %s" % filetype)
        result = pd.DataFrame(data.todense(), index=genes, columns=barcodes)
    else:
        # treat as directory
        print("Treating as directory...")
        file_list = Path(fp).glob("*")
        for filetype, filename in match_10x(file_list, retind=False).items():
            if filename.endswith(".gz"):
                opener = gzip.open
            else:
                opener = open
            with opener(filename) as f:
                if filetype == "matrix":
                    # mmread accept filepath and call gzip in fact
                    # here just for consistency
                    data = mmread(f)
                elif filetype == "genes":
                    genes = pd.read_table(
                        f,header=None)[gene_column].rename("var")
                elif filetype == "barcodes":
                    barcodes = pd.read_table(
                        f,header=None)[cell_column].rename("obs")
                else:
                    raise ValueError("Unknown file type: %s" % filetype)
        result = pd.DataFrame(data.todense(), index=genes, columns=barcodes)
    return result
# read cooler file
def read_h5_ds(group):
    """
    Read .h5 bottom level group into dataframe.
    Input:
        group: h5py.Group; must have same length datasets.
    Output:
        pd.DataFrame
    """
    data = {}
    for key in group.keys():
        if group[key].dtype.type == np.string_:
            data[key] = group[key][:].astype("U")
        else:
            data[key] = group[key][:]
    return pd.DataFrame(data)
def load_cool(cool, root="/"):
    """
    Load simple cooler file as sparsematrix
    
    Parameters
    ----------
    cool : str
        Path to the input .cool file.

    Returns
    -------
    mat : scipy coo_matrix
        Hi-C contact map in COO format.
    frags : pandas DataFrame
        Table of bins matching the matrix.
    chroms : pandas DataFrame
        Table of chromosome informations.
    """
    with h5py.File(cool,"r") as f:
        frags = read_h5_ds(f[root]["bins"])
        n_frags = frags.groupby("chrom", sort=False).count()["start"]
        n_frags.name = "n_frags"
        chroms = read_h5_ds(f[root]["chroms"])
        mat = read_h5_ds(f[root]["pixels"])
    frags["id"] = frags.groupby("chrom", sort=False).cumcount() + 1
    chroms["cumul_length"] = (
        chroms.length.shift(1).fillna(0).cumsum().astype(int)
    )
    chroms = pd.concat([chroms,n_frags],axis=1)    
    n = int(max(np.amax(mat.bin1_id), np.amax(mat.bin2_id))) + 1
    shape = (n, n)
    mat = coo_matrix((mat["count"], (mat.bin1_id, mat.bin2_id)), shape=shape)

    return mat, frags, chroms
# read schilcuster results
def parse_hicluster_res(embed, sample_table):
    """
    Read schicluster embedding result into dataframe.
    Input:
        embed: schicluster concat-cell output hdf5 file
        sample_table: input sample file of schicluster pipeline or list of sample names
    Examples:
        from hic_basic.io import parse_hicluster_res
        from hic_basic.scAB_embedding import do_umap
        hicluster_res = parse_hicluster_res(embed, sample_table)
        umap_res = do_umap(hicluster_res.values)
    """
    if isinstance(sample_table, str):
        sample_table = read_meta(sample_table)
        samples = sample_table.index
    else:
        samples = sample_table
    f = h5py.File(embed, "r")
    data = pd.DataFrame(f["data"][:], index=samples)
    f.close()
    data.columns = ["PC%d" % i for i in range(1, data.shape[1]+1)]
    return data
import h5py
import numpy as np
from scipy.sparse import csr_matrix


def schicluster2mat(filei):
    """
    Read in schicluster impute .hdf5 or np.savez_compressed .npz file to sparse matrix.

    This function supports two file formats:
    1. HDF5 files with standard schicluster sparse matrix structure
    2. NPZ files saved using np.savez_compressed with 'data' key

    Args:
        filei: Path to input file. Supported formats:
               - .hdf5: schicluster imputed HDF5 file
               - .npz: Compressed numpy archive with sparse matrix components

    Returns:
        csr_matrix: Scipy sparse matrix in CSR format

    Raises:
        ValueError: If file format is not supported or required keys are missing
        OSError: If file cannot be opened or read
    """
    if filei.endswith('.hdf5') or filei.endswith('.h5'):
        # Read HDF5 format
        with h5py.File(filei, "r") as f:
            g = f["Matrix"]
            A = csr_matrix(
                (g["data"][()], g["indices"][()], g["indptr"][()]), 
                g.attrs["shape"]
            )
        return A
    
    elif filei.endswith('.npz'):
        # Read compressed NPZ format
        npz_data = np.load(filei, allow_pickle=True)
        
        # Check if it's a sparse matrix saved with individual components
        if all(key in npz_data for key in ['data', 'indices', 'indptr', 'shape']):
            # Reconstruct CSR matrix from components
            A = csr_matrix(
                (npz_data['data'], npz_data['indices'], npz_data['indptr']),
                shape=tuple(npz_data['shape'])
            )
            return A
        
        # Otherwise, check for a single 'data' key with pickled matrix
        elif 'data' in npz_data:
            # Extract the matrix - handle 0-d array case
            mat_data = npz_data['data']
            if mat_data.shape == ():
                # 0-dimensional array, extract the actual object
                A = mat_data.item()
            else:
                # Regular array, convert to csr_matrix
                A = csr_matrix(mat_data)
            return A
        
        else:
            raise ValueError(
                f"NPZ file has unexpected structure. "
                f"Expected either CSR components (data, indices, indptr, shape) "
                f"or a 'data' key with matrix. "
                f"Available keys: {list(npz_data.keys())}"
            )
    
    else:
        raise ValueError(
            f"Unsupported file format: {filei}. "
            "Supported formats: .hdf5, .h5, .npz"
        )
def schiclusterDir2mat(hiclusterdir):
    """
    Read in schicluster imputed .hdf5 files (from a dir) to sparse matrix.
    Input:
        hicluster: schicluster imputed .hdf5 file dir. e.g. f"/shareb/ychi/repo/sperm22_schicluster/imputed_matrix/100000/{chrom}/"
    Output:
        sparse matrix
    """
    Ap = sum(map(schicluster2mat, map(lambda x:os.path.join(hiclusterdir,x),os.listdir(hiclusterdir))))
    #Apc = mat_coarsen(Ap.todense(), coarsen)
    return Ap
# write scipy sparse matrix to 3-col file
def write_triplet(sparseM,filep,max_coo=False,zipping=True):
    """
    Write sparse matrix to triplet txt, assume symmetry.
    For scipy doesn't give a way to write triplet format directly.
    Input:
        sparseM: scipy sparse matrix
        filesp: output file path
        max_coo: if True, write largest possible coordinate at first row
    Output:
        indi, indj, data; separated by tab
    """
    sparseM.maxprint = sparseM.count_nonzero()
    if zipping:
        opener = gzip.open
    else:
        opener = open
    if max_coo == True:
        with opener(filep,"wt") as f:
            f.write("%d\n" % sparseM.shape[0])
        with opener(filep,"at") as f:
            for line in StringIO(str(sparseM)):
                coord, value = line.split("\t")
                c12 = " ".join(coord.strip().strip("()").split(", "))
                f.write(" ".join([c12, value]))
    else:
        with opener(filep,"wt") as f:
            for line in StringIO(str(sparseM)):
                coord, value = line.split("\t")
                c12 = " ".join(coord.strip().strip("()").split(", "))
                f.write(" ".join([c12, value]))
# --- private formats ---
def read_meta(fp):
    """
    Read general metadata, take care of sample_name.
    """
    fp = str(fp)
    if fp.endswith(".csv") or fp.endswith(".csv.gz"):
        df = pd.read_csv(fp,index_col=0,dtype={0:"string"})
    elif fp.endswith(".tsv") or fp.endswith(".tsv.gz"):
        df = pd.read_table(fp,index_col=0,dtype={0:"string"})
    else:
        raise ValueError("Unknown file extension")
    df.index.name = "sample_name"
    df.columns = df.columns.astype("string")
    return df
def read_umi_tools(path, sep=None)->pd.DataFrame:
    """
    Read umi_tools long-form output matrix.
    A wrapper for read_expr.
    By default umi_tools file store matrix in VO format, this function will transpose it to OV.
    See read_expr for more details.
    Input:
        path: path to matrix file
        sep: separator, if None will infer from file extension
    Output:
        obs * var matrix
    """
    return read_expr(path,sep).T
def matra(file):
    """
    Read umi-tools output to anndata object.
    """
    expr = read_umi_tools(file,"\t")
    expr.columns = expr.columns.astype("str")
    expr.columns.name = "gene"
    expr.index = expr.index.astype("str")
    expr.index.name = "sample_name"
    adata = ad.AnnData(expr)
    return adata
def get_chrom_contact_counts(dump_dir):
    """
    Get per-chromosome-per-phase contacts counts of samples.
    Input:
        dump_dir: rd/dump generated by hires_utils
    Output:
        pd.DataFrame
    """
    res = pd.concat((pd.read_pickle(file) for file in [os.path.join(dump_dir, i) for i in os.listdir(dump_dir)]),axis=1)
    res.columns = [i.split("_")[0] for i in os.listdir(dump_dir)]
    return res
def align_to_df(primary_views):
    """
    Extract general dataframe from primary_view function output.
    Input:
        primary_views: complex data structure list of dict
    Output:
        dataframe
    """
    mapper = {
        "head-tail":0,
        "dorsal-ventral":1,
        "left-right":2,
    }
    dfs = [[],[],[]]
    samples = []  
    for sample in primary_views:
        samples.append(sample)
        for aname, figures in zip(
            primary_views[sample]["name_of_vectors"],
            primary_views[sample]["primary_figures"]
        ):
            for figure in figures: # 2 figure for each axis
                figure[figure == np.inf] = np.nan
                dfs[mapper[aname]].append(figure.reshape(1,-1))
    index = list(product(samples, ["A","B"]))
    #return dfs
    dfs = [
        pd.DataFrame(
            np.concatenate(df,axis=0),
            index=pd.MultiIndex.from_product([samples, ["A", "B"]])
        )
        for df in dfs
    ]
    return dfs
# --- misc ---
def dump_json(obj, filep):
    """
    Dump to json shortcut.
    """
    with open(filep,"wt") as f:
        json.dump(obj,f)
def load_json(filep):
    """
    Load from json file shortcut.
    """
    with open(filep,"rt") as f:
        return json.load(f)
def load_pickle(filep):
    """
    Load from pickle file shortcut.
    """
    with open(filep,"rb") as f:
        return pickle.load(f)
def dump_pickle(obj, filep):
    """
    Dump to pickle shortcut.
    """
    with open(filep,"wb") as f:
        pickle.dump(obj,f)
def get_ref_dir():
    """
    Return a path to src's ref
    """
    ref_dir = os.path.join(os.path.dirname(__file__), "ref")
    return os.path.join(ref_dir, "")