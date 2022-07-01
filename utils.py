import six
import re
import os
import numpy as np
import pandas as pd

from .io import parse_gtf
def two_sets(ref, check, warning=False):
    """
    Check two sample list.
    Input:
        ref: reference list, iterable.
        check: list to be checked, iterable.
    Output:
        [missing, extra]; list of sets
    """
    missing_sample = set(ref) - set(check)
    extra = set(check) - set(ref)
    if warning:
        if len(missing_sample):
            print("[two_sets] Warning: can't find value for %d samples:" % len(missing_sample))
            print(",".join(list(missing_sample)[:20]))
            if len(extra) > 20:
                print("...")
        if len(extra):
            print("[two_sets] Warning: find extra value for:")
            print(",".join(list(extra)[:20]))
            if len(extra) > 20:
                print("...")
    return missing_sample, extra
# --- util functions paste from cooler ---
def natsort_key(s, _NS_REGEX=re.compile(r"(\d+)", re.U)):
    return tuple([int(x) if x.isdigit() else x for x in _NS_REGEX.split(s) if x])
def argnatsort(array):
    array = np.asarray(array)
    if not len(array):
        return np.array([], dtype=int)
    cols = tuple(zip(*(natsort_key(x) for x in array)))
    return np.lexsort(cols[::-1])
def read_chromsizes(
    filepath_or,
    name_patterns=(r"^chr[0-9]+$", r"^chr[XY]$", r"^chrM$"),
    all_names=False,
    **kwargs
):
    """
    Parse a ``<db>.chrom.sizes`` or ``<db>.chromInfo.txt`` file from the UCSC
    database, where ``db`` is a genome assembly name.

    Parameters
    ----------
    filepath_or : str or file-like
        Path or url to text file, or buffer.
    name_patterns : sequence, optional
        Sequence of regular expressions to capture desired sequence names.
        Each corresponding set of records will be sorted in natural order.
    all_names : bool, optional
        Whether to return all contigs listed in the file. Default is
        ``False``.

    Returns
    -------
    :py:class:`pandas.Series`
        Series of integer bp lengths indexed by sequence name.

    References
    ----------
    * `UCSC assembly terminology <http://genome.ucsc.edu/FAQ/FAQdownloads.html#download9>`_
    * `GRC assembly terminology <https://www.ncbi.nlm.nih.gov/grc/help/definitions>`_

    """
    if isinstance(filepath_or, six.string_types) and filepath_or.endswith(".gz"):
        kwargs.setdefault("compression", "gzip")
    chromtable = pd.read_csv(
        filepath_or,
        sep="\t",
        usecols=[0, 1],
        names=["name", "length"],
        dtype={"name": str},
        **kwargs
    )
    if not all_names:
        parts = []
        for pattern in name_patterns:
            part = chromtable[chromtable["name"].str.contains(pattern)]
            part = part.iloc[argnatsort(part["name"])]
            parts.append(part)
        chromtable = pd.concat(parts, axis=0)
    chromtable.index = chromtable["name"].values
    return chromtable["length"]
def binnify(chromsizes, binsize):
    """
    Divide a genome into evenly sized bins.

    Parameters
    ----------
    chromsizes : Series
        pandas Series indexed by chromosome name with chromosome lengths in bp.
    binsize : int
        size of bins in bp

    Returns
    -------
    bins : :py:class:`pandas.DataFrame`
        Dataframe with columns: ``chrom``, ``start``, ``end``.

    """

    def _each(chrom):
        clen = chromsizes[chrom]
        n_bins = int(np.ceil(clen / binsize))
        binedges = np.arange(0, (n_bins + 1)) * binsize
        binedges[-1] = clen
        return pd.DataFrame(
            {"chrom": [chrom] * n_bins, "start": binedges[:-1], "end": binedges[1:]},
            columns=["chrom", "start", "end"],
        )

    bintable = pd.concat(map(_each, chromsizes.keys()), axis=0, ignore_index=True)

    bintable["chrom"] = pd.Categorical(
        bintable["chrom"], categories=list(chromsizes.index), ordered=True
    )

    return bintable
def gen_fileis(sample_table, dir_path, str_pat):
    """
    Check and generate fileis(input file tables storing file path) for downstream calculating.
    Input:
        sample_table: meta, sample_names as index, try to ensure all samples have their input file.
        dir_path: where to find those input files.
        str_pat: string pattern, gen exact file name. using {} to represent sample_name.
    Output:
        pd.DataFrame
    """
    if not sample_table.index.is_unique:
        raise ValueError("[gen_fileis] Error: sample_table index is not unique")
    a = dict(zip([str_pat.format(i) for i in sample_table.index], list(sample_table.index)))
    b = set(os.listdir(dir_path))
    hitting = a.keys() & b
    missing_sample = [a[k] for k in (a.keys() - hitting)]
    extra = list(b - hitting)
    
    if len(missing_sample):
        print("[gen_fileis] Warning: can't finds files for %d samples:" % len(missing_sample))
        print(",".join(missing_sample))
    if len(extra):
        print("[gen_fileis] Warning: find extra files for:")
        print(",".join(extra[:20]))
        if len(extra) > 20:
            print("...")
    return pd.Series([os.path.join(dir_path, k) for k in hitting], index=[a[k] for k in hitting])
def _mouse_id_name_table(ref_file):
    """
    Get a table of mouse IDs and names from a reference file.
    Input:
        ref_file
    Output:
        pd.DataFrame
    """
    ref = parse_gtf(ref_file, True, True)
    ref = ref[["gene_id", "gene_name"]].drop_duplicates(ignore_index = True)
    gene_id = ref["gene_id"].str.split(".", expand=True)
    gene_id.columns = ["main_gene_id", "gene_id_surplus"]
    id_name_table = pd.concat([gene_id["main_gene_id"], ref["gene_name"]], axis=1, join="inner")
    id_name_table = id_name_table.set_index("main_gene_id")
    if not id_name_table.index.is_unique:
        print("Warning id_name_table index is not unique")
    return id_name_table
def mouse_id2name(id_lists, ref_file, multi=False):
    """
    Convert moust gene_IDs to gene_names.
    Input:
        id_list: gene_ID list or list of gene_ID lists, do this because ref parsing is slow
            e.g. [["ENSG0000012345", "ENSG0000012346"], ["ENSG0000012347"]]
        ref_file: reference file
        multi: whether id_lists has multiple lists
    Output:
        list of gene_names
    """
    id_name_table = _mouse_id_name_table(ref_file)["gene_name"].to_dict()
    if multi:
        name_lists = []
        for id_list in id_lists:
            two_sets(id_list, id_name_table.keys(), warning=True)
            name_lists.append(np.array([id_name_table.get(id) for id in id_list]))
        return name_lists
    else:
        two_sets(id_list, id_name_table.keys(), warning=True)
        return np.array([id_name_table.get(id) for id in id_list])