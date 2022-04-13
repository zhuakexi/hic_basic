import six
import re
import numpy as np
import pandas as pd

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