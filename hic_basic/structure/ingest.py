from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd
import xarray as xr
from hires_utils.hires_io import parse_3dg
from tqdm import tqdm

from ..binnify import GenomeIdeograph

def _3dg_to_xr(_3dg_file, sample, genome=None, binsize=20000, flavor="hickit"):
    """
    Convert a .3dg file to xarray dataset.
    Input:
        _3dg_file: .3dg file path
        sample: sample name, used as a dimension name
        genome: genome name, used to determine bins
        binsize: binsize of input .3dg file
        flavor: flavor of bins of .3dg file, see GenomeIdeograph.bins
    Output:
        xarray dataset
    """
    if genome is not None:
        assert binsize is not None
    _3dg = parse_3dg(_3dg_file)
    if genome is not None:
        bins = GenomeIdeograph(genome).bins(
            binsize = binsize,
            bed=True,
            flavor=flavor
        )
        _3dg = pd.merge(
            bins,
            _3dg,
            left_on = ["chrom","start"],
            right_index = True,
            how = "left"
        )
        _3dg = _3dg.set_index(["chrom","start"]).drop(
            "end",
            axis=1
            )
    _3dg.columns.name = "features"
    _3dg = _3dg.stack().sort_index()
    _3dg_xr = xr.DataArray.from_series(
        _3dg
    )
    # add a sample name dimension
    _3dg_xr = _3dg_xr.expand_dims(sample = [sample])
    _3dg_xr_dataset = _3dg_xr.to_dataset(name = "3dg")
    return _3dg_xr_dataset
def _3dg2netcdf(_3dg_file, sample, output, genome="GRCh38", binsize=20000000, flavor="hickit", force=False):
    """
    Convert a .3dg file to xarray dataset and save it to netcdf file.
    Input:
        _3dg_file: .3dg file path
        sample: sample name, used as a dimension name
        output: output netcdf file path
        genome: genome name, used to determine bins
        binsize: binsize of input .3dg file
        flavor: flavor of bins of .3dg file, see GenomeIdeograph.bins
        force: whether to overwrite existing file
    Output:
        output: output netcdf file path
    """
    if Path(output).exists() and not force:
        return output
    _3dg_xr_dataset = _3dg_to_xr(_3dg_file, sample, genome=genome, binsize=binsize)
    _3dg_xr_dataset.to_netcdf(
        output,
        encoding = {
            "sample": {"dtype": "str"}
        }
        )
    return output
def _3dgs2netcdfs(_3dg_files:list,samples:list,outdir:str,
    genome="GRCh38",binsize=20000,flavor="hickit",force=False):
    """
    Convert 3dg files to aligned netcdf files.
    TODO: make it real multithreading
    Input:
        _3dg_files: file paths
        samples: sample name of each _3dg_file
        outdir: where to store nc files, will create if not exists
            nc files are named as "{sample}.nc"
        genome: genome name, used to determine bins
        binsize: binsize of input .3dg file
        flavor: flavor of bins of .3dg file, see GenomeIdeograph.bins
        force: whether to overwrite existing file
    Output:
        results: output netcdf file paths
    """
    outdir = Path(outdir)
    if not outdir.exists():
        outdir.mkdir(parents=True)
    outpat = str(Path(outdir) / "{sample}.nc")
    with ThreadPoolExecutor(1) as executor:
        futures = []
        for sample, _3dg_file in zip(samples, _3dg_files):
            output = outpat.format(sample=sample)
            future = executor.submit(
                _3dg2netcdf,
                _3dg_file,
                sample,
                output,
                genome = genome,
                binsize = binsize,
                flavor = flavor,
                force = force
            )
            futures.append(future)
        results = []
        for future in tqdm(as_completed(futures),desc="samples",total=len(futures)):
            # TODO: fix tqdm 0% print
            results.append(future.result())
    return results