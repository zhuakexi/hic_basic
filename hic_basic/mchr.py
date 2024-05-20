import warnings
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
from dask.distributed import LocalCluster, Client
from hires_utils.hires_io import parse_3dg
from tqdm import tqdm

from .binnify import GenomeIdeograph
from .data import chromosomes
from .genome import RegionPair

warnings.filterwarnings('ignore', category=FutureWarning)

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
    _3dg = _3dg.stack(dropna=False, sort=True)
    _3dg_xr = xr.DataArray.from_series(
        _3dg
    )
    # add a sample name dimension
    _3dg_xr = _3dg_xr.expand_dims(sample_name = [sample])
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
            "sample_name": {"dtype": "str"}
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
class Mchr:
    """
    A class to handle multiple .3dg files.
    """
    def __init__(self, _3dg_files:list, samples:list, ddir:str="tmp", ingest=True, genome=None, binsize=20000, flavor="hickit", force=True):
        """
        Initialize a Mchr object.
        Input:
            _3dg_files: file paths
            samples: sample name of each _3dg_file
            ddir: where to store intermediate files
            ingest: whether to convert .3dg files to netcdf files
            genome: genome name, used to determine bins
                e.g. "GRCh38"
            binsize: binsize of input .3dg file
            flavor: flavor of bins of .3dg file, see GenomeIdeograph.bins
        """
        assert genome is not None, "Genome must be specified."
        self._3dg_files = _3dg_files
        self.samples = samples
        self.genome = genome
        self.binsize = binsize
        self.flavor = flavor
        self.ddir = Path(ddir)
        self._3dg_ds_fp = self.ddir / "3dg.nc"
        self.in_disk = False
        if ingest:
            self._ingest_3dgs(self.ddir / "tmp", force=force)
    def _ingest_3dgs(self, outdir:str, force=False):
        """
        Convert 3dg files to aligned netcdf files.
        Input:
            outdir: where to store nc files, will create if not exists
                nc files are named as "{sample}.nc"
            force: whether to overwrite existing file
        Output:
            results: output netcdf file paths
        """
        self._3dg_xr_datasets = _3dgs2netcdfs(
            self._3dg_files,
            self.samples,
            outdir,
            genome = self.genome,
            binsize = self.binsize,
            flavor = self.flavor,
            force = force
        )
        with xr.open_mfdataset(self._3dg_xr_datasets) as ds:
            ds.attrs["genome"] = self.genome
            ds.attrs["binsize"] = self.binsize
            ds.attrs["flavor"] = self.flavor
            ds.attrs["_3dg_files"] = [str(i) for i in self._3dg_files]
            ds.to_netcdf(
                self._3dg_ds_fp,
                encoding = {
                    "sample_name" : {"dtype": "str"},
                },
                )
        self.in_disk = True
    @classmethod
    def from_netcdf(cls, _3dg_ds_fp:str, ddir=None, _3dg_files=None, samples=None, genome=None, binsize=None, flavor=None):
        """
        Initialize a Mchr object from a netcdf file.
        Input:
            _3dg_ds_fp: netcdf file path
        Output:
            Mchr object
        """
        with xr.open_dataset(_3dg_ds_fp) as ds:
            _3dg_files = ds.attrs["_3dg_files"] if _3dg_files is None else _3dg_files
            samples = ds["sample_name"] if samples is None else samples
            genome = ds.attrs["genome"] if genome is None else genome
            binsize = ds.attrs["binsize"] if binsize is None else binsize
            flavor = ds.attrs["flavor"] if flavor is None else flavor
        if ddir is None:
            ddir = Path(_3dg_ds_fp).with_suffix(".tmp")
        mchr = cls(
            _3dg_files,
            samples,
            ddir = ddir,
            genome = genome,
            binsize = binsize,
            flavor = flavor,
            ingest = False
        )
        mchr.in_disk = True
        mchr._3dg_ds_fp = _3dg_ds_fp
        return mchr
    def DM(self, samples:list, region_pair:list, min_samples=None, n_jobs=None)->pd.DataFrame:
        """
        Compute distance matrix.
        Input:
            samples: sample names to compute
            region_pair: genome region to compute DM, can include inter-chromosome regions,
                can cross chromosomes.
                [[(chrom1,pos1),(chrom1,pos2)],[(chrom2,pos1),(chrom2,pos2)]
            proximity: if not None, calculate proximity instead of distance
                count pixels within proximity
            n_jobs: number of threads for parallel computing
                if None, use a simple in-memory method
                else, use dask
        Output:
            DM: distance matrix, index and columns are 2-level multiindex
            (chrom, start)
        """
        return self._calc_DM(
            samples,
            region_pair,
            proximity = None,
            min_samples = min_samples,
            n_jobs = n_jobs
        )
    def PM(self, samples:list, region_pair:list, proximity:int, min_samples=None, n_jobs=None)->pd.DataFrame:
        """
        Compute proximity matrix.
        Input:
            samples: sample names to compute
            region_pair: genome region to compute PM, can include inter-chromosome regions,
                can cross chromosomes.
                [[(chrom1,pos1),(chrom1,pos2)],[(chrom2,pos1),(chrom2,pos2)]
            proximity: count pixels within proximity
            n_jobs: number of threads for parallel computing
                if None, use a simple in-memory method
                else, use dask
        Output:
            PM: proximity matrix, index and columns are 2-level multiindex
            (chrom, start)
        """
        assert proximity is not None, "Proximity must be specified."
        return self._calc_DM(
            samples,
            region_pair,
            proximity = proximity,
            min_samples = min_samples,
            n_jobs = n_jobs
        )
    def _calc_DM(self, samples:list, region_pair:list, proximity=None, min_samples=None, n_jobs=None)->pd.DataFrame:
        """
        Compute distance matrix.
        Input:
            samples: sample names to compute
            region_pair: genome region to compute DM, can include inter-chromosome regions,
                can cross chromosomes.
                [[(chrom1,pos1),(chrom1,pos2)],[(chrom2,pos1),(chrom2,pos2)]
            proximity: if not None, calculate proximity instead of distance
                count pixels within proximity
            n_jobs: number of threads for parallel computing
                if None, use a simple in-memory method
                else, use dask
        Output:
            DM: distance matrix, index and columns are 2-level multiindex
            (chrom, start)
        """
        if min_samples is not None:
            assert len(samples) >= min_samples, "Not enough samples to fulfill min_samples."
        if not self.in_disk:
            print("Ingesting .3dg files to disk...")
            self._3dgs2netcdfs(self.ddir / "tmp", force=True)
        if (proximity is not None) and (min_samples is None):
            # proximity requires at least 1 sample for coverage division
            min_samples = 1
        if n_jobs is None:
            with xr.open_dataset(self._3dg_ds_fp) as ds:
                if samples is not None:
                    ds = ds.sel(sample_name = samples)
                # if "chromint" not in ds.xindexes:
                #     if "chrom" in ds.xindexes:
                #         ds = ds.reset_index("chrom")
                #     if "chromint" not in ds.coords:
                #         chroms = chromosomes(self.genome).index
                #         chrom_dict = dict(zip(chroms,range(len(chroms))))
                #         ds = ds.assign_coords(
                #             chromint = ("chrom", [
                #                 chrom_dict[chrom]
                #                 for chrom in ds.coords["chrom"].to_series().values
                #             ])
                #         )
                #     ds = ds.set_xindex("chromint")
                #     # TODO: make sure no other indexes bound to dim "chrom"
                region1, region2 = RegionPair(
                    region_pair,
                    genome = self.genome,
                    binsize = self.binsize
                    ).region_pair
                ds1 = ds.rename(
                    {
                        "chrom" : "chrom1",
                        "start" : "start1"
                        }
                )
                ds1 = ds1.stack(gpos1=["chrom1","start1"])
                ds2 = ds.rename(
                    {
                        "chrom" : "chrom2",
                        "start" : "start2"
                        }
                )
                ds2 = ds2.stack(gpos2=["chrom2","start2"])
                diff = ds1.rename({"gpos1":"region1"}).sel(region1 = region1.bins) \
                    - ds2.rename({"gpos2":"region2"}).sel(region2 = region2.bins)
                dis_mat = np.sqrt((diff**2).sum(dim="features",skipna=False))
                if proximity is not None:
                    mean_dis_mat = (dis_mat < proximity).sum(dim="sample_name",skipna=True)
                    # mean_dis_mat = mean_dis_mat / cov
                    # divide by coverage after filtering out low coverage bins
                else:
                    mean_dis_mat = dis_mat.mean(dim="sample_name",skipna=True)
                if min_samples is not None:
                    # filter out low coverage bins
                    cov1 = (~np.isnan(ds1)).all(
                        dim = "features"
                    ).sum(
                        dim = "sample_name"
                    )
                    mask1 = (cov1 < min_samples)
                    cov2 = (~np.isnan(ds2)).all(
                        dim = "features"
                    ).sum(
                        dim = "sample_name"
                    )
                    # mask out low coverage bins
                    # where returns cond == True so we need to invert the mask
                    mask2 = (cov2 < min_samples)
                    mean_dis_mat = mean_dis_mat.where(~mask1,other=np.nan).where(~mask2,other=np.nan)
                if proximity is not None:
                    # calculate pixel-level coverage
                    cov1 = (~np.isnan(ds1)).all(
                        dim = "features"
                    )
                    cov2 = (~np.isnan(ds2)).all(
                        dim = "features"
                    )
                    # treat NaN as 0
                    cov = (cov1 & cov2).sum(dim="sample_name",skipna=True)
                    mean_dis_mat = mean_dis_mat / cov
                DM_df = mean_dis_mat.to_dataframe()["3dg"] # a 4-level multiindex series
        else:
            with LocalCluster(n_workers=n_jobs, threads_per_worker=1, memory_limit="2GB") as cluster, Client(cluster) as client:
                print("Computing distance matrix, dashboard available at", client.dashboard_link)
                with xr.open_dataset(
                    self._3dg_ds_fp,
                    chunks={"chrom":1,"start":100}
                    ) as ds:
                    region1, region2 = RegionPair(
                        region_pair,
                        genome = self.genome,
                        binsize = self.binsize
                        ).region_pair
                    if samples is not None:
                        ds = ds.sel(sample_name = samples)
                    # select by chromosome blocks
                    ds1 = ds.rename(
                        {
                            "chrom" : "chrom1",
                            "start" : "start1"
                            }
                    )
                    ds1 = ds1.sel(chrom1=region1.region_chroms)
                    ds2 = ds.rename(
                        {
                            "chrom" : "chrom2",
                            "start" : "start2"
                            }
                    )
                    ds2 = ds2.sel(chrom2=region2.region_chroms)
                    diff = ds1 - ds2
                    dis_mat = np.sqrt((diff**2).sum(dim="features",skipna=False))
                    if proximity is not None:
                        mean_dis_mat = (dis_mat < proximity).sum(dim="sample_name",skipna=True)
                        # divide by coverage after filtering out low coverage bins
                    else:
                        mean_dis_mat = dis_mat.mean(dim="sample_name",skipna=True)
                    if min_samples is not None:
                        # filter out low coverage bins
                        cov1 = (~np.isnan(ds1)).all(
                            dim = "features"
                        ).sum(
                            dim = "sample_name"
                        )
                        mask1 = (cov1 < min_samples)
                        cov2 = (~np.isnan(ds2)).all(
                            dim = "features"
                        ).sum(
                            dim = "sample_name"
                        )
                        mask2 = (cov2 < min_samples)
                        # mask out low coverage bins
                        # where returns cond == True so we need to invert the mask
                        mean_dis_mat = mean_dis_mat.where(~mask1,other=np.nan).where(~mask2,other=np.nan)
                    if proximity is not None:
                        # calculate pixel-level coverage
                        cov1 = (~np.isnan(ds1)).all(
                            dim = "features"
                        )
                        cov2 = (~np.isnan(ds2)).all(
                            dim = "features"
                        )
                        # treat NaN as 0
                        cov = (cov1 & cov2).sum(dim="sample_name",skipna=True)
                        mean_dis_mat = mean_dis_mat / cov
                    DM_df = mean_dis_mat.to_dataframe()["3dg"]
        print("Tidying up distance matrix...")
        DM_df = DM_df.reset_index()
        DM_df = DM_df.astype(
            {
                "chrom1" : pd.CategoricalDtype(
                    GenomeIdeograph("GRCh38").chromosomes.index,
                    ordered=True
                    ),
                "chrom2" : pd.CategoricalDtype(
                    GenomeIdeograph("GRCh38").chromosomes.index,
                    ordered=True
                    )
            }
        )
        DM_mat = DM_df.pivot(
            index = ["chrom1","start1"],
            columns = ["chrom2","start2"],
            values = "3dg"
        )
        # --- filter out non-region bins for dask mode ---
        if n_jobs is not None:
            # only slice chromosomes in dask mode, now continue to slice start positions
            DM_mat = DM_mat.loc[region1.bins]
            DM_mat = DM_mat.loc[:,region2.bins]
        return DM_mat