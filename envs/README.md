# hic_basic Environment Files

This directory contains conda environment files for the `hic_basic` project.

## Available Environments

### `hic_basic.yaml` (Recommended for most users)

**Minimal, sufficient environment** containing only the packages that are directly imported and used by the `hic_basic` codebase.

- **33 conda packages** + **2 pip packages**
- **93.4% smaller** than the full environment
- All versions tested and match `hic_basic_v096`
- Python 3.10.16

**Installation:**
```bash
micromamba env create -f envs/hic_basic.yaml
micromamba activate hic_basic
pip install -e .
```

**Core dependencies included:**
- Hi-C processing: `cooler`, `cooltools`, `bioframe`
- Data science: `numpy`, `pandas`, `scipy`, `matplotlib`
- Bioinformatics: `pysam`, `pybedtools`, `h5py`
- Single-cell: `anndata`, `scanpy`, `scvelo`
- Machine learning: `scikit-learn`, `scikit-image`, `umap-learn`
- Visualization: `plotly`, `seaborn`, `pymol-open-source`
- Utilities: `dask`, `xarray`, `tqdm`, `click`, `pytest`