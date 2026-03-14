import os
import numpy as np
import pandas as pd
import pytest
from numpy.testing import assert_allclose
from hic_basic.interpolate import _gaussian_interpolate, gp_interpolate
import anndata as ad
from hic_basic.hicio import read_meta
def test_interpolate(request):
    """
    Test the interpolate function.
    """
    orig = read_meta(os.path.join(request.fspath.dirname, "data", "origVals.csv.gz"))
    traj = read_meta(os.path.join(request.fspath.dirname, "data", "traj.csv.gz"))
    traj = traj.iloc[:,0]
    ref_res = read_meta(os.path.join(request.fspath.dirname, "data", "interpolatedVals.csv.gz"))
    res, newtraj = _gaussian_interpolate(orig, traj, 0.1, 200)
    assert_allclose(res.values, ref_res.values)
def _make_gp_adata():
    X = np.array([
        [1.0, 2.0, 0.5],
        [2.0, 2.5, 1.0],
        [3.0, 2.0, 1.5],
        [4.0, 1.5, 2.0],
        [5.0, 1.0, 2.5],
    ])
    adata = ad.AnnData(
        X=X,
        obs=pd.DataFrame(
            {"pt": [0.0, 0.25, 0.5, 0.75, 1.0]},
            index=[f"cell_{i}" for i in range(X.shape[0])]
        ),
        var=pd.DataFrame(
            {"use_gp": [True, False, True]},
            index=["gene_a", "gene_b", "gene_c"]
        )
    )
    adata.layers["scaled"] = X * 2
    return adata
def test_gp_interpolate_stores_results_in_uns():
    adata = _make_gp_adata()
    result = gp_interpolate(
        adata,
        "pt",
        gene_mask_key="use_gp",
        numPts=5,
        n_restarts_optimizer=0
    )
    assert result is adata
    assert "gp_interpolate" in adata.uns

    gp_res = adata.uns["gp_interpolate"]
    assert list(gp_res["mean"].index) == ["gene_a", "gene_c"]
    assert gp_res["mean"].shape == (2, 5)
    assert gp_res["variance"].shape == (2, 5)
    assert gp_res["hyperparameters"].shape[0] == 2
    assert gp_res["time_grid"].shape[0] == 5
    assert list(gp_res["mean"].columns) == [f"Gp_{i}" for i in range(5)]
    assert np.all(gp_res["variance"].to_numpy() >= 0)
    assert gp_res["metadata"]["num_cells_used"] == 5
    assert gp_res["metadata"]["num_cells_dropped"] == 0
    assert gp_res["metadata"]["gene_mask_key"] == "use_gp"
def test_gp_interpolate_inplace_false_returns_copy():
    adata = _make_gp_adata()
    result = gp_interpolate(
        adata,
        "pt",
        result_key="gp_test",
        numPts=4,
        inplace=False,
        n_restarts_optimizer=0
    )
    assert result is not adata
    assert "gp_test" not in adata.uns
    assert "gp_test" in result.uns
def test_gp_interpolate_drops_na_pseudotime_and_uses_layer():
    adata = _make_gp_adata()
    adata.obs.loc["cell_4", "pt"] = np.nan

    gp_interpolate(
        adata,
        "pt",
        layer="scaled",
        gene_mask_key="use_gp",
        result_key="gp_scaled",
        numPts=3,
        n_restarts_optimizer=0
    )
    gp_res = adata.uns["gp_scaled"]
    assert gp_res["metadata"]["layer"] == "scaled"
    assert gp_res["metadata"]["num_cells_used"] == 4
    assert gp_res["metadata"]["num_cells_dropped"] == 1
def test_gp_interpolate_validates_inputs():
    adata = _make_gp_adata()
    with pytest.raises(KeyError, match="missing_pt"):
        gp_interpolate(adata, "missing_pt", n_restarts_optimizer=0)

    adata.var["bad_mask"] = [1, 0, 1]
    with pytest.raises(ValueError, match="boolean"):
        gp_interpolate(adata, "pt", gene_mask_key="bad_mask", n_restarts_optimizer=0)

    adata.obs["flat_pt"] = 1.0
    with pytest.raises(ValueError, match="unique"):
        gp_interpolate(adata, "flat_pt", n_restarts_optimizer=0)
