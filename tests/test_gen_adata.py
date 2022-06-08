import os
from numpy.testing import assert_allclose
from ..wet.adtools import gen_adata
from ..io import read_meta
def test_gen_adata(request, tmp_path):
    cache_dir = tmp_path / "test_gen_adata"
    qc = read_meta(os.path.join(request.fspath.dirname, "data", "testqc.csv.gz"))
    expr = os.path.join(request.fspath.dirname, "data", "counts.gene.tsv.gz")
    adata = gen_adata(qc, cache_dir, expr = ([[expr], dict(zip(qc.index.values, qc.index.values))],{"samplelist":qc.index.values}))
    assert adata.obs.shape[0] == 3