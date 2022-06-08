import os
from numpy.testing import assert_allclose
from wet.adtools import gen_adata
def test_gen_adata(request):
    expr = os.path.join(request.fspath.dirname, "data", "expr.csv.gz")
    pass