import os
from numpy.testing import assert_allclose
from hic_basic.interpolate import _gaussian_interpolate
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