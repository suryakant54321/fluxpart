"""
References
----------
[SS08] Scanlon, T.M., and P. Sahu (2008), On the correlation structure
of water vapor and carbon dioxide in the atmospheric surface layer: A
basis for flux partitioning. Water Resour. Res., 44, W10418,
doi:10.1029/2008WR006932.

"""
from types import SimpleNamespace
import numpy.testing as npt

from fluxpart.partition import partition_from_qc_averages

def test_partition_from_qcdat():
    """TODO

    References
    ----------
    [PRV14] L Palatella, G Rana, D Vitale (2014), Towards a flux-
    partitioning procedure based on the direct use of high-frequency
    eddy-covariance data.  Boundary-Layer Meteorol. 153:327--337,
    doi:10.1007/s10546-014-9947-x
    """

    # "physical" April 7 example from [PRV14] Table 1
    qcdat = SimpleNamespace(
        var_q = 0.411163e-3 **2,
        var_c = 5.182580e-6 **2,
        wq = 0.033140e-3,
        wc = -0.472108e-6,
        corr_qc = -0.881017)
    wue = -37.158598e-3
    ns, fc = partition_from_qc_averages(qcdat, wue)
    if ns.success and ns.validroot:
        npt.assert_allclose(ns.var_cp, 5.230500e-6**2, atol=1e-14)
        npt.assert_allclose(ns.corr_cp_cr, -0.757626, atol=1e-2)
        npt.assert_allclose(fc.wqt, 0.012878e-3, atol=1e-7)
        npt.assert_allclose(fc.wqe, 0.020262e-3, atol=1e-7)
        npt.assert_allclose(fc.wcp, -0.476489e-6, atol=1e-9)
        npt.assert_allclose(fc.wcr, 0.004381e-6, atol=1e-11)
    else:
        assert False

    # "non-physical" April 5 example from [PRV14] Table 1
    qcdat = SimpleNamespace(
        var_q = 0.455994e-3**2,
        var_c = 4.544450e-6**2,
        wq = 0.062700e-3,
        wc = -0.712862e-6,
        corr_qc = -0.922292)
    wue = -24.558131e-3
    ns, fc = partition_from_qc_averages(qcdat, wue)
    #NB: valid_pair refers only to var_cp and corr_cp_cr values, not the
    # partitioned flux values, which are not checked in partition
    if ns.success and ns.validroot:
        npt.assert_allclose(ns.var_cp, 3.951510e-6**2, atol=1e-14)
        npt.assert_allclose(ns.corr_cp_cr, -0.724706, atol=1e-2)
        npt.assert_allclose(fc.wqt, 0.025416e-3, atol=1e-7)
        npt.assert_allclose(fc.wqe, 0.037284e-3, atol=1e-7)
        npt.assert_allclose(fc.wcp, -0.624172e-6, atol=1e-9)
        npt.assert_allclose(fc.wcr, -0.088690e-6, atol=1e-10)
    else:
        assert False

    # comparable to Ray Anderson peach data, 2012-06-07 1300
    qcdat = SimpleNamespace(
        var_q = 0.40639e-6,
        var_c = 7.68505e-12,
        wq = 0.1506337e-3,
        wc = -0.6254288e-6,
        corr_qc = -.9501656)
    wue = -7.060177e-3
    ns, fc = partition_from_qc_averages(qcdat, wue)

    if ns.success and ns.validroot:
        npt.assert_allclose(ns.var_cp, 18.92722e-12, atol=1e-1)
        npt.assert_allclose(fc.wqt, 0.145615e-3, atol=1e-6)
        npt.assert_allclose(fc.wqe, 0.005009e-3, atol=1e-6)
        npt.assert_allclose(fc.wcp, -1.028074e-6, atol=1e-7)
        npt.assert_allclose(fc.wcr, 0.402683e-6, atol=1e-8)
    else:
        assert False

    # comparable to Ray Anderson peach data, 2012-06-07 0230
    qcdat = SimpleNamespace(
        var_q=0.001326586e-6,
        var_c=9.948297e-12,
        wq=0.00088955655e-3,
        wc=0.07513186e-6,
        corr_qc=0.8886955)
    wue = -68.77103e-3
    ns, fc = partition_from_qc_averages(qcdat, wue)
    assert ns.success and ns.validroot

if __name__ == '__main__':
    test_partition_from_qcdat()
