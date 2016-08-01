import os
from types import SimpleNamespace
import numpy.testing as npt
from fluxpart import flux_partition 
from fluxpart.fluxpart import set_fieldsite_data

def test_flux_partition():

    testdir = os.path.dirname(os.path.realpath(__file__))
    fname = os.path.join(testdir, 
            'data/18June2012/TOA5_6843.ts_Above_2012_06_07_1300.dat')
    cols = (2,3,4,6,5,7,8)

    site = SimpleNamespace(meas_ht=7.11, canopy_ht=4.42, ppath='C3')
    
    ans = flux_partition(
        fname,
        cols=cols,
        sitedata=site,
        delimiter=",",
        skip_header=4,
        unit_convert={
            'q': 1e-3,
            'c': 1e-6,
            'P': 1e3},
        ci_mod='const_ppm',
        temper_unit='C')

    npt.assert_allclose(ans['numsoln'].var_cp, 18.9272e-12, atol=1e-12)

    wcp, wcr, wqe, wqt = (-1.02807378762793e-6, 0.402683101177712e-6, 
                          0.00500869036088203e-3, 0.145615860044424e-3)
    fc = ans['fluxes']
    npt.assert_allclose(fc.Fcp, wcp, atol=0.1e-6)
    npt.assert_allclose(fc.Fcr, wcr, atol=0.1e-6)
    npt.assert_allclose(fc.Fqe, wqe, atol=0.01e-3)
    npt.assert_allclose(fc.Fqt, wqt, atol=0.01e-3)

    (Fcp,Fcr,LEe,LEt) = (-23.3600042633021e-6, 9.14980916104776e-6,
                         12.1870591763138, 354.310004313924)
    afc = ans['fluxes']
    npt.assert_allclose(afc.Fcp_mol, Fcp, rtol=0.15)
    npt.assert_allclose(afc.Fcr_mol, Fcr, rtol=0.15)
    npt.assert_allclose(afc.LEe, LEe, rtol=0.15)
    npt.assert_allclose(afc.LEt, LEt, rtol=0.15)

    fname = os.path.join(testdir, 
            'data/18June2012/TOA5_6843.ts_Above_2012_06_07_1245.dat')
    ans = flux_partition(
        fname,
        cols=cols,
        sitedata=site,
        delimiter=",",
        skip_header=4,
        unit_convert={
            'q': 1e-3,
            'c': 1e-6,
            'P': 1e3},
        temper_unit='C',
        ci_mod='const_ppm')

def test_set_fieldsite_data():
    sd = set_fieldsite_data({'meas_ht': 3, 'canopy_ht': 2.5, 'ppath': 'C4'})
    assert_sd(sd)
    sd = set_fieldsite_data([3, 2.5, 'C4'])
    assert_sd(sd)
    sd = set_fieldsite_data((3, 2.5, 'C4'))
    assert_sd(sd)
    sd = set_fieldsite_data(
            SimpleNamespace(meas_ht=3, canopy_ht=2.5, ppath='C4'))
    assert_sd(sd)

def assert_sd(sd):
    npt.assert_allclose(sd.meas_ht, 3)
    npt.assert_allclose(sd.canopy_ht, 2.5)
    assert(sd.ppath == 'C4')
    

if __name__ == "__main__":
    test_flux_partition()
