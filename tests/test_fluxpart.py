import os
from types import SimpleNamespace
import numpy.testing as npt
from fluxpart import flux_partition
#from fluxpart.fluxpart import _set_fieldsite_data

TESTDIR = os.path.dirname(os.path.realpath(__file__))


def test_flux_partition():
    """Integrated test of highest level flux_partition function.

    In test below, the "desired" results are not exact or even
    necessarily correct. They were obtained with an independent (matlab)
    code, and there should be reasonable agreement.
    """

    cols = (2, 3, 4, 6, 5, 7, 8)
    wue_data = {'meas_ht':7.11, 'canopy_ht':4.42, 'ppath': 'C3',
                 'ci_mod': 'const_ppm'}

    # soln exists for this data without any wavelet filtering
    fname = os.path.join(TESTDIR,
                         'data/TOA5_6843.ts_Above_2012_06_07_1300.dat')

    matlab_fluxes = SimpleNamespace(
        Fcp=-1.02807378762793e-6,
        Fcr=0.402683101177712e-6,
        Fqe=0.00500869036088203e-3,
        Fqt=0.145615860044424e-3,
        Fcp_mol=-23.3600042633021e-6,
        Fcr_mol=9.14980916104776e-6,
        LEe=12.1870591763138,
        LEt=354.310004313924)

    result = flux_partition(
        fname,
        cols=cols,
        wue_params=wue_data,
        delimiter=",",
        skip_header=4,
        unit_convert={
            'q': 1e-3,
            'c': 1e-6,
            'P': 1e3},
        temper_unit='C')

    npt.assert_allclose(result['numsoln'].var_cp, 18.9272e-12, atol=1e-12)
    assert_flux_components(result['fluxes'], matlab_fluxes)

    # soln is obtained after some wavelet filtering
    fname = os.path.join(TESTDIR,
                         'data/TOA5_6843.ts_Above_2012_06_07_1245.dat')

    matlab_fluxes = SimpleNamespace(
        Fcp=-0.866856083109642e-6,
        Fcr=0.353428894620522e-6,
        Fqe=0.0124697200158396e-3,
        Fqt=0.117438136138301e-3,
        Fcp_mol=-23.1074540435422e-6,
        Fcr_mol=10.6590633820467e-6,
        LEe=35.6007693518818,
        LEt=335.283229492226)

    result = flux_partition(
        fname,
        cols=cols,
        wue_params=wue_data,
        delimiter=",",
        skip_header=4,
        unit_convert={
            'q': 1e-3,
            'c': 1e-6,
            'P': 1e3},
        temper_unit='C')

    npt.assert_allclose(result['numsoln'].var_cp, 15.2944e-12, atol=1e-12)
    assert_flux_components(result['fluxes'], matlab_fluxes)


def assert_flux_components(calc, desired):
    npt.assert_allclose(calc.Fcp, desired.Fcp, atol=0.1e-6)
    npt.assert_allclose(calc.Fcr, desired.Fcr, atol=0.1e-6)
    npt.assert_allclose(calc.Fqe, desired.Fqe, atol=0.01e-3)
    npt.assert_allclose(calc.Fqt, desired.Fqt, atol=0.01e-3)
    npt.assert_allclose(calc.Fcp_mol, desired.Fcp_mol, atol=10e-6)
    npt.assert_allclose(calc.Fcr_mol, desired.Fcr_mol, atol=10e-6)
    npt.assert_allclose(calc.LEe, desired.LEe, atol=10)
    npt.assert_allclose(calc.LEt, desired.LEt, atol=50)


if __name__ == "__main__":
    test_flux_partition()
