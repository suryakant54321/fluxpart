import os, io
import numpy as np
import numpy.testing as npt
from fluxpart.hfdata import HFData
from fluxpart.fluxpart import converter_func

def test_hfdata_from_txt():

    testdir = os.path.dirname(os.path.realpath(__file__))
    cols = (2,3,4,6,5,7,8)

    fname = os.path.join(testdir, './data/example_hfdata.dat')

    ex = HFData(
        fname,
        cols,
        converters = {
            'T': converter_func(1, 273.15),
            'q': converter_func(1e-3, 0),
            'c': converter_func(1e-6, 0),
            'P': converter_func(1e3, 0)},
        flags=(9,0),
        delimiter=",",
        skip_header=4)

    npt.assert_allclose(ex['u'][0], 2.50825)
    npt.assert_allclose(ex['v'][0], 1.352)
    npt.assert_allclose(ex['w'][0], -0.49575)
    npt.assert_allclose(ex['c'][0], 675.2817e-6)
    npt.assert_allclose(ex['q'][0], 8.287507e-3)
    npt.assert_allclose(ex['T'][0], 26.58066 + 273.15)
    npt.assert_allclose(ex['P'][0], 100.2458e3)
    npt.assert_allclose(ex['u'][-1], 0.38675)
    npt.assert_allclose(ex['v'][-1], -0.8665)
    npt.assert_allclose(ex['w'][-1], -0.9030001)
    npt.assert_allclose(ex['c'][-1], 669.1664e-6)
    npt.assert_allclose(ex['q'][-1], 9.393601e-3)
    npt.assert_allclose(ex['T'][-1], 27.66293 + 273.15)
    npt.assert_allclose(ex['P'][-1], 100.2458e3)

    npt.assert_allclose(ex['u'].mean(), 1.25613, atol=1e-4)
    npt.assert_allclose(ex['v'].mean(), -0.836322, atol=1e-4)
    npt.assert_allclose(ex['w'].mean(), 0.08104079, atol=1e-4)
    npt.assert_allclose(ex['c'].mean(), 669.3725e-6, atol=1e-10)
    npt.assert_allclose(ex['q'].mean(), 9.2696019e-3, atol=1e-7)
    npt.assert_allclose(ex['T'].mean(), 27.418948 + 273.15, atol=1e-4)
    npt.assert_allclose(ex['P'].mean(), 100.23959e3, atol=1e-1)

    toy_data = (
        'foobar baz\n'
        'asdf,0,2,3,4,5,6,7,9,0\n'
        'asdf,1,2,3,4,5,6,7,9,0\n'
        'asdf,2,2,3,4,5,6,7,9,1\n'
        'asdf,3,2,3,4,5,6,,9,0\n'
        'asdf,4,2,3,4,5,6,7,9,0\n'
        '# foo\n'
        'asdf,5,2,3,4,5,6,7,9,0\n'
        'asdf,6,2,3,4,5,6,7,xxx,0\n'
        'asdf,7,???,3,4,5,6,7,9,0\n'
        'asdf,8,2,3,4,5,6,7,9,0\n' 
        'asdf,9,2,3,4,5,6,7,9,0\n' 
        'asdf,10, 2,3,4,5,6,7,9,0\n'
        'asdf,11,-2,3,4,5,6,7,9,0\n'
        )

    ex = HFData(
        io.BytesIO(toy_data.encode()),
        cols = (1,2,3,7,6,4,5),
        comments = '#',
        skip_header = 1,
        missing_values="???",
        converters = {'q' : converter_func(10.,0)},
        flags = (9, 0),
        delimiter = ",",  # blank space already
        rd_tol = 0.1,
        ad_tol = 2,
        bounds = {'v' : (0, np.inf)}
        )

    npt.assert_allclose(ex.data_table['u'], [4,5,6])
    npt.assert_allclose(ex.data_table['v'], 3*[2,])
    npt.assert_allclose(ex.data_table['w'], 3*[3,])
    npt.assert_allclose(ex.data_table['q'], 3*[70,])
    npt.assert_allclose(ex.data_table['c'], 3*[6,])
    npt.assert_allclose(ex.data_table['T'], 3*[4,])
    npt.assert_allclose(ex.data_table['P'], 3*[5,])

if __name__ == '__main__':
    test_hfdata_from_txt()
