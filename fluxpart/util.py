from collections import namedtuple
from itertools import permutations

import numpy as np
import pywt

from fluxpart.constants import MOLECULAR_WEIGHT as MW    # kg/mol


def stats2(sarray, names=None):
    """Calculate means and (co)variances for structured array data."""

    if names is None:
        names = sarray.dtype.names
    nvar = len(names)
    data = tuple(sarray[name] for name in names)
    cov = np.cov(data)
    nondiag_cov = list(cov[i, j] for i, j in permutations(range(nvar), 2))

    names_ave = list('ave_' + name for name in names)
    names_var = list('var_' + name for name in names)
    names_cov = list(
        'cov_' + n1 + "_" + n2 for n1, n2 in permutations(names, 2))

    out = dict(zip(names_ave, np.mean(data, axis=1)))
    out.update(zip(names_var, cov.diagonal()))
    out.update(zip(names_cov, nondiag_cov))

    NamedStats = namedtuple('Stats2', names_ave + names_var + names_cov)
    return NamedStats(**out)


def qflux_mass_to_heat(massflux, Tk):                # kg/m^2/s, K
    Lv = 2.5008e6 - 2366.8 * (Tk - 273.15)           # J/kg
    return massflux * Lv                             # W/m^2


def cflux_mass_to_mol(massflux):                     # kg/m^2/s
    return 1. / MW.co2 * massflux                    # mol/m^2/s


def progressive_lowcut_series(series):
    """Progressively remove low-frequency components of 1D series.

    Yields sequence in which the low-frequency (large-scale) components
    of `series` are progressively removed.  The sequence is obtained
    from reconstrution of a multilevel discrete haar wavelet
    decompositon of `series`.

    Parameters
    ----------
    series : array_like
        1D data series with a length that is a power of 2

    Yields
    -------
    lowcut_series : array
        Sequence of progressively lowcut filtered data `series`. The
        yielded series have the same length as `series`.

    Notes
    -----
    After an N level discrete wavelet decomposition, a data series S can
    be reconstructed in terms of wavelet 'approximations' (A) and
    'details' (D):

    S = A(N) + D(N) + D(N-1) ... D(2) + D(1)
      = A(N-1) + D(N-1) + ... D(2) + D(1)        [A(N-1) = A(N) + D(N)]
        ...
      = A(1) + D(1)                              [(A(1) = A(2) + D(2)]

    where A(N) represents the 'lowest' level approximation [e.g., for
    the haar wavelet and a complete decomposition of dyadic length data,
    A(N) is equal to mean(S)]

    The sequence returned by this function is:

    S - A(N),
    S - A(N-1),
    S - A(N-2),
    ...,
    S - A(1)

    This sequence is computed by the equivalent:

    S - A(N),
    S - A(N) - D(N),
    S - A(N) - D(N) - D(N-1),
    ...,
    S - A(N) - D(N) - D(N-1) - ... - D(2),

    i.e. the details are removed in sequence

    lowcut = S - A(N)
    for j = N to 2
        lowcut = lowcut - D(j)

    N.B.
    ----
    The length of `series` should be a power of 2. Otherwise ... I don't
    know ... probably would have to do something to ensure the filtered
    series have the same length as `series`.
    """
    series_data = np.asarray(series)
    wavelet = pywt.Wavelet('haar')
    nlevel = pywt.dwt_max_level(series_data.size, wavelet.dec_len)
    decomp_coef = pywt.wavedec(series_data, wavelet=wavelet, level=nlevel)
    cAn, cD = decomp_coef[0], decomp_coef[1:]
    lowcut_series = series_data - pywt.upcoef('a', cAn, wavelet, level=nlevel)
    yield lowcut_series
    for j, cDj in enumerate(cD[:-1]):
        Dj = pywt.upcoef('d', cDj, wavelet, level=nlevel - j)
        lowcut_series -= Dj
        yield lowcut_series


if __name__ == "__main__":
    pass
