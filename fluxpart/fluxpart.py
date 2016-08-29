import numpy as np

import fluxpart.partition as fp
import fluxpart.wue as wue
import fluxpart.util as util
from fluxpart.hfdata import HFData
from fluxpart.containers import SiteData, Fluxes, WUE, Result


def flux_partition(fname, cols, sitedata, unit_convert=None, temper_unit='K',
                   bounds=None, flags=None, rd_tol=0.4, ad_tol=1024,
                   correcting_external=True, adjusting_fluxes=False,
                   meas_wue=None, ci_mod='const_ratio', label=None, **kwargs):
    """Partition CO2, H2O fluxes into stomatal & nonstomatal components.

    Read high frequency eddy covarince data, analyze water vapor and
    carbon dioxide fluxes, and partition the fluxes into stomatal
    (tranpiration, photosynthesis) and nonstomatal (evaporation,
    respiration) components using the method of [SS08]_.

    The following notation is used in variable naming and documentation
    to represent meteorological quantities:

        u, v, w = wind velocities
        q = water vapor mass concentration
        c = carbon dioxide mass concentration
        T = air temperature
        P = total air pressure

    Parameters
    ----------
    fname : str
        Name of delimited file containing high-frequency eddy covariance
        time series data.
    cols : 7*(int,)
        7-tuple of integers indicating the column numbers of `fname`
        containing series data for (u, v, w, q, c, T, P), in that
        order. Uses 0-based indexing.
    sitedata : dict or object
        Container holding field site meta data with the following keys
        (if dict) or attributes (if object):
            meas_ht : float, m
                Eddy covariance measurement height
            canopy_ht : float, m
                Plant canopy height
            ppath : {'C3' or 'C4'}
                Photosynthetic pathway
        Is required unless a value for `meas_wue` is provided, in which
        case `sitedata` can be set to None.
    unit_convert : dict, optional
        Dictionary of multiplication factors required to convert any
        u, v, w, q, c, or P data not in SI units to SI units (m/s,
        kg/m^3, Pa). (Note T is not in that list). The dictionary keys
        are the variable names. For example, if all data in `fname` are
        in SI units except P and c, which are in units of kPa and
        mg/m^3, respectively, then set:
            unit_convert = {'P': 1e3, 'c': 1e-6}
        since it is necessary to multiply the kPa pressure data by 1e3
        to obtain the SI pressure unit (Pa), and the mg/m^3 CO2 data by
        1e-6 to obtain the SI concentration unit (kg/m^3). Is optional
        only if all data are in SI units.
    temper_unit : {'K' or 'C'}, optional
        The units of the temperature data T in `fname`. Default is 'K'.
    bounds : dict, optional
        Dictionary specifying any prescribed lower and upper bounds for
        legal data. Dict entries have the form
            varname: (float, float)
        where varname is one of 'u', 'v', 'w', 'q', 'c', 'T', or 'P',
        and the 2-tuple holds values for the lower and upper bounds:
        (lower, upper).  Data records are rejected if a variable in the
        record is outside the prescribed bounds. Default is
            bounds = {'c': (0, np.inf), 'q': (0, np.inf)}
        such that data records are rejected if c or q data are not
        positive values.
    flags : 2-tuple or list of 2-tuples, optional
        Specifies that one or more columns in `fname` are used to flag
        bad data records. Each tuple is of the form (col, badval),
        where col is an int specifying the column number containing the
        flag (0-based indexing), and badval is the value of the flag
        that indicates a bad data record.
    rd_tol : float, optional
    ad_tol : int, optional
        Tolerances for rejecting the datafile. If the datafile contains
        bad records (not readable, out-of-bounds, or flagged data), the
        partitioning analysis is performed using the longest stretch of
        consecutive good data records found, unless that stretch is too
        short, in which case the analysis is aborted. The criteria for
        judging 'too short' can be specified in both relative and
        absolute terms: the datafile is rejected if the good stretch is
        a fraction of the total data that is less than `rd_tol`, and/or
        is less than `ad_tol` records long.  Defaults are `rd_tol` = 0.4
        and `ad_tol` = 1024.
    correcting_external : bool, optional
        If True (default), the water vapor and carbon dioxide series
        data are corrected for external fluctuations associated with air
        temperature and vapor density.
    adjusting_fluxes, bool, optional
        If True, the final partitioned fluxes are adjusted to ensure the
        total fluxes match exactly the fluxes indcated in the original
        data. Default is False.
    meas_wue : float, optional
        Measured leaf-level water use efficiency (kg CO2 / kg H2O). Due
        to sign conventions (photosynthesis CO2 flux < 0, transpiration
        vapor flux > 0), water use efficiency must be `meas_wue` < 0.
    ci_mod : str or [str, float] or [str, (float, float)], optional
        Method used to estimate the leaf intercellular CO2 concentration
        when calculating water use efficiency.
        :ref:`(See here) <fluxpart.wue.py>' for a full description
        of `ci_mod`.  Default is `ci_mod` = 'const_ratio'.
        Not used if meas_wue if provided.
    label : object, optional
        Identification label for the data set. Could be a str, int,
        datetime object, etc.
    kwargs
        Keyword arguments passed to numpy.genfromtxt to specify
        formatting of the delimited datafile. See numpy documentation
        for full description of available options. Among the most
        commonly required are:
    delimiter : str, int, or sequence, optional
        The string used to separate values. By default, any consecutive
        whitespaces act as delimiter. An integer or sequence of integers
        can also be provided as width(s) of each field.
    skip_header : int, optional
        The number of lines to skip at the beginning of the file.
    skip_footer : int, optional
        The number of lines to skip at the end of the file.
    comments : str, optional
        The character used to indicate the start of a comment. All the
        characters occurring on a line after a comment are discarded

    Returns
    -------
    dict, with keys:
        'label': `label`
        'result': Result namedtuple
        'fluxes': Fluxes namedtuple
        'datsumm': HFSummary namedtuple, or None
        'wue': WUE namedtuple, or None
        'numsoln': NumerSoln namedtuple, or None

    Examples
    --------
    see TODO

    """

    usecols = np.array(cols, dtype=int).reshape(7,)

    converters = None
    if unit_convert:
        converters = {
            k: converter_func(float(v), 0.) for k, v in unit_convert.items()}
    if temper_unit.upper() == 'C' or temper_unit.upper() == 'CELSIUS':
        converters = converters or {}
        converters['T'] = converter_func(1., 273.15)

    # read high frequency data
    try:
        hfdat = HFData(fname, cols=usecols, converters=converters, flags=flags,
                       bounds=bounds, rd_tol=rd_tol, ad_tol=ad_tol, **kwargs)
    except (TypeError, ValueError) as err:
        mssg = 'HFData read fail: ' + err.args[0]
        result = Result(dataread=False, valid_partition=False, mssg=mssg)
        return {'label': label,
                'result': result,
                'fluxes': Fluxes(*np.full(12, np.nan)),
                'datsumm': None,
                'wue': None,
                'numsoln': None}

    # preliminary data processing and analysis
    hfdat.truncate()
    if correcting_external:
        hfdat.qc_correct()
    hfsum = hfdat.summarize()

    # exit if water vapor flux is negative (directed downward)
    if hfsum.Fq <= 0:
        mssg = ('Negative (downward) water vapor flux, Fq = {:.4}, is '
                'incompatible with partitioning algorithm'.format(hfsum.Fq))
        result = Result(dataread=True, valid_partition=False, mssg=mssg)
        return {'label': label,
                'result': result,
                'fluxes': Fluxes(*np.full(12, np.nan)),
                'datsumm': hfsum,
                'wue': None,
                'numsoln': None}

    # get or calculate water use efficiency
    try:
        lwue = WUE(float(meas_wue), *8 * (np.nan,))
    except TypeError:
        try:
            meta = set_fieldsite_data(sitedata)
        except (TypeError, KeyError, AttributeError) as err:
            lwue = WUE(*np.full(9, np.nan))
            mssg = err.args[0]
            result = Result(dataread=True, valid_partition=False, mssg=mssg)
        else:
            lwue = wue.water_use_efficiency(hfsum, meta, ci_mod)
            if lwue.wue is np.nan or lwue.wue >= 0:
                mssg = ('wue={} must be less than zero'.format(lwue.wue))
                result = Result(dataread=True, valid_partition=False,
                                mssg=mssg)

    # exit if wue value is no good
    if lwue.wue is np.nan or lwue.wue >= 0:
        return {'label': label,
                'result': result,
                'fluxes': Fluxes(*np.full(12, np.nan)),
                'datsumm': hfsum,
                'wue': lwue,
                'numsoln': None}

    # compute partitioned fluxes and return
    pout = fp.partition_from_wqc_series(hfdat['w'], hfdat['q'], hfdat['c'],
                                        lwue.wue, adjusting_fluxes)

    result = Result(dataread=True,
                    valid_partition=pout['valid_partition'],
                    mssg=pout['partmssg'])

    if pout['valid_partition']:
        fluxes = Fluxes(
            *pout['fluxcomps'],
            LE=util.qflux_mass_to_heat(pout['fluxcomps'].wq, hfsum.T),
            LEt=util.qflux_mass_to_heat(pout['fluxcomps'].wqt, hfsum.T),
            LEe=util.qflux_mass_to_heat(pout['fluxcomps'].wqe, hfsum.T),
            Fc_mol=util.cflux_mass_to_mol(pout['fluxcomps'].wc),
            Fcp_mol=util.cflux_mass_to_mol(pout['fluxcomps'].wcp),
            Fcr_mol=util.cflux_mass_to_mol(pout['fluxcomps'].wcr))
    else:
        fluxes = Fluxes(*np.full(12, np.nan))

    return {'label': label,
            'result': result,
            'fluxes': fluxes,
            'datsumm': hfsum,
            'wue': lwue,
            'numsoln': pout['numsoln'],
            'qcdat': pout['qcdat']}


def converter_func(slope, intercept):
    """Return a function to be used to transform data when reading file.
    """
    def func(stringval):
        try:
            return slope * float(stringval.strip()) + intercept
        except ValueError:
            return np.nan
    return func


def set_fieldsite_data(data):
    """Parse and verify site `data`, return a SiteData namedtuple."""

    # parse `data`
    try:
        meas_ht = data['meas_ht']
        canopy_ht = data['canopy_ht']
        ppath = data['ppath']
    except KeyError as err:
        err.args = ('sitedata dict missing required entry: ' + err.args[0],)
        raise
    # not a dict, check attributes
    except TypeError:
        try:
            meas_ht = data.meas_ht
            canopy_ht = data.canopy_ht
            ppath = data.ppath
        # lacks required attribures, try list
        except AttributeError:
            try:
                meas_ht = data[0]
                canopy_ht = data[1]
                ppath = data[2]
            # unable to parse `data`, raise error
            except IndexError as err:
                raise TypeError('Unable to parse sitedata = {}'.format(data))

    # `data` parse is success, verify values
    try:
        ppath = str(ppath).upper()
        if ppath not in ('C3', 'C4'):
            raise ValueError
    except (TypeError, ValueError) as err:
        err.args = ("sitedata(ppath={}): ppath must be 'C3' or 'C4'"
                    "".format(ppath),)
        raise
    try:
        meas_ht = float(meas_ht)
        canopy_ht = float(canopy_ht)
    except ValueError as err:
        err.args = (("sitedata(meas_ht={}, canopy_ht={}): "
                     + err.args[0]).format(meas_ht, canopy_ht),)
        raise
    if meas_ht <= 0 or canopy_ht <= 0:
        raise ValueError("sitedata(meas_ht={}m, canopy_ht={}m): both values "
                         "must be positive".format(meas_ht, canopy_ht),)

    # data values are plausible, return SiteData tuple
    return SiteData(meas_ht=meas_ht, canopy_ht=canopy_ht, ppath=ppath)


if __name__ == "__main__":
    pass
