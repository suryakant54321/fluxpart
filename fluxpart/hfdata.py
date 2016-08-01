import numpy as np

import fluxpart.util as util
from fluxpart.containers import HFSummary
from fluxpart.constants import MOLECULAR_WEIGHT as MW
from fluxpart.constants import SPECIFIC_HEAT_CAPACITY as CP
from fluxpart.constants import SPECIFIC_GAS_CONSTANT as GC


class HFData:
    """High-frequency eddy covariance data.

    The following notation is used in variable naming and documentation
    to represent meteorological quantities (SI units):

        u, v, w = wind velocities (m/s)
        q = water vapor mass concentration (kg/m^3)
        c = carbon dioxide mass concentration (kg/m^3)
        T = air temperature (K)
        P = total air pressure (Pa)

    Attributes
    ----------
    data_table : structured array
        Table of high frequency data series with columns 'u', 'v', 'w',
        'q', 'c', 'T', and 'P'
    """

    def __init__(self, fname, cols, converters=None, flags=None, bounds=None,
                 rd_tol=0.5, ad_tol=1024, **kwargs):
        """Read high frequency eddy covarince data, apply some QA/QC.

        Parameters
        ----------
        fname : str
            Name of delimited file containing high-frequency eddy covariance
            time series data.
        cols : 7*(int,)
            7-tuple of integers indicating the column numbers of `fname`
            containing series data for (u, v, w, q, c, T, P), in that
            order. Uses 0-based indexing.
        converters : dict, optional
            Dictionary of functions used to convert any non-SI data to
            SI units.  Dict keys are 'u', 'v', 'w', 'q', 'c', 'T', or
            'P', and the funcs take string single argument, e.g.
            ``converters = {'P': lambda s: 1e3*float(s)}``.
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
        """

        var_names = ['u', 'v', 'w', 'q', 'c', 'T', 'P']
        self.names = var_names.copy()
        all_names = var_names.copy()
        self._qc_already_corrected = False
        usecols = np.array(cols, dtype=int).reshape(7,)

        # make sure no duplicate kws
        disallowed_kws = ('usecols', 'names', 'dtype', 'converters',
                          'filling_values')
        for kw in disallowed_kws:
            delete_if_present = kwargs.pop(kw, False)

        var_dtype = list(zip(var_names, 7 * (float, )))
        all_dtype = var_dtype.copy()

        # add flag columns to data table
        if flags:
            if isinstance(flags, list):
                fcols, goodvals = list(zip(*flags))
            else:
                fcols, goodvals = [flags[0], ], [flags[1], ]
            flag_names = list('flg' + str(i) for i in range(len(fcols)))

            usecols = np.append(usecols, np.array(fcols, dtype=int))
            all_dtype += list(zip(flag_names, [type(v) for v in goodvals]))
            all_names += flag_names

            # for later convenience
            flag_dict = dict(zip(flag_names, goodvals))

        data = np.genfromtxt(fname, usecols=usecols, dtype=all_dtype,
                             names=all_names, converters=converters, **kwargs)

        # 1D bool mask is True if any `data` field is missing (=nan)
        mask = np.isnan(data[var_names].view(float).reshape(-1, 7)).any(axis=1)

        # Also mask records with out-of-bounds or flagged data
        dbounds = {'c': (0, np.inf), 'q': (0, np.inf)}
        if bounds:
            dbounds.update(bounds)
        for var, (low, high) in dbounds.items():
            mask[(data[var] <= low) | (data[var] >= high)] = True
        if flags:
            for flgname, goodval in flag_dict.items():
                mask[data[flgname] != goodval] = True

        # Find longest span of valid (unmasked) data
        dummy_ma = np.ma.array(np.zeros([data.shape[0], ]), mask=mask)
        unmasked_slices = np.ma.clump_unmasked(dummy_ma) or [slice(0, 0), ]
        max_indx = np.argmax([s.stop - s.start for s in unmasked_slices])
        max_slice = unmasked_slices[max_indx]
        len_max_slice = max_slice.stop - max_slice.start

        # verify sufficient data length
        data_frac = len_max_slice / data.shape[0]
        if data_frac < rd_tol or len_max_slice < ad_tol:
            self.data_table = None
            raise ValueError(
                'HF Data read but rejected because the longest continuous run '
                'of valid data was too short on a relative (length data / '
                'total length) = {:.4} < dtol = {:.4}) AND/OR absolute basis '
                '(data length = {} < {})'
                ''.format(data_frac, rd_tol, len_max_slice, ad_tol))

        self.data_table = data[var_names][max_slice].copy()
        return

    def __getitem__(self, name):
        """Column-wise get without specifying data_table attribute"""
        return self.data_table[name]

    def __setitem__(self, name, value):
        """Column-wise set without specifying data_table attribute"""
        self.data_table[name] = value

    def truncate(self):
        """Truncate length of data_table view to largest power of 2."""
        truncate_len = 2 ** int(np.log2(self.data_table.shape[0]))
        self.data_table = self.data_table[:truncate_len]

    def qc_correct(self):
        """Adjust q and c data series to correct for external effects.

        References
        ----------
        M. Detto, G.G. Katul (2007). Simplified expression for adjusting
        higher-order turbulent statistics obtained from open path gas
        analyzers. Boundary-Layer Meterol, 122:205-216.
        DOI:10.1007/s10546-006-9105-1.

        E.K. Webb, G.I. Pearman, R. Leuning (1980). Correction of flux
        measurements for density effects due to head and water vapor
        transfer. Quart J R Met Soc, 106:85-100.
        """

        if self._qc_already_corrected:
            return
        ave_vapor = self['q'].mean()
        ave_co2 = self['c'].mean()
        ave_T = self['T'].mean()
        dev_vapor = self['q'] - ave_vapor
        dev_T = self['T'] - ave_T

        Pdryair = self['P'].mean() - ave_vapor * GC.vapor * ave_T
        rho_totair = ave_vapor + Pdryair / GC.dryair / ave_T

        specific_vapor = ave_vapor / rho_totair
        specific_co2 = ave_co2 / rho_totair
        mu = MW.dryair / MW.vapor
        muq = mu * specific_vapor
        muc = mu * specific_co2

        self['q'] += muq * dev_vapor + (1 + muq) * ave_vapor * dev_T / ave_T
        self['c'] += muc * dev_vapor + (1 + muq) * ave_co2 * dev_T / ave_T
        self._qc_already_corrected = True
        return

    def summarize(self):
        """Return statisitical summary of high frequency data table.

        Returns
        -------
        HFSummary container
        """

        hfs = util.stats2(self.data_table, self.names)
        Pvap = hfs.ave_q * GC.vapor * hfs.ave_T
        rho_dryair = (hfs.ave_P - Pvap) / GC.dryair / hfs.ave_T
        rho_totair = rho_dryair + hfs.ave_q
        Cp = CP.dryair * (1 + 0.84 * hfs.ave_q / rho_totair)

        return HFSummary(
            T=hfs.ave_T,
            P=hfs.ave_P,
            Pvap=Pvap,
            ustar=(hfs.cov_w_u**2 + hfs.cov_w_v**2)**0.25,
            wind_w=hfs.ave_w,
            var_w=hfs.var_w,
            rho_vapor=hfs.ave_q,
            var_vapor=hfs.var_q,
            rho_co2=hfs.ave_c,
            var_co2=hfs.var_c,
            H=rho_totair * Cp * hfs.cov_w_T,
            Fq=hfs.cov_w_q,
            LE=util.qflux_mass_to_heat(hfs.cov_w_q, hfs.ave_T),
            Fc=hfs.cov_w_c,
            Fc_mol=util.cflux_mass_to_mol(hfs.cov_w_c),
            rho_dryair=rho_dryair,
            rho_totair=rho_totair,
            cov_w_T=hfs.cov_w_T,
            N=self.data_table.shape[0])


if __name__ == "__main__":
    pass
