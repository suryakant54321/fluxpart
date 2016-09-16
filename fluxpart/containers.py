"""Namedtuples used to group data and results."""

from collections import namedtuple

class SiteData(namedtuple('SiteData', 'meas_ht canopy_ht ppath')):
    """Field site meta data.

    Attributes
    ----------
    meas_ht, canopy_ht : float
        Eddy covariance measurement height and canopy height, m
    ppath : {'C3' or 'C4'}
        Photosynthetic pathway

    """

    __slots__ = ()

class FluxComponents(namedtuple('FluxComponents', 'wq wqt wqe wc wcp wcr')):
    """Vapor and CO2 flux components.

    Attributes
    ----------
    wq, wqt, wqe : float
        Total, transpiration, and evaporation H2O fluxes, kg/m^2/s
    wc, wcp, wcr : float
        Total, photosynthesis, and respiration CO2 fluxes, kg/m^2/s

    """

    __slots__ = ()

    def __str__(self):
        # Print with common mass units instead of SI
        wqs = [1e3 * self.wq, 1e3 * self.wqt, 1e3 * self.wqe]
        wcs = [1e6 * self.wc, 1e6 * self.wcp, 1e6 * self.wcr]
        return ('FluxComponents(\n'
                '    wq = {:.4} g/m^2/s,\n'
                '    wqt = {:.4} g/m^2/s,\n'
                '    wqe = {:.4} g/m^2/s,\n'
                '    wc = {:.4} mg/m^2/s,\n'
                '    wcp = {:.4} mg/m^2/s,\n'
                '    wcr = {:.4} mg/m^2/s)'
                ''.format(*(wqs + wcs)))


class Fluxes(namedtuple('Fluxes', 'Fq Fqt Fqe Fc Fcp Fcr LE LEt LEe Fc_mol '
                        'Fcp_mol Fcr_mol')):

    """Water vapor and CO2 fluxes.

    Attributes
    ----------
    Fq, Fqt, Fqe : float
        Total, transpiration, and evaporation H2O fluxes, kg/m^2/s.
    LE, LEt, LEe : float
        The same vapor fluxes expressed as latent heat, W/m^2.
    Fc, Fcp, Fcr : float
        Total, photosynthesis, and respiration CO2 fluxes, kg/m^2/s.
    Fc_mol, Fcp_mol, Fcr_mol : float
        The same CO2 fluxes expressed as mol/m^2/s.

    """

    __slots__ = ()

    def __str__(self):

        # For some fields, print common units instead of SI
        Fq, Fc, LE, Fc_mol = list(self[i:i + 3] for i in range(0, 12, 3))
        Fq = list(1e3 * f for f in Fq)
        Fc = list(1e6 * f for f in Fc)
        Fc_mol = list(1e6 * f for f in Fc_mol)

        return ('Fluxes(\n'
                '    Fq = {:.4} g/m^2/s,\n'
                '    Fqt = {:.4} g/m^2/s,\n'
                '    Fqe = {:.4} g/m^2/s,\n'
                '    Fc = {:.4} mg/m^2/s,\n'
                '    Fcp = {:.4} mg/m^2/s,\n'
                '    Fcr = {:.4} mg/m^2/s,\n'
                '    LE = {:.4} W/m^2,\n'
                '    LEt = {:.4} W/m^2,\n'
                '    LEe = {:.4} W/m^2,\n'
                '    Fc_mol = {:.4} umol/m^2/s,\n'
                '    Fcp_mol = {:.4} umol/m^2/s,\n'
                '    Fcr_mol = {:.4} umol/m^2/s)'
                ''.format(*(Fq + Fc + LE + Fc_mol)))


class HFSummary(namedtuple('HFSummary', 'T P Pvap ustar wind_w var_w '
                           'rho_vapor rho_co2 var_vapor var_co2 Fq Fc H LE '
                           'Fc_mol rho_dryair rho_totair cov_w_T N')):

    """Summary of high frequency eddy covariance data.

    Attributes
    ----------
    T : float
        Mean air temperature, K.
    P, Pvap : float
        Mean total atmospheric (`P`) and vapor (`Pvap`) pressure, Pa.
    ustar : float
        Mean friction velocity, m/s.
    wind_w : float
        Mean vertical wind velocity, m/s.
    var_w : float
        Variance of vertical wind velocity, (m/s)^2.
    rho_vapor, rho_co2 : float
        Mean H2O vapor and CO2 concentrations, kg/m^3.
    var_vapor, var_co2 : float
        Variance of H2O vapor and CO2 concentrations, (kg/m^3)^2.
    Fq, Fc : float
        Water vapor (q) and CO2 (c) mass fluxes, kg/m^2/s.
    H, LE : float
        Sensible and latent heat fluxes, W/m^2.
    Fc_mol : float
        Carbon dioxide mass flux expressed in mol/m^2/s.
    rho_dryair, rho_totair : float
        Dry and moist air densities, kg/m^3.
    cov_w_T : float
        Covariance of temperature and vertical wind velocity, K m/s.
    N : int
        Length of data series.

    """

    __slots__ = ()

    def __str__(self):

        # For some fields, print common units instead of SI
        dum = self._replace(
            T=self.T - 273.15,
            P=1e-3 * self.P,
            Pvap=1e-3 * self.Pvap,
            rho_vapor=1e3 * self.rho_vapor,
            rho_co2=1e6 * self.rho_co2,
            var_vapor=1e6 * self.var_vapor,
            var_co2=1e12 * self.var_co2,
            Fq=1e3 * self.Fq,
            Fc=1e6 * self.Fc,
            Fc_mol=1e6 * self.Fc_mol)

        return ('HFSummary(\n'
                '    T = {:.4} C,\n'
                '    P = {:.4} kPa,\n'
                '    Pvap = {:.4} kPa,\n'
                '    ustar = {:.4} m/s,\n'
                '    wind_w = {:.4} m/s,\n'
                '    var_w = {:.4} (m/s)^2,\n'
                '    rho_vapor = {:.4} g/m^3,\n'
                '    rho_co2 = {:.4} mg/m^3,\n'
                '    var_vapor = {:.4} (g/m^3)^2,\n'
                '    var_co2 = {:.4} (mg/m^3)^2,\n'
                '    Fq = {:.4} g/m^2/s,\n'
                '    Fc = {:.4} mg/m^2/s,\n'
                '    H = {:.4} W/m^2,\n'
                '    LE = {:.4} W/m^2,\n'
                '    Fc_mol = {:.4} umol/m^2/s,\n'
                '    rho_dryair = {:.4} kg/m^3,\n'
                '    rho_totair = {:.4} kg/m^3,\n'
                '    cov_w_T = {:.4} C m/s,\n'
                '    N = {})'
                ''.format(*dum))


class NumerSoln(namedtuple('NumerSoln', 'corr_cp_cr var_cp sig_cr co2soln_id '
                           'validroot validmssg init success mssg nfev')):
    """Results of numerical root finding.

    The sought root is (corr_cp_cr, var_cp).

    Attributes
    ----------
    corr_cp_cr : float
        Correlation coefficient for CO2 concentrations connected with
        photosynthesis (cp) and respiration (cr).
    var_cp : float
        Variance of CO2 concentration connected with photosynthesis,
        (kg/m^3)^2.
    sig_cr : float
        Standard deviation of CO2 concentration connected with
        respiration, kg/m^3.
    co2soln_id : {0 or 1}
        Indicates the solution used for the quadratic Eq. 13b of
        [SS08]_. Equal to 1 for the '+' root, equal to 0 for the '-'
        root.
    validroot : bool
        Indicates whether the obtained root (`corr_cp_cr`, `var_cp) is
        physically plausible.
    validmssg : str
        Possibly informative message if `validroot` = False.
    init : (float, float)
        (corr_cp_cr, var_cp) pair used to initialize root finding
        routine.
    success, mssg, nfev : bool, str, int
        Status information passed from scipy.optimize.root_.


    .. _scipy.optimize.root:
        http://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.root.html

    """

    __slots__ = ()

    def __str__(self):

        # For some fields, print common units instead of SI
        dum = self._replace(
            var_cp=1e12 * self.var_cp,
            sig_cr=1e6 * self.sig_cr,
            init=(self.init[0], 1e12 * self.init[1]))

        return ('NumerSoln(\n'
                '    corr_cp_cr = {0:.4},\n'
                '    var_cp = {1:.4} (mg/m^3)^2,\n'
                '    sig_cr = {2:.4} mg/m^3,\n'
                '    co2soln_id = {3},\n'
                '    validroot = {4},\n'
                '    validmssg = {5},\n'
                '    init = ({6[0]:.4}, {6[1]:.4}),\n'
                '    success = {7},\n'
                '    mssg = {8},\n'
                '    nfev = {9})'
                ''.format(*dum))


class QCData(namedtuple('QCData', 'var_q var_c corr_qc wq wc wave_lvl')):

    """Summary stats for water vapor and CO2 series data.

    Attributes
    ----------
    var_q, var_c : float
        Variance of vapor (q) and CO2 (c) concentrations, (kg/m^3)^2.
    corr_qc : float
        Correlation coefficient for vapor and CO2 concentrations.
    wq, wc : float
        Mean vapor (`wq`) and CO2 (`wc`) fluxes, kg/m^2/s.
    wave_lvl : (int, int)
        2-tuple indicating the level of filtering applied (number of
        components removed from the series data). The second int is the
        maximum possible wavelet decompostion level given the length of
        the data. The first is the number of components remaining in the
        data.  So when the first number is equal to the second, no
        components have been removed (no filtering applied). When the
        first number is 1, the maximum level of filtering was applied.

    """

    __slots__ = ()

    def __str__(self):

        # For some fields, print common units instead of SI
        dum = self._replace(
            var_q=1e6 * self.var_q,
            var_c=1e12 * self.var_c,
            wq=1e3 * self.wq,
            wc=1e6 * self.wc)

        return ('QCData(\n'
                '    var_q = {:.4} (g/m^3)^2,\n'
                '    var_c = {:.4} (mg/m^3)^2,\n'
                '    corr_qc = {:.4},\n'
                '    wq = {:.4} g/m^2/s,\n'
                '    wc = {:.4} mg/m^2/s,\n'
                '    wave_lvl = {})'
                ''.format(*dum))


class Result(namedtuple('Result', 'dataread valid_partition mssg')):
    """Overall outcome of partitioning.

    Attributes
    ----------
    dataread, valid_partition : bool
        Indicates success or failure in reading high frequency data and
        obtaining a valid partioning solution.
    mssg : str
        Possibly informative message if `dataread` or `valid_partition`
        are False

    """

    __slots__ = ()

    def __str__(self):
        return ('Result(\n'
                '    dataread = {},\n'
                '    valid_partition = {},\n'
                '    mssg = {})'
                ''.format(*self))


class WUE(namedtuple('WUE',
                     'wue inter_h2o inter_co2 ambient_h2o ambient_co2 vpd '
                     'ci_mod ppath meas_ht canopy_ht')):
    """Summary of leaf-level water use efficiency calculation.

    Attributes
    ----------
    wue : float
        Leaf-level water use efficiency, kg CO2 / kg H2O.
    inter_h2o, inter_co2, ambient_h2o, ambient_co2 : float
        Concentrations of intercellular water vapor, intercellular CO2,
        ambient atmospheric water vapor, and ambient CO2, kg/m^3.
    vpd : float
        Atmospheric vapor pressure deficit, Pa.
    ci_mod : (str, float) or (str, (float,float))
        The str is the name of the model used to estimate the
        intercellular CO2 concentration, and float or (float, float) are
        the model parameter values used.
    ppath : {'C3' or 'C4'}
        Photosynthetic pathway.
    meas_ht, canopy_ht : float
        Eddy covariance measument height and plant canopy height, m.

    """

    __slots__ = ()

    def __str__(self):

        # For some fields, print common units instead of SI
        dum = self._replace(
            wue=1e3 * self.wue,
            inter_h2o=1e3 * self.inter_h2o,
            inter_co2=1e6 * self.inter_co2,
            ambient_h2o=1e3 * self.ambient_h2o,
            ambient_co2=1e6 * self.ambient_co2,
            vpd=1e-3 * self.vpd)

        return ('WUE(\n'
                '    wue = {:.4} mg/g,\n'
                '    inter_h2o = {:.4} g/m^3,\n'
                '    inter_co2 = {:.4} mg/m^3,\n'
                '    ambient_h2o = {:.4} g/m^3,\n'
                '    ambient_co2 = {:.4} mg/m^3,\n'
                '    vpd = {:.4} kPa,\n'
                '    ci_mod = {},\n'
                '    ppath = {},\n'
                '    meas_ht = {} m,\n'
                '    canopy_ht = {} m)'
                ''.format(*dum))
