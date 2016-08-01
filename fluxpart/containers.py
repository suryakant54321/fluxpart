from collections import namedtuple

SiteData = namedtuple('SiteData', 'meas_ht canopy_ht ppath')
SiteData.__doc__ += "Field site meta data."
SiteData.meas_ht.__doc__ = "Eddy covariance measurement height, m"
SiteData.canopy_ht.__doc__ = "Canopy height, m"
SiteData.ppath.__doc__ = "Photosynthetic pathway {'C3' or 'C4'}"


FluxComponents = namedtuple('FluxComponents', 'wq wqt wqe wc wcp wcr')
FluxComponents.__doc__ += "Water vapor and carbon dioxide flux components"
FluxComponents.wq.__doc__ = "Total water vapor flux : float, kg/m^2/s"
FluxComponents.wqt.__doc__ = "Transpiration vapor flux : float, kg/m^2/s"
FluxComponents.wqe.__doc__ = "Evaporation vapor flux : float, kg/m^2/s"
FluxComponents.wc.__doc__ = "Total CO2 flux : float, kg/m^2/s"
FluxComponents.wcp.__doc__ = "Photosynthesis CO2 flux : float, kg/m^2/s"
FluxComponents.wcr.__doc__ = "Respiration CO2 flux : float, kg/m^2/s"


class Result(namedtuple('Result', 'dataread valid_partition mssg')):
    """Overall outcome of partitioning.

    Attributes
    ----------
    dataread, valid_partition : bool
        Indicates success/failure in reading high frequency data
        (`dataread`) and obtaining a valid partioning solution 
        (`valid_partition`)
    mssg : str
        Possibly informative message if `dataread` or`valid_partition`
        are False

    """

    slots = ()

    def __str__(self):
        return ('Result(\n'
                '    dataread = {},\n'
                '    valid_partition = {},\n'
                '    mssg = {})'
                ''.format(*self))


class NumerSoln(namedtuple('NumerSoln', 'corr_cp_cr var_cp sig_cr co2soln_id '
                           'validroot validmssg init success mssg nfev')):
    """Results of numerical root finding.
    
    The sought root is (corr_cp_cr, var_cp).

    Attributes
    ----------
    corr_cp_cr : float
        Correlation coefficient for CO2 concentrations connected with
        photosynthesis (cp) and respiration (cr).
    var_cp : float, (kg/m^3)^2
        Variance of CO2 concentration connected with photosynthesis.
    sig_cr : float, kg/m^3
        Standard deviation of CO2 concentration connected with
        respiration.
    co2soln_id : {0 or 1}
        Indicates the solution used for the quadratic Eq. 13b of [SS08].
        Equal to 1 for the '+' root, equal to 0 for the '-' root.
    validroot : bool
        Indicates whether the obtained root (`corr_cp_cr`, `var_cp) is
        consistent with the physical model underlying the partitioning
        algorithm.
    validmssg : str
        Possibly informative message if `validroot` = False.
    init : (float, float)
        Value used to initialize root finding routine.
    success, mssg, nfev : bool, str, int
        Status information passed from scipy.optimize.root.

    """

    slots = ()

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


class WUE(namedtuple('WUE',
                     'wue inter_h2o inter_co2 ambient_h2o ambient_co2 vpd '
                     'ci_mod ppath meas_ht canopy_ht')):
    """Summary of leaf-level water use efficiency calculation.
    
    Attributes
    ----------
    wue : float, kg CO2 / kg H2O
        Leaf-level water use efficiency.
    inter_h2o, inter_co2, ambient_h2o, ambient_co2 : float, kg/m^3
        Concentrations of intercellular water vapor (`inter_h2o`),
        intercellular CO2 (`inter_co2`), ambient atmospheric water vapor
        (`ambient_h2o`), and ambient CO2 (`ambient_co2`).
    vpd : float, Pa
        Atmospheric vapor pressure deficit.
    ci_mod : (str, float) or (str, (float,float))
        str = {'const_ppm', 'const_ratio', 'linear', 'sqrt'} indicates
        the model used to estimate the intercellular CO2 concentration. 
        float or (float, float) are the model parameter values used.
    ppath : {'C3' or 'C4'}
        Photosynthetic pathway.
    meas_ht, canopy_ht : float, m
        Eddy covariance measument (`meas_ht`) and plant canopy
        (`canopy_ht`) heights.

    """

    slots = ()

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


class HFSummary(namedtuple('HFSummary', 'T P Pvap ustar wind_w var_w '
                           'rho_vapor rho_co2 var_vapor var_co2 Fq Fc H LE '
                           'Fc_mol rho_dryair rho_totair cov_w_T N')):

    """Summary of high frequency eddy covariance data.

    Attributes
    ----------
    T : float, K
        Mean air temperature.
    P, Pvap = float, Pa
        Mean total atmospheric (`P`) and vapor (`Pvap`) pressure.
    ustar : float, m/s
        Mean friction velocity.
    wind_w : float, m/s
        Mean vertical wind velocity.
    var_w : float, (m/s)^2
        Variance of vertical wind velocity
    rho_vapor, rho_co2 : float, kg/m^3
        Mean vapor (`rho_vapor`) and CO2 (`rho_co2`) concentrations.
    var_vapor, var_co2 : float, (kg/m^3)^2
        Variance of vapor (`var_vapor`) and CO2 (`var_co2`)
        concentrations.
    Fq, Fc = float, kg/m^2/s
        Water vapor (`Fq`) and CO2 (`Fc`) mass fluxes.
    H, LE : float, W/m^2
        Sensible (`H`) and latent (`LE`) heat fluxes.
    Fc_mol = mol/m^2/s
        Carbon dioxide mass flux expressed in mols.
    rho_dryair, rho_totair : float, kg/m^3
        Dry (`rho_dryair`) and moist (`rho_totair`) air densities.
    cov_w_T : float, K m/s
        Covariance of temperature and vertical wind velocity.
    N : int
        Length of data series.

    """

    slots = ()

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


class Fluxes(namedtuple('Fluxes', 'Fq Fqt Fqe Fc Fcp Fcr LE LEt LEe Fc_mol '
                        'Fcp_mol Fcr_mol')):

    """Water vapor and CO2 fluxes.

    Attributes
    ----------
    Fq, Fqt, Fqe : float, kg/m^2/s
        Water vapor mass fluxes: `Fq` = total, `Fqt` = transpiration,
        and`Fqe` = evaporation.
        flux.
    LE, LEt, LEe : float, W/m^2
        The same vapor fluxes expressed as latent heat.
    Fc, Fcp, Fcr : float, kg/m^2/s
        CO2 mass fluxes: `Fc` = total; `Fcp` = photosynthesis, and
        `Fcr` = respiration.
    Fc_mol, Fcp_mol, Fcr_mol : float, mol/m^2/s
        The same CO2 fluxes expressed as mols.

    """

    slots = ()

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
                ''.format(*Fq, *Fc, *LE, *Fc_mol))


class QCData(namedtuple('QCData', 'var_q var_c corr_qc wq wc wave_lvl')):

    """Summary stats for water vapor and CO2 series data.

    Attributes
    ----------
    var_q, var_c : float, (kg/m^3)^2"
        Variance of vapor (`var_q`) and CO2 (`var_c`) concentrations.
    corr_qc : float
        Correlation coefficient for vapor and CO2 concentrations.
    wq, wc : float, kg/m^2/s"
        Mean vapor (`wq`) and CO2 (`wc`) fluxes.
    wave_lvl : (int, int)
        2-tuple indicating the level of filtering applied. The second
        int is the maximum possible wavelet decompostion level given
        the length of the data. In other words, it is the number of
        'details' that could be removed from the data. The first int is
        a counter that starts equal the second int and decrements as
        details are removed. So when the first number is equal to the
        second, no details have been removed (no filtering applied).
        When the first number is 1, the maximum level of filtering was
        applied.

    """

    slots = ()

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
