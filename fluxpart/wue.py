from math import log, exp, sqrt

from fluxpart.constants import GRAVITY, VON_KARMAN
from fluxpart.constants import MOLECULAR_WEIGHT as MW
from fluxpart.constants import SPECIFIC_GAS_CONSTANT as Rgas
from fluxpart.containers import WUE


# Defalut parameter values intercellular CO2 (ci) models
CI_DEFAULT_PARAMS = {
    'C3': {'const_ppm': 280.,    # ppm
           'const_ratio': 0.7,
           'linear': (1, 1.6e-4),
           'sqrt': 22e-9},       # kg-co2 / m^3 / Pa

    'C4': {'const_ppm': 130.,
           'const_ratio': 0.44,
           'linear': (1, 2.7e-4)}}


def water_use_efficiency(hfs, meas_ht, canopy_ht, ppath, ci_mod, 
                         ci_mod_param=None, meas_leaf_temp=None,
                         leaf_temp_corr=0):

    """Estimate leaf-level water use efficiency.

    Parameters
    ----------
    hfs : :class:`~fluxpart.containers.HFSummary` tuple or equivalent namespace
        High frequency summary statistics, must possesses the following
        attributes (all floats):
        `rho_vapor`, mean vapor density (kg/m^3);
        `rho_co2`, mean carbon dioxide concentration (kg/m^3);
        `T`, mean air temperature (K);
        `P`, mean atmospheric pressure (Pa);
        `Fq`, mean water vapor flux (kg/m^2/s);
        `Fc`, mean carbon dioxide flux (kg/m^2/s);
        `cov_w_T`, covariance of wind and temperature (K m/s);
        `ustar`, friction velocity (m/s);
        `rho_totair`, moist air density (kg/m^3).
    meas_ht, canopy_ht : float 
        `meas_ht` = eddy covariance measurement height (m); 
        `canopy_ht` = vegetation canopy height (m);
    ppath : {'C3', 'C4'}
         photosynthetic pathway
    ci_mod : {'const_ratio', 'const_ppm', 'linear', 'sqrt'}``
        Specifies the model to be used to determine the leaf
        intercellular CO2 concentration. See Notes below for model
        descriptions.
    ci_mod_param : float or 2-tuple of floats, optional
        Paramter values to be used with `ci_mod`. The number of
        parameters required depends on the model (see Notes). The
        default is ci_mod_param=None, in which case default values are
        used
    meas_leaf_temp : float, optional
        Measured canopy tempeature in degrees K. If None (defalut), the 
        leaf temperature is taken to be equal to the air temperature in
        `hfs`.
    leaf_temp_corr : float, optional
        Optional adjustment to leaf temperature. The temperature used to
        calculate intercelluat vapor and CO2 concentrations is
        leaf_T + leaf_temp_corr, where leaf_T is `meas_leaf_temp` if
        provided, and the air temperature in `hfs` otherwise. Default is
        leaf_temp_corr = 0.

    Returns
    -------
    namedtuple
        :class:`~fluxpart.containers.WUE`

    Notes
    -----
    Leaf-level water use efficiency is estimated as [CN98]_ (pg. 238)::

        wue = 0.7 * (ca - ci) / (qa - qi)

    where::

        ca, ci = ambient and intercellular CO2 concentration, resp.
        qa, qi = ambient and intercellular vapor concentration, resp.

    ca, qa, and qi are caclulated as discussed by [SS08]_.

    To estimate ci, the following models are available:

    'const_ppm'
        ci (kg/m^3) is determined from a specified constant ppm value,
        ci_ppm::

            ci = f(ci_ppm; temperature, pressure)

        Default parameter values for ci_ppm are 280 ppm when ppath='C3',
        and 130 ppm when ppath='C4'. [CN98]_ (pg. 237).

    'const_ratio'
        The ratio ci/ca is constant::

            ci/ca = K.

        Default parameter values are K=0.7 for C3 plants and K=0.44 for
        C4.

    'linear'
        ci/ca is a linear function of vpd (the atmospheric vapor
        pressure deficit, which is calculated internally)::

            ci/ca = b - m * vpd

        b is dimensionless with a value of ~1 while m has units of 1/Pa.
        The parameter pair (b, m) defaluts to (1, 1.6e-4) for C3 plants
        and (1, 2.7e-4) for C4. See [MG83]_.

    'sqrt'
        ci/ca is a function of sqrt(vpd/ca) [KPO09]_::

            ci/ca = 1 - sqrt(1.6 * lambd *  vpd / ca)

        The paramater lambd has units of kg-CO2 / m^3 / Pa, and defaults
        to 22e-9 for C3 plants. The sqrt model is not enabled for C4
        plants.

    """

    # Assume zero-plane and roughness params for vapor and CO2 are the same.
    # d0 = Zero-plane displacement height (L), Eq. 5.2 of [CN98]
    d0 = 0.65 * canopy_ht
    # zv = Roughness parameter (L), Eqs. 5.3 and 7.19 of [CN98]
    zv = 0.2 * (0.1 * canopy_ht)

    # Virtual temperature flux
    Qflux = hfs.cov_w_T + 0.61 * hfs.T * hfs.Fq / hfs.rho_totair  # K m/s
    # Obukhov length
    obukhov_len = -hfs.T * hfs.ustar ** 3 / (VON_KARMAN * GRAVITY * Qflux)  # m

    # Stability correction
    zeta = (meas_ht - d0) / obukhov_len
    # Unstable
    if zeta < -0.04:
        psi_v = 2. * log((1 + (1 - 16. * zeta) ** 0.5) / 2)
    # Neutral
    elif zeta <= 0.04:
        psi_v = 0.
    # Stable
    else:
        psi_v = -5. * zeta

    # Ambient concentrations (kg/m^3)
    arg = (log((meas_ht - d0) / zv) - psi_v) / VON_KARMAN / hfs.ustar
    ambient_h2o = (hfs.rho_vapor + hfs.Fq * arg)
    ambient_co2 = (hfs.rho_co2 + hfs.Fc * arg)

    leaf_T = (meas_leaf_temp or hfs.T) + leaf_temp_corr
    # Intercellular saturation vapor pressure `esat`
    Tr = 1 - 373.15 / leaf_T
    esat = 101325. * (
        exp(13.3185 * Tr - 1.9760 * Tr**2 - 0.6445 * Tr**3 - 0.1299 * Tr**4))
    # Intercellular vapor pressure deficit `vpd`
    vpd = esat - hfs.rho_vapor * Rgas.vapor * leaf_T

    # Intercellular vapor density.
    eps = MW.vapor / MW.dryair
    inter_h2o = hfs.rho_totair * eps * esat / (hfs.P - (1 - eps) * esat)

    # Intercellular CO2 concentration, aka ci (kg/m^3)
    if isinstance(ci_mod, str):
        ci_mod = (ci_mod, None)

    ci_mod_name = ci_mod[0]
    if ci_mod_name == 'sqrt' and ppath == 'C4':
        err = "Combination of 'sqrt' ci model and 'C4' ppath not enabled"
        raise ValueError(err)
    ci_mod_params = ci_mod[1] or CI_DEFAULT_PARAMS[ppath][ci_mod_name]

    ci_dispatch = {
        'const_ppm': _ci_const_ppm(hfs.P, leaf_T, Rgas.co2),
        'const_ratio': _cica_const_ratio(ambient_co2),
        'linear': _cica_linear(ambient_co2, vpd),
        'sqrt': _cica_sqrt(ambient_co2, vpd)}
    inter_co2 = ci_dispatch[ci_mod_name](ci_mod_params)

    wue = 0.7 * (ambient_co2 - inter_co2) / (ambient_h2o - inter_h2o)

    return WUE(
        wue=wue,
        inter_h2o=inter_h2o,
        inter_co2=inter_co2,
        ambient_h2o=ambient_h2o,
        ambient_co2=ambient_co2,
        vpd=vpd,
        ppath=ppath,
        meas_ht=meas_ht,
        canopy_ht=canopy_ht,
        ci_mod=(ci_mod_name, ci_mod_params))


def _ci_const_ppm(pressure, temperature, Rco2):
    """ci is a fixed ppm value."""
    def ci_func(ci_ppm):
        """Return ci = intercellular CO2 concentration, kg/m^3."""
        return ci_ppm * 1e-6 * pressure / Rco2 / temperature
    return ci_func


def _cica_const_ratio(ambient_co2):
    """ci/ca is constant."""
    def ci_func(const):
        """Return ci = intercellular CO2 concentration, kg/m^3."""
        return const * ambient_co2
    return ci_func


def _cica_linear(ambient_co2, vpd):
    """ci/ca is a decreasing linear function of vapor pressure deficit."""
    def ci_func(linear_params):
        """Return ci = intercellular CO2 concentration, kg/m^3."""
        # b is unitless with a value of ~1, and m (> 0) has units of Pa^-1
        b, m = linear_params
        return ambient_co2 * (b - m * vpd)
    return ci_func


def _cica_sqrt(ambient_co2, vpd):
    """ci/ca is a function of sqrt(`vpd`/ca)."""
    def ci_func(lambd):
        """Return ci = intercellular CO2 concentration, kg/m^3."""
        # lambd has units of kg-co2 / m^3 / Pa"""
        return ambient_co2 * (1 + sqrt(1.6 * lambd * vpd / ambient_co2))
    return ci_func
