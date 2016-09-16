from math import log, exp, sqrt

from fluxpart.constants import GRAVITY, VON_KARMAN
from fluxpart.constants import MOLECULAR_WEIGHT as MW
from fluxpart.constants import SPECIFIC_GAS_CONSTANT as Rgas
from fluxpart.containers import WUE


# Defalut parameter values intercellular CO2 (ci) models
CI_DEFAULT_PARAMS = {
    'C3': {'const_ppm': 280.,    # ppm
           'const_ratio': 0.7,
           'linear': (1, 2e-2),
           'sqrt': 22e-9},       # kg-co2 / m^3 / Pa

    'C4': {'const_ppm': 130.,
           'const_ratio': 0.44,
           'linear': (1, 2e-2),
           'sqrt': 22e-9}}


def water_use_efficiency(hfs, site, ci_mod='const_ratio'):

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
    site : :class:`~fluxpart.containners.SiteData` tuple or equivalent namespace
        Container holding field-site meta data, possessing the following
        attributes:
        `meas_ht`, eddy covariance measurement height (m);
        `canopy_ht`, vegetation canopy height (m);
        `ppath`, photosynthetic pathway ('C3' or 'C4')
    ci_mod : str or [str, float] or [str, (float, float)], optional
        Specifies the model to be used to determine the leaf
        intercellular CO2 concentration. The str is the model name
        and is one fo the following (see Notes below for definitions),
        ``str = {'const_ratio', 'const_ppm', 'linear', 'sqrt'}``
        If `ci_mod` is a 2-list (or 2-tuple), the first element is str
        and the second is the model parameter value(s) (see Notes). If
        the model has more than one parameter, the second element is a
        tuple holding the values. If only str is provided, default
        parameter values are used (see Notes). Default is `ci_mod` =
        'const_ratio'.

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
        ci is determined from a specified constant ppm value, ci_ppm::

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
        The parameter pair (b, m) defaluts to (1, 2e-2) for C3 plants
        and (1, TODO) for C4. See [MG83]_.

    'sqrt'
        ci/ca is a function of sqrt(vpd/ca)::

            ci/ca = 1 - sqrt(1.6 * lambd *  vpd / ca)

        The paramater lambd has units of kg-CO2 / m^3 / Pa, and defaults
        to 22e-9 for C3 plants and TODO for C4. See [KP09]_.

    """

    # Assume zero-plane and roughness params for vapor and CO2 are the same.
    # d0 = Zero-plane displacement height (L), Eq. 5.2 of [CN98]
    d0 = 0.65 * site.canopy_ht
    # zv = Roughness parameter (L), Eqs. 5.3 and 7.19 of [CN98]
    zv = 0.2 * (0.1 * site.canopy_ht)

    # Obukhov length
    Qflux = hfs.cov_w_T + 0.61 * hfs.T * hfs.Fq / hfs.rho_totair  # K m/s
    obukhov_len = -hfs.T * hfs.ustar ** 3 / (VON_KARMAN * GRAVITY * Qflux)  # m

    # Stability correction
    zeta = (site.meas_ht - d0) / obukhov_len
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
    arg = (log((site.meas_ht - d0) / zv) - psi_v) / VON_KARMAN / hfs.ustar
    ambient_h2o = (hfs.rho_vapor + hfs.Fq * arg)
    ambient_co2 = (hfs.rho_co2 + hfs.Fc * arg)

    # Saturation vapor pressure `esat` & vapor pressure deficit `vpd`
    Tr = 1 - 373.15 / hfs.T
    esat = 101325. * (
        exp(13.3185 * Tr - 1.9760 * Tr**2 - 0.6445 * Tr**3 - 0.1299 * Tr**4))
    vpd = esat - hfs.rho_vapor * Rgas.vapor * hfs.T

    # Intercellular vapor density. Assume leaf temperature = air temperature.
    eps = MW.vapor / MW.dryair
    inter_h2o = hfs.rho_totair * eps * esat / (hfs.P - (1 - eps) * esat)

    # Intercellular CO2 concentration, aka ci (kg/m^3)
    if isinstance(ci_mod, str):
        ci_mod = (ci_mod, None)

    ci_mod_name = ci_mod[0]
    ci_mod_params = ci_mod[1] or CI_DEFAULT_PARAMS[site.ppath][ci_mod_name]
    ci_dispatch = {
        'const_ppm': _ci_const_ppm(hfs.P, hfs.T, Rgas.co2),
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
        ppath=site.ppath,
        meas_ht=site.meas_ht,
        canopy_ht=site.canopy_ht,
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
