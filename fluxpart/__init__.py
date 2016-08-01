"""
Module for partitioning fluxes into stomatal & nonstomatal components.

The partitioning method is from [SS08] and is based on flux
variance similarity assumptions.

The method uses as input interval averaged (e.g., 30 min) EC data for:
    <w'q'> : water vapor flux [g m^-2 s^-1]
    <w'c'> : carbon dioxide flux [mg m^-2 s^-1]
    var_q  : variance of water vapor concentration [(g m^-3)^2]
    var_c  : variance of carbon dioxide concentration [(mg m^-3)^2]
    wue    : water use efficiency (stomatal exchange) [mg g^-1].

A successful outcome obtains the components of the water vapor and
carbon dioxide fluxes that are associated with stomatal and nonstomatal
mechanisms:
    <w'qe'> : nonstomatal H2O component (evaporation) [g m^-2 s^-1]
    <w'qt'> : stomatal H2O component (transpiration) [g m^-2 s^-1]
    <w'cr'> : nonstomatal CO2 component (repsiration) [mg m^-2 s^-1]
    <w'cp'> : stomatal CO2 component (photosynthesis) [mg m^-2 s^-1]

A solution is found by solving two simultaneous nonlinear equations in
which the two free parameters are (i.) the variance of the carbon
dioxide concentration associated with photosynthesis, var_cp, and
(ii.) the the correlation coefficient for carbon dioxide concentrations
associated with photosynthesis and respiration, corr_cp_cr. See [XXXX].


Notes
-----
Basic physical relations, definitions, and requirements:
    q' = qe' + qt'
    c' = cr' + cp'
    cp' = wue * qt'
    <w'q'> = <w'qe'> + <w'qt'>
    <w'c'> = <w'cr'> + <w'cp'>
    <w'cp'>  < 0
    <w'cr'>  > 0
    <w'qt'>  > 0
    <w'qe'>  > 0
    <w'cr'> / <w'cp'>  < 0
    <w'qe'> / <w'qt'>  > 0
    wue < 0

References
----------
[SS08] Scanlon, T.M., and P. Sahu (2008), On the correlation structure
of water vapor and carbon dioxide in the atmospheric surface layer: A
basis for flux partitioning. Water Resour. Res., 44, W10418,
doi:10.1029/2008WR006932.

[PA14] Palatella, L., G. Rana, D. Vitale (2014), Towards a flux-
partitioning procedure based on the direct use of high-frequency eddy-
covariance data.  Boundary-Layer Meteorol. 153:327--337, 
doi:10.1007/s10546-014-9947-x
"""

from fluxpart.fluxpart import flux_partition
