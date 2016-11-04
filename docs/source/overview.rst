.. _fluxpart-overview:

========
Overview
========

The eddy covariance method is commonly used for measuring gas fluxes over
agricultural fields and natural ecosystems.
For many applications, it is desirable to partition the measured fluxes into
constitutive components: the water vapor flux into transpiration and direct
evaporation components, and the carbon dioxide flux in photosynthesis and
respiration components.
:cite:t:`Scanlon+Sahu:2008` devised a partitioning method based on flux
variance similarity relationships and correlation analyses of high-frequency
eddy covariance measurements.
Fluxpart is a free and open-source Python module that implements a version
of the :cite:t:`Scanlon+Sahu:2008` flux partitioning procedure.
High frequency eddy covariance data are processed and estimates of the various
flux compoents are obtained.
