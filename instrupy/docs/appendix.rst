Appendix
*********

.. _satellite_to_target_viewing_geometry:

Satellite to Target viewing geometry
=============================================

.. figure:: space_viewing_geometry_SandiaReport.png

*   :math:`\mathbf{R = T - S}`
*   :math:`\gamma = \cos^{-1}(\mathbf{\dfrac{R}{|R|}} \cdot \mathbf{\dfrac{-S}{|S|}})`
*   :math:`\theta_i = \sin^{-1}(\sin\gamma  \hspace{1mm}  \dfrac{R_E + h}{R_E})`

    .. note:: This formula assumes spherical Earth of radius :math:`R_E`

where,

* :math:`\mathbf{S}`: Position vector of the satellite in the Earth-Centered-Inertial frame (equatorial-plane)
* :math:`\mathbf{T}`: Position vector of the target ground-point in the Earth-Centered-Inertial frame (equatorial-plane)
* :math:`\mathbf{R}`: Range vector from satellite to target ground point
* :math:`\gamma`:  Look-angle to target ground point from satellite
* :math:`\theta_i`: Incidence angle at the target ground point
* :math:`R_E`: Nominal equatorial radius of Earth
* :math:`h`: altitude of satellite


.. _manuv_desc:

Manuverability description
=============================================


.. todo: Add Solar elevation illustration

  
Glossary of terms used in the package
======================================

Pixels vs Detectors

Pixels: Refer to ground pixels imaged. Dimensions vary according to imaging geometry.

Detectors: Refer to physical detector elements on the imaging aperture.

Access vs Coverage

satellite, spacecraft

target

observation incidence angle for the case of non-radars


Coding Conventions
===================

* variables denoting physical quantities, unless otherwise indicated are always in S.I. units.



