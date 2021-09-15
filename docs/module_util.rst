``instrupy.util`` --- Utility classes, functions
************************************************

Description
============

This module contains the common classes and functions used by the instrument models. 
The ``Orientation``, ``SphericalGeometry``, ``ViewGeometry``, ``Maneuver``, ``Antenna`` and ``SyntheticDataConfiguration`` classes are
purposed for handling common instrument parameters. Objects of these classes can be obtained from json-strings or python-dictionaries. 
The expected key/value pairs is described :ref:`here<common_instru_params>`. 

:class:`instrupy.util.Orientation`
----------------------------------
Class to store and handle orientation of instruments or spacecrafts.

:class:`instrupy.util.SphericalGeometry`
-----------------------------------------
The :code:`SphericalGeometry` object can be used to handle the spherical-geometry specifications of field-of-view/ scene-field-of-view/ field-of-regard.

:class:`instrupy.util.ViewGeometry` 
--------------------------------------
Container class which congregates the :class:`SphericalGeometry` and :class:`Orientation` objects and can be used to model
the Field of View (FOV) or Scene FOV or Field of Regard (FOR) of a sensor. 

The :code:`SphericalGeometry` member of the container describes the spherical geometry (polygon/ circle) in the SENSOR_BODY_FIXED frame
with the Z-axis as the pointing axis. This can be paired with an Orientation object (which describes the orientation of the sensor (hence the SENSOR_BODY_FIXED frame)
with respect to a reference frame) to obtain the position of the spherical geometry in any desired reference frame.

.. note:: In the current :class:`instrupy` implementation when used to model the FOR, the Orientation is always defined with respect to the 
             NADIR_POINTING reference frame. 

:class:`instrupy.util.Maneuver` 
--------------------------------------
Class handling the maneuverability of the satellite-sensor system. The maneuverability is always specified with reference to the *NADIR_POINTING* frame. The maneuver specifications 
describe the angular-space where the pointing axis of the sensor can be positioned.

This class includes the function ``calc_field_of_regard`` which can be used to obtain the Field-Of-Regard in terms of a *proxy-FOV setup*. 
The proxy-sensor setup is characterized by orientation (wrt the *NADIR_POINTING* frame) of the proxy-sensor (hence the proxy-sensor *SENSOR_BODY_FIXED* frame)
and a spherical geometry (polygon/circle) specification of the proxy-sensor's field-of-view. This allows to calculate all coverage opportunities
by the satellite-sensor pair, taking into account the satellite and/or sensor maneuverability. 

API
======

.. rubric:: Classes

.. autosummary::
   :nosignatures:
   :toctree: generated/
   :template: classes_template.rst
   :recursive:

   instrupy.util.Entity
   instrupy.util.EnumEntity
   instrupy.util.Constants
   instrupy.util.Orientation
   instrupy.util.SphericalGeometry
   instrupy.util.ViewGeometry
   instrupy.util.Maneuver
   instrupy.util.Antenna
   instrupy.util.SyntheticDataConfiguration
   instrupy.util.SyntheticDataInterpolator
   instrupy.util.MathUtilityFunctions
   instrupy.util.GeoUtilityFunctions
   instrupy.util.FileUtilityFunctions


