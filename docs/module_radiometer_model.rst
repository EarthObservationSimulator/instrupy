.. _radiometer_model_module:

``instrupy.radiometer_model`` --- Radiometer Model
*****************************************************

Description
^^^^^^^^^^^^^
This module provides the ``RadiometerModel`` class (and supporting classes) which can be used to represent microwave radiometers.
The module consists of separate classes for handling each of the radiometric systems (listed in ``SystemType``), all with identical interface functions.
Similarly there are separate classes to handle different scan techniques (listed in ``ScanTech``), all with identical interface functions.  
By having identical interface functions, the functions can be invoked without consideration of the underlying radiometric system class.

The ``RadiometerModel`` object can be obtained from json/ dict representations using the ``from_json`` or ``from_dict`` functions. 
Please refer to :ref:`radiometer model<radiometer_model_desc>` for description of the instrument model and the expected key/value pairs.

Refer to the :ref:`base module<base_module>` and OrbitPy coverage calculator module (:class:`orbitpy.coveragecalculator`) to know how to use the ``RadiometerModel``
object for coverage calculations.

The ``calc_data_metrics`` function can be used to calculate the data-metrics given the spacecraft-state and observed-location (target).
The expected key/value pairs for the spacecraft-state input argument (of type :class:`dict`) are:

.. csv-table:: 
   :header: Column, Data type, Units, Description
   :widths: 30,10,10,40

   *time [JDUT1]*, float, days, Time in Julian Day UT1. Corresponds to the time of observation. 
   "*x [km]*, *y [km]*, *z [km]*", float, km, Cartesian coordinates of satellite in EARTH_CENTERED_INERTIAL frame at the time of observation.
   "*vx [km/s]*, *vy [km/s]*, *vz [km/s]*", float, km/s, Velocity of spacecraft in EARTH_CENTERED_INERTIAL frame at the time of observation.

The expected key/value pairs for the observed-location input argument (of type :class:`dict`) are:

.. csv-table:: 
   :header: Column, Data type, Units, Description
   :widths: 30,10,10,40

   "*lat [deg]*, *lon [deg]*", float, degrees, "Ground-point accessed (geocentric latitude, geocentric longitude)."



Examples
^^^^^^^^^


API
^^^^^

.. rubric:: Classes

.. autosummary::
   :nosignatures:
   :toctree: generated/
   :template: classes_template.rst
   :recursive:

   instrupy.radiometer_model.SystemType
   instrupy.radiometer_model.TotalPowerRadiometerSystem
   instrupy.radiometer_model.UnbalancedDikeRadiometerSystem
   instrupy.radiometer_model.BalancedDikeRadiometerSystem
   instrupy.radiometer_model.NoiseAddingRadiometerSystem
   instrupy.radiometer_model.ScanTech
   instrupy.radiometer_model.FixedScan
   instrupy.radiometer_model.CrossTrackScan
   instrupy.radiometer_model.ConicalScan
   instrupy.radiometer_model.RadiometerModel

.. rubric:: Functions

.. autosummary::
   :nosignatures:
   :toctree: generated/
   :template: functions_template.rst
   :recursive:

   instrupy.radiometer_model.PredetectionSectionParams
   instrupy.radiometer_model.SystemParams