.. _basic_sensor_model_module:

``instrupy.basic_sensor_model`` --- Basic Sensor Model
========================================================

Description
^^^^^^^^^^^^^

This module provides the ``BasicSensorModel`` class which can be used to represent basic properties of an instrument.
The ``BasicSensorModel`` object can be obtained from json/ dict representations using the ``from_json`` or ``from_dict`` functions. 
Please refer to :ref:`basic sensor model<basic_sensor_model_desc>` for description of the instrument model and the expected key/value pairs.

Refer to the :ref:`base module<base_module>` and OrbitPy coverage calculator module (:class:`orbitpy.coveragecalculator`) to know how to use the ``BasicSensorModel``
object for coverage calculations.

The ``calc_data_metrics`` function can be used to calculate the data-metrics given the spacecraft-state and observed-location (target).
The expected key/value pairs for the spacecraft-state input argument (of type :class:`dict`) are:

.. csv-table:: 
   :header: Column, Data type, Units, Description
   :widths: 30,10,10,40

   *time [JDUT1]*, float, JDUT1, Time in Julian Day UT1. Corresponds to the time of observation. 
   "*x [km]*, *y [km]*, *z [km]*", float, km, Cartesian coordinates of satellite in EARTH_CENTERED_INERTIAL frame at the time of observation.
   "*vx [km/s]*, *vy [km/s]*, *vz [km/s]*", float, km/s, Velocity of spacecraft in EARTH_CENTERED_INERTIAL frame at the time of observation.

The expected key/value pairs for the observed-location input argument (of type :class:`dict`) are:

.. csv-table:: 
   :header: Column, Data type, Units, Description
   :widths: 30,10,10,40

   "*lat [deg]*, *lon [deg]*", float, degrees, "Ground-point accessed (geocentric latitude, geocentric longitude)."

Examples
^^^^^^^^^
1. Initializing a basic-sensor with 5 deg circular FOV geometry and aligned to spacecraft-body. Maneuver is possible within a 10 deg circular angular region.
   Unique identifier "bs1" is assigned.

   .. code-block:: python

      from instrupy.basic_sensor_model import BasicSensorModel
        
      bs1 = BasicSensorModel.from_json('{"name": "Atom", "mass":10, "volume":12.45, "dataRate": 40, "bitsPerPixel": 8, "power": 12, \
                                            "orientation": {"referenceFrame": "SC_BODY_FIXED", "convention": "REF_FRAME_ALIGNED"}, \
                                            "fieldOfViewGeometry": {"shape": "CIRCULAR", "diameter":5 }, \
                                            "maneuver":{"maneuverType": "CIRCULAR", "diameter":10} \
                                            "@id": "bs1"
                                                        }')

2. This example assumes the observer and the target location, and calculates the resulting viewing geometry. Note that none of the sensor parameters
   affect calculation of the viewing geometry directly. The parameters are used to determine if the target location can be viewed or not.

   .. code-block:: python

      from instrupy.basic_sensor_model import BasicSensorModel

      bs1 = BasicSensorModel.from_json('{}')  
      epoch_JDUT1 =  2458543.06088 # 2019 Feb 28 13:27:40 is time at which the ECEF and ECI frames approximately align, hence ECEF to ECI rotation is identity. See <https://www.celnav.de/longterm.htm> online calculator of GMST.
        
      SpacecraftOrbitState = {'time [JDUT1]':epoch_JDUT1, 'x [km]': 6878.137, 'y [km]': 0, 'z [km]': 0, 'vx [km/s]': 0, 'vy [km/s]': 7.6126, 'vz [km/s]': 0} # altitude 500 km
      TargetCoords = {'lat [deg]': 0, 'lon [deg]': 0}
      obsv_metrics = bs1.calc_data_metrics(SpacecraftOrbitState, TargetCoords)

      print(obsv_metrics)

      >> {'observation range [km]': 500.0, 'look angle [deg]': 0.03, 'incidence angle [deg]': 0.03, 'solar zenith [deg]': 20.33}

API
^^^^^

.. rubric:: Classes

.. autosummary::
   :nosignatures:
   :toctree: generated/
   :template: classes_template.rst
   :recursive:

   instrupy.basic_sensor_model.BasicSensorModel