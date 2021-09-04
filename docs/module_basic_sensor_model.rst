.. _basic_sensor_model_module:

``instrupy.basic_sensor_model`` --- Basic Sensor Model
========================================================

Description
^^^^^^^^^^^^^


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