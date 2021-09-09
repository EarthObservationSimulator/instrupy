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

``instru_look_angle_from_target_inc_angle`` flag
---------------------------------------------------

In addition to the above mentioned input-arguments to the ``calc_data_metrics`` function, there is an additional optional flag which can be 
specified to instruct how to obtain the instrument look-angle. If the ``instru_look_angle_from_target_inc_angle`` flag is set to ``False`` (also default value),
the instrument look-angle is obtained from the instrument orientation specification. If the flag is set to ``True``, the instrument look-angle is obtained from the target location incidence-angle. 

Examples
^^^^^^^^^

1. Define a total-power radiometer (with component-level predetection-stage specifications) with *FIXED* scan configuration and calculate metrics.
   The ``instru_look_angle_from_target_inc_angle`` flag is set to ``False`` in this example and hence it implies that the instrument look-angle is derived
   from the instrument orientation.
   Beam-efficiency for circular aperture antennas is not currently supported and hence is output as ``nan``. 

   .. note:: The swath-width is calculated more precisely as compared to the pixel-resolution calculations. This leads to a slight difference
             as can be seen in the below results. 
             (For *FIXED* scan, the entire radiometric image corresponds to a single pixel and hence pixel size is equal to the swath-width.)
   
   .. code-block:: python

      from instrupy.radiometer_model import RadiometerModel

      # instrument specifications as a dict
      radio1_json = {"@type": "Radiometer", "name": "ray1", "mass": 50, "volume": 3, "power": 10,
                     "orientation": {"referenceFrame": "SC_BODY_FIXED", "convention": "REF_FRAME_ALIGNED"},
                     "bitsPerPixel": 16,
                     "operatingFrequency": 1.25e9,
                     "antenna": {"shape": "CIRCULAR", "diameter": 1, "apertureExcitationProfile": "UNIFORM",
                                 "radiationEfficiency": 0.8, "phyTemp": 300},
                     "system": {"tlLoss": 0.5, "tlPhyTemp": 290, 
                                 "rfAmpGain": 30, "rfAmpInpNoiseTemp": 200, 
                                 "rfAmpGainVariation": 10, "mixerGain": 23, "mixerInpNoiseTemp": 1200,
                                 "mixerGainVariation": 2, "ifAmpGain": 30, "ifAmpInputNoiseTemp": 100,
                                 "ifAmpGainVariation": 10, "integratorVoltageGain": 1, "integrationTime": 100e-3,
                                 "bandwidth": 10e6, "@type": "TOTAL_POWER"},
                     "scan": {"@type": "FIXED"},
                     "targetBrightnessTemp": 345
                     }
      
      # define the viewing-geometry (observer state and target location)
      epoch_JDUT1 =  2458543.06088 # 2019 Feb 28 13:27:40 is time at which the ECEF and ECI frames approximately align, hence ECEF to ECI rotation is identity. See <https://www.celnav.de/longterm.htm> online calculator of GMST.
      SpacecraftOrbitState = {'time [JDUT1]':epoch_JDUT1, 'x [km]': 6878.137, 'y [km]': 0, 'z [km]': 0, 'vx [km/s]': 0, 'vy [km/s]': 7.6126, 'vz [km/s]': 0} # altitude 500 km
      TargetCoords = {'lat [deg]': 0, 'lon [deg]': 0.5} 

      radio1 = RadiometerModel.from_json(self.radio1_json)
      data_metrics = radio1.calc_data_metrics(sc_orbit_state=SpacecraftOrbitState, target_coords=TargetCoords, instru_look_angle_from_target_inc_angle=False)
      print(data_metrics)

      >> {'ground pixel along-track resolution [m]': 120708.29, 'ground pixel cross-track resolution [m]': 121567.92, 'swath-width [m]': 120565.56, 
          'sensitivity [K]': 17.94, 'incidence angle [deg]': 6.82, 'beam efficiency': nan}

   In the below snippet, the instrument look angle is considered from the incidence angle at the target ground-point (6.82 deg).
   Note that the calculated swath-width is different.
      
   .. code-block:: python

      data_metrics = o.calc_data_metrics(sc_orbit_state=SpacecraftOrbitState, target_coords=TargetCoords, instru_look_angle_from_target_inc_angle=True)
      print(data_metrics)

      >> {'ground pixel along-track resolution [m]': 120708.29, 'ground pixel cross-track resolution [m]': 121567.92, 'swath-width [m]': 122255.0, 
            'sensitivity [K]': 17.94, 'incidence angle [deg]': 6.82, 'beam efficiency': nan}



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