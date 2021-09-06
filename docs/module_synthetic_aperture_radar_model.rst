.. _synthetic_aperture_radar_model_module:

``instrupy.synthetic_aperture_radar_model`` --- Synthetic Aperture Radar Model
******************************************************************************

Description
^^^^^^^^^^^^^

This module provides the ``SyntheticApertureRadarModel`` class which can be used to represent instruments imaging at optical/near-optical wavelengths.
The ``SyntheticApertureRadarModel`` object can be obtained from json/ dict representations using the ``from_json`` or ``from_dict`` functions. 
Please refer to :ref:`synthetic aperture radar model<synthetic_aperture_radar_model_desc>` for description of the instrument model and the expected key/value pairs.

Refer to the :ref:`base module<base_module>` and OrbitPy coverage calculator module (:class:`orbitpy.coveragecalculator`) to know how to use the ``SyntheticApertureRadarModel``
object for coverage calculations.

The ``calc_data_metrics`` function can be used to calculate the data-metrics given the viewing-geometry. There are two ways in which the 
viewing-geometry can be specified:

**1. With spacecraft-state and target location:**

      Here the spacecraft-state and observed-location (target) are specified.
      The expected key/value pairs for the spacecraft-state input argument (of type :class:`dict`) are:

      .. csv-table:: 
            :header: Column, Data type, Units, Description
            :widths: 30,10,10,40

            *time [JDUT1]*, float, Julian Day UT1, Time in Julian Day UT1. Corresponds to the time of observation. 
            "*x [km]*, *y [km]*, *z [km]*", float, km, Cartesian coordinates of satellite in EARTH_CENTERED_INERTIAL frame at the time of observation.
            "*vx [km/s]*, *vy [km/s]*, *vz [km/s]*", float, km/s, Velocity of spacecraft in EARTH_CENTERED_INERTIAL frame at the time of observation.

      The expected key/value pairs for the observed-location input argument (of type :class:`dict`) are:

      .. csv-table:: 
            :header: Column, Data type, Units, Description
            :widths: 30,10,10,40

            "*lat [deg]*, *lon [deg]*", float, degrees, "Ground-point accessed (geocentric latitude, geocentric longitude)."

**2. With incidence-angle specification:**

      Here the following input arguments is expected to the ``calc_data_metrics`` function:

      .. csv-table:: 
            :header: Column, Data type, Units, Description
            :widths: 30,10,10,40

            alt_km, float, km, Spacecraft altitude.
            sc_speed_kmps, float, km/s, Spacecraft speed in orbit.
            sc_gnd_speed_kmps, float, km/s, Spacecraft ground-speed/ antenna-footprint speed.
            inc_angle_deg, float, degrees, Incidence angle at the observed location.


``instru_look_angle_from_target_inc_angle`` flag
---------------------------------------------------

In addition to the above mentioned input-arguments to the ``calc_data_metrics`` function, there is an additional optional flag which can be 
specified to instruct how to obtain the instrument look-angle. If the ``instru_look_angle_from_target_inc_angle`` flag is set to ``False`` (also default value),
the instrument look-angle is obtained from the instrument orientation specification. If the flag is set to ``True``, the instrument look-angle is obtained from the target location incidence-angle. 

The instrument look-angle influences the calculation of the swath-width, (valid) PRF and hence the NESZ.

Examples
^^^^^^^^^

Example model parameters derived from real-world instruments can be found in the ``examples`` folder. Note that the models do not emulate the
real-world instrument entirely, but are meant to be approximate versions.

1. Example of MicroXSAR (X-band SAR). The default polarization (*SINGLE*), swath (*FULL*) and scan-technique (*Stripmap*) is considered.
   The view geometry is specified with the spacecraft state and ground-point geo-coordinates. 
   The instrument look angle considered is from the specified instrument orientation. 
   
      .. code-block:: python

            import numpy as np
            from instrupy.synthetic_aperture_radar_model import SyntheticApertureRadarModel

            RE = 6378.137
            def orbital_speed(alt_km):
            return np.sqrt(398600.5/(RE + alt_km))

            microxsar_dict = {"@type": "Synthetic Aperture Radar",
                              "name": "MiroXSAR",
                              "mass": 130,
                              "volume": 0.343,
                              "power": 1100,
                              "orientation": {
                                    "referenceFrame": "SC_BODY_FIXED",
                                    "convention": "SIDE_LOOK",
                                    "sideLookAngle": 30
                              },
                              "dataRate": 2000,
                              "bitsPerPixel": 16,
                              "pulseWidth": 31e-6,
                              "antenna":{"shape": "RECTANGULAR", "height": 4.9, "width": 0.7, 
                                    "apertureEfficiency": 0.5, 
                                    "apertureExcitationProfile": "UNIFORM"},
                              "operatingFrequency": 9.65e9,
                              "peakTransmitPower": 1e3,
                              "chirpBandwidth": 75e6,
                              "minimumPRF": 3000,
                              "maximumPRF": 8000,
                              "radarLoss": 3.5,
                              "sceneNoiseTemp": 290,
                              "systemNoiseFigure": 4.3
                              }
            microxsar = SyntheticApertureRadarModel.from_dict(microxsar_dict)
            epoch_JDUT1 = 2451623.999630# 2000 3 20 11 59 28.000
            # lat = 0, lon = 0 corresponds to [6378, 0, 0] km in ECI for observer position, check using Matlab function: eci2lla([6378, 0, 0] ,[2000 3 20 11 59 28.000])
            sc_orbit_state = {'time [JDUT1]':epoch_JDUT1, 'x [km]': 6378.137 + 600, 'y [km]': 0, 'z [km]': 0, 'vx [km/s]': 0, 'vy [km/s]': orbital_speed(600), 'vz [km/s]': 0} # equatorial orbit, altitude 600 km
            target_coords = {'lat [deg]': 2, 'lon [deg]': 0} # incidence angle to the target is 28.9855 deg
            obsv_metrics = microxsar.calc_data_metrics(sc_orbit_state, target_coords, instru_look_angle_from_target_inc_angle=False)
            print(obsv_metrics)

            >> {'ground pixel along-track resolution [m]': 2.09, 'ground pixel cross-track resolution [m]': 4.95, 
            'NESZ [dB]': -12.89, 'incidence angle [deg]': 28.99, 'swath-width [km]': 37.3, 'PRF [Hz]': 3751}

      In the below snippet, the instrument look angle is considered from the incidence angle at the target ground-point (28.99 deg).
      Note that the calculated swath-width is smaller. A different PRF is evaluated (higher PRF) which results in a different (improved) NESZ.
      
      .. code-block:: python

            obsv_metrics = microxsar.calc_data_metrics(sc_orbit_state, target_coords, instru_look_angle_from_target_inc_angle=True)
            print(obsv_metrics)

            >> {'ground pixel along-track resolution [m]': 2.09, 'ground pixel cross-track resolution [m]': 4.95, 
                'NESZ [dB]': -14.12, 'incidence angle [deg]': 28.99, 'swath-width [km]': 34.4, 'PRF [Hz]': 4976}

2. Example of dual-polarization SAR and *FIXED* swath configuration. View geometry is specified using the altitude, spacecraft speed, footprint speed and incidence angle specifications.
   The ``instru_look_angle_from_target_inc_angle`` flag is set to ``True``.

      .. code-block:: python

            import numpy as np
            from instrupy.synthetic_aperture_radar_model import SyntheticApertureRadarModel

            RE = 6378.137
            def orbital_speed(alt_km):
            return np.sqrt(398600.5/(RE + alt_km))

            alt_km = 500
            orb_speed = orbital_speed(alt_km)

            test_sar_dict = { "@type": "Synthetic Aperture Radar",
                              "name": "DSHIELD L-Band SAR",
                              "orientation": {
                                    "referenceFrame": "SC_BODY_FIXED",
                                    "convention": "SIDE_LOOK",
                                    "sideLookAngle": 45
                              },
                              "pulseWidth": 14.16e-6,
                              "antenna":{"shape": "RECTANGULAR", "height": 14.38, "width": 1.48, 
                                    "apertureEfficiency": 0.6, "apertureExcitationProfile": "UNIFORM"},
                              "operatingFrequency": 1280e6, 
                              "peakTransmitPower": 1000, 
                              "chirpBandwidth": 0.86e6,      
                              "minimumPRF": 1, 
                              "maximumPRF": 20000, 
                              "radarLoss": 2, 
                              "systemNoiseFigure": 2,
                              "swathConfig": {
                                    "@type": "fixed",
                                    "fixedSwathSize": 25
                              },
                              "polarization": {
                                    "@type": "dual",
                                    "pulseConfig": {
                                    "@type": "SMAP"   
                              } 
                              }                                                   
                        }
            test_sar = SyntheticApertureRadarModel.from_dict(test_sar_dict)

            inc_deg = 35
            obsv_metrics = test_sar.calc_data_metrics(alt_km = alt_km, sc_speed_kmps = orb_speed, sc_gnd_speed_kmps = orb_speed*(RE/(RE+alt_km)), inc_angle_deg = inc_deg, 
                                                      instru_look_angle_from_target_inc_angle = True)

            print(obsv_metrics)
            >> {'ground pixel along-track resolution [m]': 6.67, 'ground pixel cross-track resolution [m]': 364.66, 
                  'NESZ [dB]': -40.69, 'incidence angle [deg]': 35.0, 'swath-width [km]': 25.0, 'PRF [Hz]': 2666}


            inc_deg = 45
            obsv_metrics = test_sar.calc_data_metrics(alt_km = alt_km, sc_speed_kmps = orb_speed, sc_gnd_speed_kmps = orb_speed*(RE/(RE+alt_km)), inc_angle_deg = inc_deg, 
                                                      instru_look_angle_from_target_inc_angle = True)

            print(obsv_metrics)
            >> {'ground pixel along-track resolution [m]': 6.67, 'ground pixel cross-track resolution [m]': 295.79, 
                  'NESZ [dB]': -37.41, 'incidence angle [deg]': 45.0, 'swath-width [km]': 25.0, 'PRF [Hz]': 2279}

            inc_deg = 55
            obsv_metrics = test_sar.calc_data_metrics(alt_km = alt_km, sc_speed_kmps = orb_speed, sc_gnd_speed_kmps = orb_speed*(RE/(RE+alt_km)), inc_angle_deg = inc_deg, 
                                                      instru_look_angle_from_target_inc_angle = True)

            print(obsv_metrics)
            >> {'ground pixel along-track resolution [m]': 6.67, 'ground pixel cross-track resolution [m]': 255.33, 
                  'NESZ [dB]': -32.87, 'incidence angle [deg]': 55.0, 'swath-width [km]': 25.0, 'PRF [Hz]': 1578}

            








