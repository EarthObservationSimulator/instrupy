
.. _base_module:

``instrupy.base`` --- Base
********************************

Description
============
This module provides the ``Instrument`` class which serves as a base class (and an alternate means) to work with the InstruPy package. 
Instrument models can be initialized using the ``from_dict`` / ``from_json`` functions, where the expected key/value pairs depend on the
specific instrument-model (see :ref:`here<instrument_models_desc>`). 

A ``Instrument`` object can be *added* to a :class:`orbitpy.util.Spacecraft` object and coverage calculations can be performed.
Refer to the OrbitPy coverage calculator module (:class:`orbitpy.coveragecalculator`) to know how to use the ``Instrument``
object for coverage calculations.

The ``calc_data_metrics`` function can be used to calculate the data-metrics given the details of the access-event (commonly the spacecraft-state and observed-location (target)).
For synthetic-aperture-radars the access-event specification could alternatively involve specifying the incidence angle at the target.
Please refer to the description of the respective instrument-model for additional details.

The ``Instrument`` class also provides the functionality to specify *modes*. Please refer to :ref:`mode_json_obj` for details. Each mode corresponds to effectively a different 
model with the respective mode-parameters.

Example
========

1. A ``Instrument`` object based on the *Basic Sensor* model. There is no mode specification, which implies that the instrument has one mode 
   with the instrument model parameterized according to the input specifications. ``mode_id`` is entered as ``None``.

   .. code-block:: python

      from instrupy.base import Instrument

      instru1 = Instrument.from_json('{"name": "Alpha", "mass":10, "volume":12.45, "dataRate": 40, "bitsPerPixel": 8, "power": 12, \
                                 "orientation": {"referenceFrame": "SC_BODY_FIXED", "convention": "REF_FRAME_ALIGNED"}, \
                                 "fieldOfViewGeometry": {"shape": "CIRCULAR", "diameter":2.5 }, \
                                 "sceneFieldOfViewGeometry": {"shape": "CIRCULAR", "diameter":5 }, \
                                 "maneuver":{"maneuverType": "CIRCULAR", "diameter":10}, \
                                 "pointingOption": [{"referenceFrame": "NADIR_POINTING", "convention": "XYZ", "xRotation":0, "yRotation":2.5, "zRotation":0}, \
                                                      {"referenceFrame": "NADIR_POINTING", "convention": "XYZ", "xRotation":0, "yRotation":-2.5, "zRotation":0}  \
                                                      ], \
                                 "numberDetectorRows":5, "numberDetectorCols":10, "@id":"bs1", "@type":"Basic Sensor" \
                                 }')

      # define access-event
      epoch_JDUT1 =  2458543.06088 # 2019 Feb 28 13:27:40 is time at which the ECEF and ECI frames approximately align, hence ECEF to ECI rotation is identity. See <https://www.celnav.de/longterm.htm> online calculator of GMST.
      SpacecraftOrbitState = {'time [JDUT1]':epoch_JDUT1, 'x [km]': 6878.137, 'y [km]': 0, 'z [km]': 0, 'vx [km/s]': 0, 'vy [km/s]': 7.6126, 'vz [km/s]': 0} # altitude 500 km
      TargetCoords = {'lat [deg]': 0, 'lon [deg]': 0}

      obsv_metrics = instru1.calc_data_metrics(mode_id=None, sc_orbit_state=SpacecraftOrbitState, target_coords=TargetCoords)

      print(obsv_metrics)

      >> {'observation range [km]': 500.0, 'look angle [deg]': 0.03, 'incidence angle [deg]': 0.03, 'solar zenith [deg]': 20.33}

2. Instrument with multiple modes
   
   Consider a *Synthetic Aperture Radar* instrument which can operate with either Stripmap or ScanSAR scanning techniques. Such an instrument is considered
   to be made up of two modes and the scan specifications are included in a list under the ``mode`` key. Such an instrument specification
   effectively makes two synthetic aperture radar models with common parameters specified by the key/value pairs other then the ``mode`` key/value. 
   One model shall use Stripmap scan and the other model uses ScanSAR scan. 

   .. code-block:: python

            from instrupy.base import Instrument

            specs = {   "@type": "Synthetic Aperture Radar",
                        "pulseWidth": 12e-6,
                        "antenna":{"shape": "RECTANGULAR", "height": 4.9, "width": 1, "apertureEfficiency": 0.6, "apertureExcitationProfile": "UNIFORM"},
                        "operatingFrequency": 9.65e9,
                        "peakTransmitPower": 1000,
                        "chirpBandwidth": 75e6,
                        "minimumPRF": 1,
                        "maximumPRF": 20000,
                        "radarLoss": 3.5,
                        "systemNoiseFigure": 4.5,
                        "polarization": {
                           "@type": "dual",
                           "pulseConfig": {
                              "@type": "SMAP",
                              "pulseSeparation": 2e-6
                        }
                        },
                        "atmosLoss": 1,
                        "mode":[{
                           "@id": "stripmap-mode",
                           "scanTechnique": "Stripmap"
                        },
                        {   
                           "@id": "scansar-mode",
                           "scanTechnique": "ScanSAR",
                           "numSubSwaths": 2
                        }
                        ]         
                                                                           
                        }

            # obtain the instrument model instance from the python dict
            instru = Instrument.from_dict(specs)

            # define the access event (time of observation, position, velocity of observer and position of target)
            inc_deg = 30
            RE = 6378
            h = 600
            orb_speed = 7.559
            data_metrics = instru.calc_data_metrics(mode_id="stripmap-mode", alt_km=h, sc_speed_kmps=orb_speed, sc_gnd_speed_kmps=orb_speed*(RE/(RE+h)), inc_angle_deg=inc_deg, 
                                                      instru_look_angle_from_target_inc_angle=True)
            print(data_metrics)
            >> {'ground pixel along-track resolution [m]': 2.24, 'ground pixel cross-track resolution [m]': 4.8, 'NESZ [dB]': -16.23, 
               'incidence angle [deg]': 30.0, 'swath-width [km]': 24.5, 'PRF [Hz]': 6269}


            data_metrics = instru.calc_data_metrics(mode_id="scansar-mode", alt_km=h, sc_speed_kmps=orb_speed, sc_gnd_speed_kmps=orb_speed*(RE/(RE+h)), inc_angle_deg=inc_deg, 
                                                      instru_look_angle_from_target_inc_angle=True)
            print(data_metrics)
            >> {'ground pixel along-track resolution [m]': 4.48, 'ground pixel cross-track resolution [m]': 4.8, 'NESZ [dB]': -16.03, 
                'incidence angle [deg]': 30.0, 'swath-width [km]': 49.1, 'PRF [Hz]': 5996}
               
API
======

.. rubric:: Classes

.. autosummary::
   :nosignatures:
   :toctree: generated/
   :template: classes_template.rst
   :recursive:

   instrupy.base.InstrumentModelFactory
   instrupy.base.Instrument
