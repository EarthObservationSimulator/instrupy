Run
=====

There are two steps involved in using the InstruPy package for calculation of data-metrics associated with an instrument and an access-event:

1. Initialize the relevant instrument model from a python dictionary or json-string. Description of the expected model-parameters from the user 
   can be obtained from the :ref:`Instrument Models<instrument_models_desc>` section. Specifications of example instruments (in the required JSON format)
   is available in the :code:`instrupy/examples/` directory.
2. Invoke the ``calc_data_metrics`` function from the model instance by passing the access-event information.
   The access-event information is most commonly the observation-time, state (position/velocity) of the observer and the location of the target. Another type of access-event specification
   involve specifying the incidence angle at the target, which is available for the synthetic aperture radar model.
   For more details please refer to the (:ref:`API Reference <api_reference>`) section of the respective modules. 

All the instrument-models have identical interface functions, i.e.:

1. The model-instances can be obtained by using the ``from_dict``/ ``from_json`` functions in which the model parameters are passed as dict/json.
2. The ``to_dict`` function can be used to obtain a dictionary representation of the model.
3. The ``calc_data_metrics`` function can be used to calculate data-metrics by passing the access-event information as the input argument.

**Example:**

*A passive optical scanner instrument model.*

.. code-block:: python
      
      from instrupy.passive_optical_scanner_model import PassiveOpticalScannerModel

      # instrument specifications in a python dictionary
      firesat_dict =  { "@type": "Passive Optical Scanner",
                        "name": "FireSat",
                        "mass": 28,
                        "volume": 0.12,
                        "power": 32,
                        "fieldOfViewGeometry": {
                          "shape": "RECTanGULAR",
                          "angleHeight": 0.628,
                          "angleWidth": 115.8
                        },
                        "scanTechnique": "WhiskBROOM",
                        "orientation": {
                          "referenceFrame": "SC_BODY_FIXED",
                          "convention": "SIDE_loOK",
                          "sideLookAngle": 0
                        },
                        "dataRate": 85,
                        "numberDetectorRows": 256,
                        "numberDetectorCols": 1,
                        "detectorWidth": 30e-6,
                        "focalLength": 0.7,
                        "operatingWavelength": 4.2e-6,
                        "bandwidth": 1.9e-6,
                        "quantumEff": 0.5,
                        "targetBlackBodyTemp": 290,
                        "bitsPerPixel": 8,
                        "opticsSysEff": 0.75,
                        "numOfReadOutE": 25,
                        "apertureDia": 0.26,
                        "Fnum": 2.7,
                        "atmosLossModel": "LOWTRAN7"}
      
      # obtain the instrument model instance from the python dict
      firesat = PassiveOpticalScannerModel.from_dict(firesat_dict)

      # define the access event (time of observation, position, velocity of observer and position of target)
      epoch_JDUT1 = 2451623.999630
      sc_orbit_state = {'time [JDUT1]':epoch_JDUT1, 'x [km]': 7078.137, 'y [km]': 0, 'z [km]': 0, 'vx [km/s]': 0, 'vy [km/s]': 7.5, 'vz [km/s]': 0} # equatorial orbit, altitude about 700 km
      target_coords = {'lat [deg]': 0, 'lon [deg]': 0} # lat = 0, lon = 0 corresponds to [6378, 0, 0] km in ECI for observer position, check using Matlab function: eci2lla([6378, 0, 0] ,[2000 3 20 11 59 28.000])
      
      # execute the data metrics calculation
      data_metrics = firesat.calc_data_metrics(sc_orbit_state, target_coords)
      print(data_metrics)

      >> {'ground pixel along-track resolution [m]': 31.34, 'ground pixel cross-track resolution [m]': 32.91, 
          'SNR': 148.15, 'dynamic range': 902.27, 'noise-equivalent delta T [K]': 0.30609}

``Instrument`` class
-----------------------
The :class:`instrupy.base.Instrument` class can be used to model a instrument with several operating *modes*. The ``Instrument`` class
has the same interface functions as that of the instrument-model classes (i.e. ``from_dict``/ ``from_json``, ``to_dict`` and ``calc_data_metrics``).

Each mode corresponds to a specific operating point and the *mode-only* specifications are stored as a list under the ``mode`` key.
A mode-identifier can be specified by the user with which the corresponding mode can be referenced. Also see :ref:`here<mode_json_obj>`.

.. note:: The ``Instrument`` class may be used without specifying any modes. This may be desirable if one wishes to avoid using different classes
          for each instrument-model. Since the ``Instrument`` class provides the same interface functions as the instrument-model classes, the
          usage is identical. Refer to :ref:`base module<base_module>` API reference.

**Example**

Consider a *Synthetic Aperture Radar* instrument which can operate in single-pol or dual-pol. Such an instrument is considered
to be made up of two modes and the polarization specifications are included in a list under the ``mode`` key. Such an instrument specification
effectively makes two synthetic aperture radar models with common parameters specified by the key/value pairs other then the ``mode`` key/value. 
One model shall use single-pol and the other model uses dual-pol. 

Note the 3dB difference in the NESZ results of the below example. 

.. code-block:: python

               from instrupy.base import Instrument

               specs = {   "@type": "Synthetic Aperture Radar",
                            "orientation": {
                                "referenceFrame": "SC_BODY_FIXED",
                                "convention": "SIDE_LOOK",
                                "sideLookAngle": 30
                            },
                            "mode":[{
                                "@id": "dual-pol",
                                "polarization": {
                                "@type": "dual",
                                "pulseConfig": {"@type": "AIRSAR"} 
                                }
                            },
                            {   
                                "@id": "single-pol",
                                "polarization": {
                                "@type": "single"
                                }
                            }
                            ],
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
                                "fixedSwathSize": 50
                            }            
                                                                              
                        }

                # obtain the instrument model instance from the python dict
                instru = Instrument.from_dict(specs)

                # define the access event (time of observation, position, velocity of observer and position of target)
                epoch_JDUT1 = 2451623.999630# 2000 3 20 11 59 28.000
                # lat = 0, lon = 0 corresponds to [6378, 0, 0] km in ECI for observer position, check using Matlab function: eci2lla([6378, 0, 0] ,[2000 3 20 11 59 28.000])
                sc_orbit_state = {'time [JDUT1]':epoch_JDUT1, 'x [km]': 6378.137 + 600, 'y [km]': 0, 'z [km]': 0, 'vx [km/s]': 0, 'vy [km/s]': 7.559, 'vz [km/s]': 0} # equatorial orbit, altitude 600 km
                target_coords = {'lat [deg]': 2, 'lon [deg]': 0} # incidence angle to the target is 28.9855 deg

                # execute the data metrics calculation
                data_metrics = instru.calc_data_metrics(mode_id="single-pol", sc_orbit_state=sc_orbit_state, target_coords=target_coords)
                print(data_metrics)

                >> {'ground pixel along-track resolution [m]': 6.13, 'ground pixel cross-track resolution [m]': 431.62, 'NESZ [dB]': -39.58, 
                    'incidence angle [deg]': 28.99, 'swath-width [km]': 50.0, 'PRF [Hz]': 2494}


                data_metrics = instru.calc_data_metrics(mode_id="dual-pol", sc_orbit_state=sc_orbit_state, target_coords=target_coords)
                print(data_metrics)

                >> {'ground pixel along-track resolution [m]': 6.13, 'ground pixel cross-track resolution [m]': 431.62, 'NESZ [dB]': -36.57, 
                    'incidence angle [deg]': 28.99, 'swath-width [km]': 50.0, 'PRF [Hz]': 2494}
 