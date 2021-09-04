.. _passive_optical_scanner_model_module:

``instrupy.passive_optical_scanner_model`` --- Passive Optical Scanner Model
******************************************************************************

Description
^^^^^^^^^^^^


  
Examples
^^^^^^^^^

Example model parameters derived from real-world instruments can be found in the ``examples`` folder. Note that the models do not emulate the
real-world instrument entirely, but are meant to be approximate versions.

1. Example from SMAD 3rd edition Chapter 9: Firesat (*WHISKBROOM*).
   
   .. code-block:: python

         from instrupy.passive_optical_scanner_model import PassiveOpticalScannerModel

         firesat_dict =  { "@type": "Passive Optical Scanner",
                           "name": "FireSat",  
                           "mass": 28, 
                           "volume": 0.12, 
                           "power": 32, 
                           "fieldOfViewGeometry": {
                              "shape": "RECTANGULAR",
                              "angleHeight": 0.628,
                              "angleWidth": 115.8
                           },
                           "scanTechnique": "WHISKBROOM",
                           "orientation": {
                              "referenceFrame": "SC_BODY_FIXED",
                              "convention": "REF_FRAME_ALIGNED"
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
                           "atmosLossModel": "LOWTRAN7",
                           "_references": {
                              "name": "Space Mission Analysis and Design",
                              "edition": 3,
                              "chapter": 9
                           }
                        }
         firesat = PassiveOpticalScannerModel.from_dict(firesat_dict)

2. Model based on the Landsat TIRS instrument (*PUSHBROOM*). Only Band1 is modeled. The target location is the ground-point at the satellite nadir.
   A sceneFOV geometry specification is included since the instrument has a narrow along-track FOV geometry (0.0081 deg). This shall allows for
   a smaller propagation time-step and hence faster coverage calculations (at the expense of some loss in accuracy).

   .. code-block:: python

         from instrupy.passive_optical_scanner_model import PassiveOpticalScannerModel

         landsat_tirs_band1_dict = {"@type": "Passive Optical Scanner",
                                    "name": "Landsat 8 TIRS Band1",
                                    "mass": 236,
                                    "volume": 0.261, 
                                    "power": 380, 
                                    "fieldOfViewGeometry": {
                                        "shape": "RECTANGULAR",
                                        "angleHeight": 0.0081,
                                        "angleWidth": 15
                                    },
                                    "sceneFieldOfViewGeometry": {
                                        "shape": "RECTANGULAR",
                                        "angleHeight": 2,
                                        "angleWidth": 15
                                    },
                                    "scanTechnique": "PUSHBROOM",
                                    "orientation": {
                                        "referenceFrame": "SENSOR_BODY_FIXED",
                                        "convention": "REF_FRAME_ALIGNED"
                                    },
                                    "dataRate":  384,
                                    "numberDetectorRows": 1,
                                    "numberDetectorCols": 1850,
                                    "detectorWidth": 25e-6,
                                    "focalLength": 0.178,
                                    "operatingWavelength": 10.9e-6,
                                    "bandwidth": 0.6e-6,
                                    "quantumEff": 0.025,
                                    "targetBlackBodyTemp": 290,
                                    "bitsPerPixel": 12,
                                    "opticsSysEff": 0.60 ,
                                    "numOfReadOutE":  20,
                                    "apertureDia":  0.1085366,
                                    "Fnum":  1.64,
                                    "maxDetectorExposureTime": 3.49e-3,
                                    "atmosLossModel": "LOWTRAN7",
                                    "_comments": ["Above is Total payload data-rate not just off the TIRS.",
                                                "numReadOutE is guessed."]
                                   }
        landsat_tirs_band1 = PassiveOpticalScannerModel.from_dict(landsat_tirs_band1_dict)
        # landsat 8 orbit at 10 Apr 2021 14:24:17.819 UTC            
        sc_orbit_state = {'time [JDUT1]':2459315.100208333,  'x [km]': -7012.215259847972,    'y [km]': 981.6284579029395,    'z [km]': 16.62328546479549, 
                                                            'vx [km/s]': 0.1664588472531363, 'vy [km/s]': 1.055747095699285, 'vz [km/s]': 7.426472416008381 }
        target_coords = {'lat [deg]': 0.01942147899019397 , 'lon [deg]': 117.1899962481559} # nadir position of satellite
        obsv_metrics = landsat_tirs_band1.calc_data_metrics(sc_orbit_state, target_coords)
        print(obsv_metrics)

        >> {'ground pixel along-track resolution [m]': 98.78, 'ground pixel cross-track resolution [m]': 98.92, 'SNR': 1507.48, 
            'dynamic range': 113645.23, 'noise-equivalent delta T [K]': 0.04162}

3. Model based on CCAM (*MATRIX_IMAGER*). Note that SNR Is 0 since the time of observation is during the night.
   
   CCAM is referenced from the following paper: E. Allthorpe-Mullis et al., Cubesat camera: *A low cost imaging system for cubesat platforms*, in 7th Interplanetary CubeSat Workshop, 2018.
   
   .. code-block:: python

         from instrupy.passive_optical_scanner_model import PassiveOpticalScannerModel
         ccam_blue_band_dict = {
                                "@type": "Passive Optical Scanner",
                                "name": "CCAM",
                                "fieldOfViewGeometry": {
                                    "shape": "RECTANGULAR",
                                    "angleHeight": 1.2,
                                    "angleWidth": 1.2
                                },
                                "scanTechnique": "MATRIX_IMAGER",
                                "numberDetectorRows": 2048,
                                "numberDetectorCols": 2048,
                                "detectorWidth": 5.5e-6,
                                "focalLength": 520e-3,
                                "operatingWavelength": 470e-9,
                                "bandwidth": 150e-9,
                                "quantumEff": 0.40,
                                "targetBlackBodyTemp": 290,
                                "opticsSysEff": 0.6,
                                "numOfReadOutE": 13,
                                "apertureDia": 94.6e-3,
                                "Fnum": 5.5,
                                "maxDetectorExposureTime": 678e-6,
                                "atmosLossModel": "LOWTRAN7"
                            }
         ccam_blue_band = PassiveOpticalScannerModel.from_dict(ccam_blue_band_dict)
         # Aqua orbit at 10 Apr 2021 15:07:56.800 UTC  (NIGHT time)                                                                          
         sc_orbit_state = {'time [JDUT1]':2459315.130520833,  'x [km]': -5054.315202286442,    'y [km]': -4878.491479401228,    'z [km]': 883.5310463297755, 
                                                            'vx [km/s]': -1.417318347731835, 'vy [km/s]': 0.1319708892386859, 'vz [km/s]': -7.367383505358474 }
         target_coords = {'lat [deg]': 7.127116160568699 , 'lon [deg]': 158.1924750010043} # nadir position of satellite
         obsv_metrics = ccam_blue_band.calc_data_metrics(sc_orbit_state, target_coords)
         print(obsv_metrics)

         >> {'ground pixel along-track resolution [m]': 7.43, 'ground pixel cross-track resolution [m]': 7.44, 
            'SNR': 0.0, 'dynamic range': 0.0, 'noise-equivalent delta T [K]': 2302356852773662.0}


API
^^^^^

.. rubric:: Classes

.. autosummary::
   :nosignatures:
   :toctree: generated/
   :template: classes_template.rst
   :recursive:

   instrupy.passive_optical_scanner_model.ScanTech
   instrupy.passive_optical_scanner_model.AtmosphericLossModel
   instrupy.passive_optical_scanner_model.PassiveOpticalScannerModel


