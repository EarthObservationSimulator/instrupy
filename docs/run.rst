Run
=====

There are two steps involved in using the InstruPy package for calculation of data-metrics associated with an instrument and an access-event:

1. Initialize the relevant instrument model from a python dictionary or json-string. Description of the expected model-parameters from the user 
   can be obtained from the :ref:`here<instrument_models_desc>`. 
2. Invoke the ``calc_data_metrics`` function from the model instance by passing the access-event information (state of the observer and target).
   The access-event information is most commonly the state of the observer and the location of the target. Another type of access-event specification
   involve specifying the incidence angle at the target, which is available for the synthetic aperture radar model.
   For more details please refer to the API description of the respective modules (:ref:`here<api_reference>`). 

**Example:**

*A passive optical scanner instrument model.*

.. code-block:: python
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
      obsv_metrics = firesat.calc_data_metrics(sc_orbit_state, target_coords)

The :class:`instrupy.base.Instrument` class can be


Specifications of example instruments (in the required JSON format) is present in the :code:`instrupy/examples/` directory.
