Package Overview
********************

InstruPy is a python package to calculate (satellite) observation data metrics for a given instrument and associated viewing geometry (access events). 

.. figure:: instrupy_block_diagram.png
    :scale: 75 %
    :align: center

    The high-level function of the InstruPy package is shown in the figure. The package ingests 
    access data i.e. the remote-sensing satellite states (time, position, velocity) information, and instrument
    specifications (for which the access data was generated), and outputs data metrics of observation.

For a detailed description see the following articles: 

1. V. Ravindra, R. Ketzner, S. Nag, *"Earth Observation Simulator (EO-SIM): An Open-Source Software for Observation Systems Design",* IEEE International Geoscience and Remote Sensing Symposium, Brussels Belgium, July 2021.
2. V. Ravindra, S. Nag, *"Instrument Data Metrics Evaluator for Tradespace Analysis of Earth Observing Constellations",* IEEE Aerospace Conference, Big Sky, Montana, March 2020. 
3. S. Nag, V. Ravindra, J.J. LeMoigne *"Instrument Modeling Concepts for Tradespace Analysis of Satellite Constellations",* IEEE Sensors Conference, Delhi, India, October 2018.

Currently there are four instrument models supported:

1. Basic Sensor (:ref:`basic_sensor_model_desc`, :ref:`basic_sensor_model_module`)
2. Passive Optical Sensor (:ref:`passive_optical_scanner_model_desc`, :ref:`passive_optical_scanner_model_module`)
3. Synthetic Aperture Radar (:ref:`synthetic_aperture_radar_model_desc`, :ref:`synthetic_aperture_radar_model_module`)
4. Radiometer (:ref:`radiometer_model_desc`, :ref:`radiometer_model_module`)

Glossary of terms used in the package
======================================

* Satellite and spacecraft are synonymous.
  
* Instrument, payload and sensor are synonymous.

* Grid-point, ground-point and target-point synonymous.

* Pixels vs Detectors
      
      *Pixels:* Refer to ground pixels imaged. Dimensions vary according to imaging geometry.
      *Detectors:* Refer to physical detector elements on the imaging aperture.

* Access vs Coverage

      * Sometimes access and coverage are used synonymously.

      * Other times access refers to a target falling under a sensor FOV while coverage includes an additional condition that the satellite
        should be able to be make an observation. For example for an optical sensor a ground-point may fall under the sensor FOV 
        (which is regarded as an access event), but may not be observable because it is night-time (no coverage).

Coding Conventions
===================

* Variables denoting physical quantities, unless otherwise indicated are always in S.I. units.