Package Overview
********************

InstruPy is a python package to calculate observation data metrics for a given instrument and associated access events. 

.. figure:: instrupy_block_diagram.png
    :scale: 75 %
    :align: center

    The high-level function of the InstruPy package is shown in the figure. The package ingests 
    access data and satellite states (time, position, velocity) information, and instrument
    specifications (for which the access data is generated), and outputs typical data metrics of observations
    made during the access.

For a detailed description see the following articles: 

1. V. Ravindra, S. Nag, *"Instrument Data Metrics Evaluator for Tradespace Analysis of Earth Observing Constellations",* IEEE Aerospace Conference, Big Sky, Montana, March 2020. 
2. S. Nag, V. Ravindra, J.J. LeMoigne *"Instrument Modeling Concepts for Tradespace Analysis of Satellite Constellations",* IEEE Sensors Conference, Delhi, India, October 2018.

Glossary of terms used in the package
======================================

* Instrument, payload and sensor are synonymous.

* Grid-point, ground-point are target-point synonymous.

* Pixels vs Detectors
      
      Pixels: Refer to ground pixels imaged. Dimensions vary according to imaging geometry.
      Detectors: Refer to physical detector elements on the imaging aperture.

* Access vs Coverage

      While access refers to a target falling under a sensor FOV, coverage includes an additional condition that the satellite
      should be able to be make an observation. 

* Satellite, spacecraft are synonymous.


Coding Conventions
===================

* Variables denoting physical quantities, unless otherwise indicated are always in S.I. units.