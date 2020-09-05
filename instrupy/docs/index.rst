.. instrupy documentation master file, created by
   sphinx-quickstart on Sat Jan  5 10:39:21 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to InstruPy's documentation!
************************************

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

Install
========

Requires: :code:`python 3.8`, :code:`pip`, :code:`gfortran`

1. Navigate to the :code:`instruments/instrupy/` directory and run :code:`make`. 
2. Run tests using the :code:`make runtest` command and get the *OK* message.

Find the documentation in: :code:`instruments/instrupy/docs/_build/html/index.html`

Examples
=========

* Specifications of example instruments (in the required JSON format) is present in the 
  :code:`instrupy/examples/example_instrument_specs/` directory.

* The folder :code:`instrupy/examples/` contains few examples generated upon execution of :code:`instrupy/bin/instrument_module` which
  has been deprecated.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   instruments_description
   common_json_objects
   api_reference
   miscellaneous
   tests

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

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

Questions?
==========
Please contact Vinay (vinay.ravindra@nasa.gov)

