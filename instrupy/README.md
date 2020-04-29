# InstruPy

InstruPy is a python package to calculate observation data metrics for a given instrument and associated access events. 

For a detailed description see the following articles: 

1. V. Ravindra, S. Nag, *"Instrument Data Metrics Evaluator for Tradespace Analysis of Earth Observing Constellations",* IEEE Aerospace Conference, Big Sky, Montana, March 2020. 

2. S. Nag, V. Ravindra, J.J. LeMoigne *"Instrument Modeling Concepts for Tradespace Analysis of Satellite Constellations",* IEEE Sensors Conference, Delhi, India, October 2018.


## Install

Requires: `python 3.8`, `gfortran`

1. Navigate to the `instruments/instrupy/` directory and run `make`. 
2. Run tests using the `make runtest` command and get the *OK* message.

Find the documentation in: `instruments/instrupy/docs/_build/html/user_json_desc.html`

## Examples

* Specifications of example instruments (in the required JSON format) is present in the 
  `instrupy/examples/example_instrument_specs/` directory.

* The folder `instrupy/examples/` contains few examples generated upon execution of `instrupy/bin/instrument_module` which
  has been deprecated.

## Questions?

Please contact Vinay (vinay.ravindra@nasa.gov)

