# InstruPy

InstruPy is a python package to calculate observation data metrics for a given instrument and associated access events. 

For a detailed description see the following articles: 

1. V. Ravindra, S. Nag, *"Instrument Data Metrics Evaluator for Tradespace Analysis of Earth Observing Constellations",* IEEE Aerospace Conference, Big Sky, Montana, March 2020. 

2. S. Nag, V. Ravindra, J.J. LeMoigne *"Instrument Modeling Concepts for Tradespace Analysis of Satellite Constellations",* IEEE Sensors Conference, Delhi, India, October 2018.


## Install

Requires: `python 3.8`, `gfortran`, `pip`

1. Navigate to the `instruments/instrupy/` directory and run `make`. 
2. Run tests using the `make runtest` command and get the *OK* message.

Find the documentation in: `instruments/instrupy/docs/_build/html/user_json_desc.html`

## Examples

* Specifications of example instruments (in the required JSON format) is present in the 
  `instrupy/examples/example_instrument_specs/` directory.

* The folder `instrupy/examples/` contains few examples generated upon execution of `instrupy/bin/instrument_module` which
  has been deprecated.

## License and Copyright

Copyright 2021 Bay Area Environmental Research Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

## Credits and Acknowledgments

This work has been funded by grants from the National Aeronautics and Space Administration (NASA) Earth Science Technology Office (ESTO) through the Advanced Information Systems Technology (AIST) Program.

## Questions?

Please contact Vinay (vinay.ravindra@nasa.gov)

