# InstruPy

InstruPy is a python package to calculate observation data metrics for a given instrument and associated access events. 

For a detailed description see the following articles: 

1. V. Ravindra, S. Nag, *"Instrument Data Metrics Evaluator for Tradespace Analysis of Earth Observing Constellations",* IEEE Aerospace Conference, Big Sky, Montana, March 2020. 

2. S. Nag, V. Ravindra, J.J. LeMoigne *"Instrument Modeling Concepts for Tradespace Analysis of Satellite Constellations",* IEEE Sensors Conference, Delhi, India, October 2018.


## Install

Requires: `python 3.8`, `pip`, `gfortran`,  `lowtran`

1. Run `make` from the root repo directory.
2. Run tests using the `make runtest` command and get the *OK* message.

Find the documentation in: `instrupy/docs/_build/html/user_json_desc.html`

### Lowtran

Lowtran python package allows the execution of the LOWTRAN7 model. This is assumed to estimate the atmospheric losses for the visible and near-visible spectrum. 

The package is is available publicly here:
https://pypi.org/project/lowtran/

An backup copy is present in the `third_party\lowtran-2.4.1` folder.

This package requires the `gfortran` Fortran compiler. 

If a Fortran compiler is not already installed, install `gfortran` as follows:

* Linux: `apt install gfortran`

* Mac: `brew install gcc`

* Windows: https://www.scivision.dev/windows-gcc-gfortran-cmake-make-install/

## Examples

* Specifications of example instruments (in the required JSON format) is present in the 
  `instrupy/examples/example_instrument_specs/` directory.

* The folder `instrupy/examples/` contains few examples generated upon execution of `instrupy/bin/instrument_module` which
  has been deprecated.

## License and Copyright

Copyright 2020 Bay Area Environmental Research Institute

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

