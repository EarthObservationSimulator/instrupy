# InstruPy

InstruPy is a python package to calculate observation data metrics for a given instrument and associated access events. 

For a detailed description see the following articles: 

1. V. Ravindra, S. Nag, *"Instrument Data Metrics Evaluator for Tradespace Analysis of Earth Observing Constellations",* IEEE Aerospace Conference, Big Sky, Montana, March 2020. 

2. S. Nag, V. Ravindra, J.J. LeMoigne *"Instrument Modeling Concepts for Tradespace Analysis of Satellite Constellations",* IEEE Sensors Conference, Delhi, India, October 2018.


## Install

Requires: Unix-like operating system (Linux (Ubuntu, CentOS...), Mac), `python 3.8`, `pip`, `gfortran`

The installation can be carried out in a `conda` environment using the below steps.

1. Install `gfortran`. See [Resource](https://fortran-lang.org/learn/os_setup/install_gfortran).

*   Linux: `sudo apt gfortran`
*   Mac: `brew install gcc`

2. Have `conda` installed using the [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or [Anaconda](https://www.anaconda.com/products/individual) distributions.

3. Create and activate a new conda environment with python. Install `pip` in the environment.
```
conda create --name py38env python=3.8
conda activate py38env
conda install pip
```

4. Run `make` from the root repo directory.

5. Run tests using the `make runtest` command and get the *OK* message.

Find the documentation in: `instrupy/docs/_build/html/user_json_desc.html`

All the dependencies are automatically installed. If any errors are encountered please check that the following dependencies are 
installed:

* `numpy`
* `pandas`
* `scipy`
* `lowtran` (requires gfortran)
* `sphinx`
* `sphinx_rtd_theme==0.5.2`
* `metpy`
* `netCDF4`
* `astropy`

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

## Directory structure

```
C:.
├───docs
├───examples (example specs)
│   ├───example_instrument_specs
│
├───instrupy (primary source files)
├───tests (test scripts)
├───third_party
│   └───lowtran-2.4.1
├───TBD (old files)
```

## Examples

Specifications of example instruments (in the required JSON format) is present in the `instrupy/examples/example_instrument_specs/` directory.

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

Please contact Vinay (vinay.ravindra@nasa.gov, vravindra@baeri.org)

