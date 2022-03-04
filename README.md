# InstruPy

InstruPy is a python package to calculate observation data metrics for a given instrument and associated viewing geometry. 

For a detailed description see the following articles (available in the `literature` folder): 

1. V. Ravindra, S. Nag, *"Instrument Data Metrics Evaluator for Tradespace Analysis of Earth Observing Constellations",* IEEE Aerospace Conference, Big Sky, Montana, March 2020. 

2. V. Ravindra, R. Ketzner and S. Nag, *"Earth Observation Simulator (EO-Sim): An Open-Source Software for Observation Systems Design,"* 2021 IEEE International Geoscience and Remote Sensing Symposium IGARSS, 2021.

3. S. Nag, V. Ravindra, J.J. LeMoigne *"Instrument Modeling Concepts for Tradespace Analysis of Satellite Constellations",* IEEE Sensors Conference, Delhi, India, October 2018.

## Install

Requires: Unix-like operating system (Linux (Ubuntu, CentOS...), Mac), `python 3.8`, `pip`, `gfortran`, `make`

The installation can be carried out in a `conda` environment using the below steps.

1. Install `gfortran`. and `make` See [here](https://fortran-lang.org/learn/os_setup/install_gfortran).

*   Linux: 

    ```
    sudo apt update
    sudo apt install build-essential
    ```

    `sudo apt install build-essential` installs additional useful utilities like `gcc`, `g++` and `make`

    If `gfortran` is not installed with `build-essential`, please run `sudo apt install gfortran`

*   Mac: `brew install gcc` and `brew install make`

2. Have `conda` installed using the [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or [Anaconda](https://www.anaconda.com/products/individual) distributions.

3. Create and activate a new conda environment with python. Install `pip` in the environment.
```
conda create --name foo python=3.8
conda activate foo
conda install pip
```

4. Run `make` from the root repo directory.

    All the dependencies are automatically installed. If any errors are encountered please check that the following dependencies are 
    installed correctly.

    * `numpy`
    * `pandas`
    * `scipy`
    * `lowtran==2.4.1` (requires `gfortran`)
    * `sphinx`
    * `sphinx_rtd_theme==0.5.2`
    * `metpy`
    * `netCDF4`
    * `astropy`

5. Run tests using the `make runtest` command and get the *OK* message.

6. Find the documentation in: `instrupy/docs/_build/html/index.html`

The present version of OrbitPy has been tested on Ubuntu 18.04.3.

---

If using a Windows system, one may consider:

1. Setting up a virtual-machine with Ubuntu (or similar) OS. 

    See tutorial here: [https://ubuntu.com/tutorials/how-to-run-ubuntu-desktop-on-a-virtual-machine-using-virtualbox#1-overview](https://ubuntu.com/tutorials/how-to-run-ubuntu-desktop-on-a-virtual-machine-using-virtualbox#1-overview)

    Virtual Machines:

        *   VMware Workstation player is available free for non-commercial, personal or home use. VMWare tools may need to be installed separately  after the player installation. 
        [https://www.vmware.com/products/workstation-player/workstation-player-evaluation.html](https://www.vmware.com/products/workstation-player/workstation-player-evaluation.html)
    
        *   Another option is Oracle Virtual Box.
            [https://www.virtualbox.org/](https://www.virtualbox.org/)

2. Using Windows Subsytem for Linux (WSL)

    [https://ubuntu.com/wsl](https://ubuntu.com/wsl)
    
    [https://docs.microsoft.com/en-us/windows/wsl/about](https://docs.microsoft.com/en-us/windows/wsl/about)


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
├───literature (resources)
├───TBD (old files)
```

## Examples

InstruPy contains models of 'basic', passive-optical, synthetic-aperture-radar and radiometer instruments. Each type of instrument is associated
with different set of instrument parameters whose description can be found in the HTML documentation. 
Example specifications and related literature of instruments (in the required JSON format) is present in the `instrupy/examples/` directory.

The directory also contains the following python scripts which can be executed after the InstruPy package has been installed.

* SAR_example.py: Illustrates the synthetic-aperture-radar models with different possible set of configurations and the usage of the InstruPy functions to      evaluate the data-metrics.

## License and Copyright

Copyright 2022 Bay Area Environmental Research Institute

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

Please contact Vinay (vravindra@baeri.org or vinay.ravindra@nasa.gov)

