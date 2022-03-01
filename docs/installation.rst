Installation
==============
Requires: Unix-like operating system (Linux (Ubuntu, CentOS...), Mac), :code:`python 3.8`, :code:`pip`, :code:`gfortran`

The installation can be carried out in a :code:`conda` environment using the below steps.

1. Install :code:`gfortran`. See `here <https://fortran-lang.org/learn/os_setup/install_gfortran>`_.

*   Linux: :code:`sudo apt gfortran`
*   Mac: :code:`brew install gcc`

2. Have :code:`conda` installed using the `Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_ or `Anaconda <https://www.anaconda.com/products/individual>`_ distributions.

3. Create and activate a new conda environment with python. Install :code:`pip` in the environment.

    .. code-block:: bash

                    conda create --name foo python=3.8
                    conda activate foo
                    conda install pip


4. Run :code:`make` from the root repo directory.

    All the dependencies are automatically installed. If any errors are encountered please check that the following dependencies are 
    installed:

    * :code:`numpy`
    * :code:`pandas`
    * :code:`scipy`
    * :code:`lowtran==2.4.1` (requires :code:`gfortran`)
    * :code:`sphinx`
    * :code:`sphinx_rtd_theme==0.5.2`
    * :code:`metpy`
    * :code:`netCDF4`
    * :code:`astropy`

5. Run tests using the :code:`make runtest` command and get the *OK* message.

6. Find the documentation in: :code:`instrupy/docs/_build/html/user_json_desc.html`

Windows
--------

If using a Windows system, one may consider:

1. Using Windows Subsytem for Linux (WSL)

    `https://ubuntu.com/wsl <https://ubuntu.com/wsl>`_
    
    `https://docs.microsoft.com/en-us/windows/wsl/about <https://docs.microsoft.com/en-us/windows/wsl/about>`_

2. Setting up a virtual-machine with Ubuntu (or similar) OS. 

    `VMware Workstation player <https://www.vmware.com/products/workstation-player/workstation-player-evaluation.html>`_  is available free for non-commercial, personal or home use. VMWare tools may need to be installed separately after the player installation.
    
    Another option is `Oracle Virtual Box <https://www.virtualbox.org/>`_.

The present version of OrbitPy has been tested on Ubuntu 18.04.3.