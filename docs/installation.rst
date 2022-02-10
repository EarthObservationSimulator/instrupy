Installation
==============
Requires: Unix-like operating system (Linux (Ubuntu, CentOS...), Mac), :code:`python 3.8`, :code:`pip`, :code:`gfortran`

The installation can be carried out in a :code:`conda` environment using the below steps.

1. Install :code:`gfortran`. See `Resource <https://fortran-lang.org/learn/os_setup/install_gfortran>`_.

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


