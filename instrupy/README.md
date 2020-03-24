# instrument_module

Python library package to calculate typical metrics of quality of observation data. 

## Installation

This project was designed for Python 3.7. It requires the following non-standard (to Anaconda) packages:
 - `lowtran`
To make the `instrupy` library visible to the Python interpreter, from the project root directory (containing `setup.py`), run:
```shell
pip install -e .
```
where the trailing period (`.`) indicates to install from the current directory.

## Documentation

Documentation of this package is managed using [Sphinx] 
(https://www.sphinx-doc.org/en/master/index.html).

A html version can be generated locally by following the below steps:
<ol>
<li> Install Python 3.8 (a dependency fpr Sphinx) </li>
<li> Install [sphinx package](https://www.sphinx-doc.org/en/master/usage/installation.html).</li>
<li> Add the Sphinx directory path to <code>PATH</code>. The Sphinx directory is located in subfolder named <code>Scripts</code> in the Python folder.  </li>
<li> Fork the instrument_module git repository into a local directory. </li>
<li> Navigate to <code>/docs</code> folder of the local repository. Run the command <code>make html </code>.</li>
<li> Navigate to <code>/docs/_build/html</code>. Open <code>index.html</code>.</li>
</ol>

## Quick run example

```shell
python bin/instrument_module.py example/arch-1
```
where `example/arch-1` is the path to the architecture directory from where the inputs/ outputs are read/ written. 

Outputs: (Shown inside  * *)
```
|-- bin/
|-- example/
    |-- arch-1/
        |-- arch.json
        |-- sat1_accessInfo.csv
        |-- sat1_accessInfo.json
        |-- sat2_accessInfo.csv
        |-- sat2_accessInfo.json
        |-- sat3_accessInfo.csv
        |-- sat3_accessInfo.json
        |-- *sat1_level0_data_metrics.csv*
        |-- *sat2_level0_data_metrics.csv*
        |-- *sat3_level0_data_metrics.csv*
        |-- *level1_data_metrics.csv*
        |-- *level2_data_metrics.csv*
        |-- *level1_coverage_metrics.csv*
        |-- *level2_coverage_metrics.csv*
|-- instrupy/
```
Index of examples:

* `arch-1` Basic Sensor type example
* `arch-2` Firesat, Thermal whiskbroom, 1 satellite, 20000 ground-points
* `arch-3` Landsat-8 TIRS pushbroom, Thermal (Band-1), 15 days, 5 satellites Homogeneous Walker, 20000 ground-points
* `arch-4` Landsat-8 OLI pushbroom, Visible (Blue Band), 5 days, 5 satellites Homogeneous Walker, 20000 ground-points
* `arch-4` MicroXSAR, X-Band, 5 days, 5 satellites Homogeneous Walker, 65000 Ground-points


The `arch.json` files in the respective `arch` directories contain only the instrument related json fields. 

## Testing

This project includes unit tests. To run using `nosetests` module, run from the command line:

```shell
python -m nose
```

## Producing OC understandable Coverage specifications

This package allows to extract OC understandable coverage specifications (a json string) from json input instrument specifications. 

Input: JSON formatted string of the instrument specifications (See `InstruPy` documentation for details).

Output: JSON formatted string containing the instrument coverage specifications.

example usage: 

```shell
python bin/instrupy_get_coverage_specs.py "{\"@type\": \"Basic Sensor\", \"name\": \"Atom\",\"acronym\":\"At\",\"mass\":10,\"volume\":12.45, \"dataRate\": 40, \"power\": 12, \"orientation\":{\"convention\": \"XYZ\",\"xRotation\":0,\"yRotation\":0,\"zRotation\":0}, \"fieldOfView\": {\"sensorGeometry\": \"CONICAL\", \"fullConeAngle\": 10 }}"

{"Orientation": {"eulerAngle1": 0.0, "eulerAngle2": 0.0, "eulerAngle3": 0.0, "eulerSeq1": 1, "eulerSeq2": 2, "eulerSeq3": 3}, "fieldOfView": {"geometry": "CONICAL", "coneAnglesVector": [5.0], "clockAnglesVector": null, "AlongTrackFov": 10.0, "CrossTrackFov": 10.0}}
```

*Note:* The character '\' needs to be prefixed before every double-quote in the json formatted string input (as tested in Windows).
    

### To Do: 
- [ ] For non-technical users who want access to documentation directly, integration with <code>readthedocs</code> platform can be done. 