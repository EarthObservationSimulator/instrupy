.. _reflectometer_model_desc:

Reflectometer Model
********************

The reflectometer model is a simple model used to represent instruments operating on the concept of reflectometery, i.e. instruments which
receive and process reflected signals of opportunity (e.g. GNSS signals) of the Earths surface. A circular microwave antenna and its orientation 
is required to be specified based on which the field-of-view is defined. Pointing-options can also be defined to evaluate the instrument at different orientations.
Alternatively the maneuver can be defined for continuous range of instrument orientations.

The sceneFOV/ FOR is used in the coverage calculations (using the OrbitPy package) to find the locations accessed (specular locations) on the ground.


