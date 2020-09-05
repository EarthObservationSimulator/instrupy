Common Input JSON Objects
**************************

This page contains description of the user-configurable JSON objects shared by the different instruments. 

.. _mode_json_obj:

:code:`mode` JSON object format
================================
Several modes (in a list) maybe specified within a single instrument. Each mode corresponds to a specific operating point. For example, 
consider a *Synthetic Aperture Radar* type instrument which operates in both L-band and C-band. Such an instrument is considered
to be made up of two *sub-sensors* with each sub-sensor of instrument type *Synthetic Aperture Radar* and operating at L-band
and C-band. A mode-identifier is to be specified by the user with which the subsensor can be referenced.

.. csv-table:: Input parameter description 
   :header: Parameter, Type, Units, Description
   :widths: 10,10,10,40

   @id, string,, Unique identifier
   @type, string,, For SAR specify Stripmap or ScanSar. For others leave empty.

The parameters outside the mode block are used as the common parameters for all the sub-sensors, while the parameters specified
within a mode list entry are specific to the particular sub-sensor.

Example: The example below is that of a *Basic Sensor* type instrument with two modes. The common parameters for both the modes
are outside the :code:`mode` block. The `NadirObservationMode` has a nadir viewing geometry with a 35 deg Conical FOV. The `SideObservationMode`
has a side looking geometry (both the side) at 25 deg conical FOV.
 
.. code-block:: python

               {        
                  "@type": "Basic Sensor",
                  "name": "Atom",
                  "@id": "senX",  
                  "mass": 28, 
                  "volume": 0.12, 
                  "power": 32, 
                  "bitsPerPixel": 8, 
                  "mode":[{
                        "@id": "NadirObservationMode"
                        "fieldOfView": {
                              "sensorGeometry": "CONICAL",
                              "fullConeAngle": 35
                        },
                        "orientation": {
                              "convention": "NADIR"
                        }      
                  },
                  {
                        "@id": "SideObservationMode"
                        "fieldOfView": {
                           "sensorGeometry": "CONICAL",
                           "fullConeAngle": 25
                        },
                        "orientation": {
                           "convention": "SIDE_LOOK",
                           "sideLookAngle": 30
                        },
                        "maneuverability":{
                           "@type": "YAW180"
                        }        
                        }
                  ]
               }


.. _orientation_json_obj:

:code:`orientation` JSON object format
========================================
At zero rotation, the instrument is aligned to the satellite body-frame which in turn is aligned to the Nadir-frame nominally. 
Thus the instrument orientation with respect to the satellite body-frame is the same as it's orientation with respect to the Nadir-frame 
in the scenario that the satellite body-frame is aligned to the Nadir-frame. It is also assumed that the instrument imaging axis is
along the instrument z-axis.


The definition of Nadir-frame is as follows:

*Nadir-frame*

* :math:`\bf X_{nadir}` axis: :math:`-({\bf Z_{nadir}} \times {\bf V})`, where :math:`\bf V` is the Velocity vector of satellite in EarthFixed frame)
* :math:`\bf Y_{nadir}` axis: :math:`({\bf Z_{nadir}} \times {\bf X_{nadir}})`
* :math:`\bf Z_{nadir}` axis: Aligned to Nadir vector (vector from Satellite to center of Earth in EarthFixed frame)

The first subfield of the :code:`orientation` JSON object is the :code:`convention` subfield.

.. csv-table:: Input parameter description 
   :header: Parameter, Type, Units, Description
   :widths: 10,10,10,40

   convention, string,, "Accepted values are *NADIR*, *SIDE_LOOK* or *XYZ*."

According to the specified :code:`convention`, other subfields materialize as follows:

1. :code:`"convention": "NADIR"`

If the orientation is aligned to the Nadir-frame:

Example:

.. code-block:: python

               "orientation": {
                                "convention": "NADIR"
                              }

2. :code:`"convention": "SIDE_LOOK"`

If the orientation is to be specified via a side-look-angle, the following subfields apply:

.. csv-table:: Input parameter description 
   :header: Parameter, Type, Units, Description
   :widths: 10,10,10,40

   sideLookAngle, number, degrees, Commonly called as nadir/ off-nadir angle. 

.. note:: A positive SIDE_LOOK corresponds to anti-clockwise rotation applied around the Nadir frame y-axis.

Example:

.. code-block:: python

               "orientation": {
                                "convention": "SIDE_LOOK",
                                "sideLookAngle":10
                              }

 
3. :code:`"convention": "XYZ"`

Here the orientation is to be specified via set of three rotation angles about the instrument primary axis (which is first aligned to the
Nadir-frame in nominal case). 
The order of (intrinsic) rotations is: (1) rotation about instrument X-axis, (2) rotation about instrument Y-axis and last 
(3) rotation about instrument Z-axis.

.. csv-table:: Input parameter description 
   :header: Parameter, Type, Units, Description
   :widths: 10,10,10,40

   xRotation, number, degrees, rotation about instrument X-axis
   yRotation, number, degrees, rotation about instrument Y-axis
   zRotation, number, degrees, rotation about instrument Z-axis

Example:

.. code-block:: python

               "orientation": {
                                "convention": "XYZ",
                                "xRotation":10,
                                "yRotation":20,
                                "zRotation":0
                              }

.. _fieldOfView_json_obj:

:code:`fieldOfView` JSON object format
========================================
The :code:`fieldOfView` can be specified in three ways, according to the parameter :code:`sensorGeometry` definition.

1. :code:`"sensorGeometry": "CONICAL"`

    .. csv-table:: Input parameter description 
        :header: Parameter, Type,Description
        :widths: 10,10,10,40

        fullConeAngle, number, degrees, Full cone angle of the instrument FOV. 

    Example:

    .. code-block:: python

                "fieldOfView": {
                                    "sensorGeometry": "CONICAL",
                                    "fullConeAngle":10
                                }

2. :code:`"sensorGeometry": "RECTANGULAR"`

    .. csv-table:: Input parameter description 
        :header: Parameter, Type, Units, Description
        :widths: 10,10,10,40

        alongTrackFieldOfView, number, degrees, (full) along-track fov. 
        crossTrackFieldOfView, number, degrees, (full) cross-track fov.

    Example:

    .. code-block:: python

                "fieldOfView": {
                                    "sensorGeometry": "RECTANGULAR",
                                    "alongTrackFieldOfView":10,
                                    "crossTrackFieldOfView":30
                                }

    .. warning:: The along-track FOV and cross-track FOV specs are assigned assuming the instrument is in nominal orientation, i.e. the instrument is aligned to nadir-frame.
                 If the instrument is rotated about the satellite body frame (by specifying non-zero orientation angles in the instrument json specs file), the actual along-track
                 and cross-track fovs simulated maybe different.

3. :code:`"sensorGeometry": "CUSTOM"`

    In this case the field-of-view is specified in terms of clock ,cone angles. The definition of the clock, cone angles is the 
    same as used in Orbit and Coverage module, i.e.

    Cone angles are angles measured from +Z sensor axis. If (:math:`xP`, :math:`yP`, :math:`zP`) is a unit vector describing a FOV point, then the 
    cone angle for the point is :math:`\pi/2 - \sin^{-1} zP`

    Clock angles are angles (right ascensions) measured anti-clockwise from the + X-axis (of instrument).  If :math:`xP`, :math:`yP`, :math:`zP` is a unit vector describing a FOV point, then the 
    cone angle for the point is :math:`atan2(y,x)`

    .. csv-table:: Input parameter description 
        :header: Parameter, Type, Units, Description
        :widths: 10,10,10,40

        customConeAnglesVector, string, degrees, array of cone angle (angle from Nadir vector) values separated by commas
        customClockAnglesVector, string, degrees, array of clock values separated by commas

    .. note:: The number of values in :code:`customConeAnglesVector` and :code:`customClockAnglesVector` should be the same (or) the number of 
              values in :code:`customConeAnglesVector` should be just one and no values in :code:`customClockAnglesVector`.


Example:

.. code-block:: python

               "fieldOfView": {
                                "sensorGeometry": "CUSTOM",
                                "customConeAnglesVector": [10,10,10,10],
                                "customClockAnglesVector": [30, 120, 180, 280]
                              }

.. _maneuverability_json_object:

:code:`maneuverability` JSON object
========================================
Total maneuverability of payload pointing (combining satellite and payload maneuverability). Four types of 
maneuverability are accepted: `Fixed`, `Cone`, `RollOnly`, `Yaw180`, `Yaw180Roll` and should be indicated in the 
:code:`@type` name, value pair. Please refer to :ref:`manuv_desc` for a complete description of the options.

1. :code:`"@type":"Fixed"`

This option indicates that the payload shall be fixed at it's nominal orientation (specified inside the :code:`instrument`
JSON object). There is no maneuverability.

Example:

.. code-block:: javascript
   
   "maneuverability":{
        "@type":"Fixed"
   }

2. :code:`"@type":"Cone"`

This option indicates that the payload pointing axis can be manuvered inside a conical region of full-cone angle as indicated
by the :code:`fullConeAngle` name, value pair. The axis of the cone is aligned to the nominal orientation of the instrument specified
in the :code:`instrument` JSON object.

.. csv-table:: Expected parameters
   :header: Parameter, Data type, Units, Description
   :widths: 10,10,5,40

   fullConeAngle, float, degrees, Full cone angle of the maneuverability conical region

Example:

.. code-block:: javascript
   
   "maneuverability":{
        "@type":"Cone",
        "fullConeAngle": 25
   }

3. :code:`"@type":"RollOnly"`

This option indicates that the payload can be manuevered only along the roll axis (about the satellite velocity vector in Inertial frame).
Such an option is expected for instruments which require a pure-side-looking target geometry.
At a :math:`roll = 0` deg, the payload shall point at the nominal orientation specified in the :code:`instrument` JSON object. 
The range of possible roll is indicated by the :code:`rollMin` and :code:`rollMax` name, value pairs.

.. csv-table:: Expected parameters
   :header: Parameter, Data type, Units, Description
   :widths: 10,10,5,40

   rollMin, float, degrees, minimum roll angle
   rollMax, float, degrees, maximum roll angle

Example:

.. code-block:: javascript
   
   "maneuverability":{
        "@type":"RollOnly",
        "rollMin": -5,
        "rollMax": 5
   }

4. :code:`"@type":"Yaw180"`

This option allows for a 180 deg manuver option about the yaw axis. 

Example:

.. code-block:: javascript
   
   "maneuverability":{
        "@type":"Yaw180"
   }

5. :code:`"@type":"Yaw180Roll"`

This option is similar to the :code:`RollOnly` option, but also includes 180 deg manuver option about the yaw axis. 
Such an option is expected for instruments which require a pure-side-looking target geometry.
At a :math:`roll = 0` deg, the payload shall point at the nominal orientation specified in the :code:`instrument` JSON object. 
The range of possible roll is indicated by the :code:`rollMin` and :code:`rollMax` name, value pairs.

.. csv-table:: Expected parameters
   :header: Parameter, Data type, Units, Description
   :widths: 10,10,5,40

   rollMin, float, degrees, minimum roll angle
   rollMax, float, degrees, maximum roll angle

Example:

.. code-block:: javascript
   
   "maneuverability":{
        "@type":"Yaw180Roll",
        "rollMin": -5,
        "rollMax": 5
   }