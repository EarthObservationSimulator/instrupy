Common Input JSON Objects
**************************

This page contains description of the user-configurable JSON objects shared by the different instruments. 

.. _mode_json_obj:

:code:`mode` JSON object format
================================
Several modes (in a list) maybe specified within a single instrument. Each mode corresponds to a specific operating point. For example, 
consider a *Synthetic Aperture Radar* type instrument which operates in both L-band and C-band. Such an instrument is considered
to be made up of two modes with one mode operating at L-band and the other mode at C-band. 
A mode-identifier can be specified by the user with which the corresponding mode can be referenced.

.. csv-table:: Input parameter description 
   :header: Parameter, Type, Units, Description
   :widths: 10,10,10,40

   @id, string,, Unique identifier of mode.

The parameters outside the mode block are used as the common parameters for all the modes, while the parameters specified
within a mode list entry are specific to the particular mode.

Example: The example below is that of a *Basic Sensor* type instrument with two modes. The common parameters for both the modes
are outside the :code:`mode` block. The `NadirObservationMode` has a nadir orientation while the `SideObservationMode`
has an off-nadir orientation.
 
.. code-block:: python

               specs = '{        
                           "@type": "Basic Sensor",
                           "name": "Atom",
                           "@id": "senX",  
                           "mass": 28, 
                           "volume": 0.12, 
                           "power": 32, 
                           "bitsPerPixel": 8, 
                           "fieldOfViewGeometry": {
                                       "sensorGeometry": "CIRCULAR",
                                       "diameter": 35
                                 },
                           "mode":[{
                                    "@id": "NadirObservationMode",                            
                                    "orientation": {
                                          "referenceFrame": "SC_BODY_FIXED",
                                          "convention": "REF_FRAME_ALIGNED"
                                    }      
                                 },
                                 {
                                    "@id": "SideObservationMode",
                                    "orientation": {
                                       "referenceFrame": "SC_BODY_FIXED",
                                       "convention": "SIDE_LOOK",
                                       "sideLookAngle": 30
                                 }       
                                 }
                           ]
                        }'

               x = Instrument.from_json(specs) 

.. _orientation_json_obj:

:code:`orientation` JSON object format
========================================
Orientation is parameterized as intrinsic rotations specified by Euler angles and sequence with respect to 
an user-specified reference frame. The definition of the Euler angle rotation is identical to the 
one used in the orbitpy->propcov->extern->gmatutil->util->AttitudeUtil, AttitudeConversionUtility C++ classes. 

A Euler sequence = 123 implies the following rotation: R = R3.R2.R1, where Ri is the rotation matrix about the ith axis.
A positive angle corresponds to an anti-clockwise rotation about the respective axis. Each rotation matrix rotates the 
coordinate system (not the vector).
See:

* https://mathworld.wolfram.com/RotationMatrix.html

The first subfield of the :code:`orientation` JSON object is the :code:`referenceFrame` subfield. 
See :ref:`reference_frames_desc` for description about the reference frames.

.. csv-table:: Input parameter description 
   :header: Parameter, Type, Units, Description
   :widths: 10,10,10,40

   referenceFrame, string,, "Accepted values are *EARTH_CENTERED_INERTIAL*, *EARTH_FIXED*, *NADIR_POINTING* or *SC_BODY_FIXED*."

The second subfield of the :code:`orientation` JSON object is the :code:`convention` subfield.

.. csv-table:: Input parameter description 
   :header: Parameter, Type, Units, Description
   :widths: 10,10,10,40

   convention, string,, "Accepted values are *REF_FRAME_ALIGNED*, *SIDE_LOOK*, *XYZ* or *EULER*."

According to the specified :code:`convention`, other subfields materialize as follows:

1. :code:`"convention": "REF_FRAME_ALIGNED"`

Aligned with respective to the underlying reference frame.

Example:

.. code-block:: python

               "orientation": {
                                "referenceFrame": "NADIR_POINTING",
                                "convention": "REF_FRAME_ALIGNED"
                              }

2. :code:`"convention": "SIDE_LOOK"`

If the orientation is to be specified via a side-look-angle (which corresponds to rotation about the y-axis only), the following subfields apply:

.. csv-table:: Input parameter description 
   :header: Parameter, Type, Units, Description
   :widths: 10,10,10,40

   sideLookAngle, float, degrees, Side-look angle

Example:

.. code-block:: python

               "orientation": {
                                "referenceFrame": "NADIR_POINTING",
                                "convention": "SIDE_LOOK",
                                "sideLookAngle":10
                              }

 
3. :code:`"convention": "XYZ"`

Here the orientation is to be specified via set of three rotation angles about the X, Y and Z axis.
The order of (intrinsic) rotations is: (1) rotation about instrument X-axis, (2) rotation about instrument Y-axis and last 
(3) rotation about instrument Z-axis.

.. csv-table:: Input parameter description 
   :header: Parameter, Type, Units, Description
   :widths: 10,10,10,40

   xRotation, float, degrees, rotation about instrument X-axis
   yRotation, float, degrees, rotation about instrument Y-axis
   zRotation, float, degrees, rotation about instrument Z-axis

Example:

.. code-block:: python

               "orientation": {
                                "referenceFrame": "NADIR_POINTING",
                                "convention": "XYZ",
                                "xRotation":10,
                                "yRotation":20,
                                "zRotation":0
                              }

4. :code:`"convention": "EULER"`

Here the orientation is to be specified via set of Euler angles and sequence.

.. csv-table:: Input parameter description 
   :header: Parameter, Type, Units, Description
   :widths: 10,10,10,40

   eulerAngle1, float, degrees, Rotation angle corresponding to the first rotation.
   eulerAngle2, float, degrees, Rotation angle corresponding to the second rotation.
   eulerAngle3, float, degrees, Rotation angle corresponding to the third rotation.
   eulerSeq1, int, Axis-number corresponding to the first rotation.
   eulerSeq2, int, Axis-number corresponding to the second rotation.
   eulerSeq3, int, Axis-number corresponding to the third rotation.

Example:

.. code-block:: python

               "orientation": {
                                "referenceFrame": "NADIR_POINTING",
                                "convention": "EULER",
                                "eulerAngle1":10,
                                "eulerAngle2":20,
                                "eulerAngle3":0,
                                "eulerSeq1": 3,
                                "eulerSeq2": 1,
                                "eulerSeq3": 3
                              }

.. _fieldOfViewGeometry_json_obj:

:code:`fieldOfViewGeometry` JSON object format
========================================
The :code:`fieldOfViewGeometry` is characterized by the key :code:`shape` definition. Three values are allows :code:`"CIRCULAR"`, :code:`RECTANGULAR`
and :code:`CUSTOM`.

1. :code:`"shape": "CIRCULAR"`

    .. csv-table:: Input parameter description 
        :header: Parameter, Type,Description
        :widths: 10,10,10,40

        diameter, number, degrees, Diameter (2 times the cone angle)

    Example:

    .. code-block:: python

                "fieldOfViewGeometry": {
                                          "shape": "CIRCULAR",
                                          "diameter":10
                                       }

2. :code:`"shape": "RECTANGULAR"`

    .. csv-table:: Input parameter description 
        :header: Parameter, Type, Units, Description
        :widths: 10,10,10,40

        angleHeight, number, degrees, Angular height (about sensor X-axis)
        angleWidth, number, degrees, Angular width (about sensor Y-axis)
    
    angleHeight and angleWith correspond to the along-track and cross-track FOVs respectively in case the sensor-frame is
    aligned to the NADIR_POINTING frame.

    Example:

    .. code-block:: python

                "fieldOfViewGeometry": {
                                          "shape": "RECTANGULAR",
                                          "angleHeight":10,
                                          "angleWidth":30
                                       }

3. :code:`"shape": "CUSTOM"`

    In this case the field-of-view is specified in terms of clock, cone angles. The definition of the clock, cone angles is the 
    same as used in Orbit and Coverage module, i.e.

    Cone angles are angles measured from +Z sensor axis. If (:math:`xP`, :math:`yP`, :math:`zP`) is a unit vector describing a FOV point, then the 
    cone angle for the point is :math:`\pi/2 - \sin^{-1} zP`

    Clock angles are angles (right ascensions) measured anti-clockwise from the + X-axis (of instrument).  If :math:`xP`, :math:`yP`, :math:`zP` is a unit vector describing a FOV point, then the 
    cone angle for the point is :math:`atan2(y,x)`

    .. csv-table:: Input parameter description 
        :header: Parameter, Type, Units, Description
        :widths: 10,10,10,40

        customConeAnglesVector, string, degrees, array of cone angle values separated by commas
        customClockAnglesVector, string, degrees, array of clock values separated by commas

    .. note:: The number of values in :code:`customConeAnglesVector` and :code:`customClockAnglesVector` should be the same (or) the number of 
              values in :code:`customConeAnglesVector` should be just one and no values in :code:`customClockAnglesVector`.


Example:

.. code-block:: python

               "fieldOfViewGeometry": {
                                          "shape": "CUSTOM",
                                          "customConeAnglesVector": [10,10,10,10],
                                          "customClockAnglesVector": [30,120,180,280]
                                       }

.. _maneuver_json_object:

:code:`maneuver` JSON object
========================================
Total maneuverability of sensor pointing (combining satellite and sensor maneuverability). Three types of 
maneuvers are accepted: `Circular`, `Single_Roll_Only` and `Double_Roll_Only`. This should be indicated in the 
:code:`maneuverType` name, value pair. Please refer to :ref:`maneuv_desc` for a complete description of the options.

1. :code:`"maneuverType":"Circular"`

This option indicates that the instrument pointing axis can be maneuvered about the nadir vector inside a circular region of diameter as indicated
by the :code:`diameter` name, value pair.

.. csv-table:: Expected parameters
   :header: Parameter, Data type, Units, Description
   :widths: 10,10,5,40

   diameter, float, degrees, Diameter

Example:

.. code-block:: python
   
   "maneuver":{
        "maneuverType":"Circular",
        "diameter": 25
   }

2. :code:`"maneuverType":"Single_Roll_Only"`

This option indicates that the instrument can be maneuvered only about the roll axis (of the nadir-pointing frame).
Such an option is expected for instruments which require a pure-side-looking target geometry.
The range of possible roll is indicated by the :code:`rollMin` and :code:`rollMax` name, value pairs. Note that these angles are
defined with respect to the NADIR_POINTING frame.

.. csv-table:: Expected parameters
   :header: Parameter, Data type, Units, Description
   :widths: 10,10,5,40

   A_rollMin, float, degrees, minimum roll angle
   A_rollMax, float, degrees, maximum roll angle

Example:

.. code-block:: python
   
   "maneuver":{
        "maneuverType":"Single_Roll_Only",
        "A_rollMin": 5,
        "A_rollMax": 15
   }

3. :code:`"maneuverType":"Double_Roll_Only"`

This option is similar to the :code:`Single_Roll_Only` option, except that it allows for definition of two set of roll-ranges (labelled as A and B).
This option is useful to model manuever by purely side-looking (look at the nadir is prohibited) instruments which may be pointed on either 'side' (i.e. positive roll region
and the negative roll region) of the nadir-pointing frame. 

.. csv-table:: Expected parameters
   :header: Parameter, Data type, Units, Description
   :widths: 10,10,5,40

   A_rollMin, float, degrees, minimum roll angle of roll region A
   A_rollMax, float, degrees, maximum roll angle of roll region A
   B_rollMin, float, degrees, minimum roll angle of roll region B
   B_rollMax, float, degrees, maximum roll angle of roll region B

Example:

.. code-block:: python
   
   "maneuver":{
        "maneuverType":"Double_Roll_Only",
        "A_rollMin": 5,
        "A_rollMax": 15,
        "B_rollMin": -15,
        "B_rollMax": -5
   }

.. _pointing_opt_json_obj:

:code:`pointingOption` JSON object
========================================
List of orientations to which the instrument axis can be manuevered. Only the NADIR_POINTING reference frame is supported.
This input specification is required to perform coverage calculations involving pointing-options.

Example:

.. code-block:: python
   
   "pointingOption":[{
      "referenceFrame": "NADIR_POINTING",
      "convention": "XYZ",
      "xRotation":0,
      "yRotation":20,
      "zRotation":0
   },
   {
      "referenceFrame": "NADIR_POINTING",
      "convention": "XYZ",
      "xRotation":0,
      "yRotation":40,
      "zRotation":0
   }]

.. _syntheticDataConfig_json_obj:

:code:`syntheticDataConfig` JSON object
================================================

This JSON object is used to describe the configuration of the synthetic data to be produced by the instrument models. 

.. csv-table:: Expected parameters
   :header: Parameter, Data type, Units, Description
   :widths: 10,10,5,40

   sourceFilePaths, list str,, List of absolute filepaths of the science-data files in NetCDF format. Each file corresponds to a specific (forecast/analysis) time.
   geophysicalVar, str,, Geophysical variable (name as present in the source NetCDF file) to be used for the synthetic data.
   interpolMethod, str,, Interpolation method to be employed while interpolating the source data onto the pixel-positions.

Example:

.. code-block:: python
   
   "syntheticDataConfig":{
        "sourceFilePaths": ["C:/workspace/gfs_forecast_data/gfs.t12z.pgrb2.0p25.f000.nc", 
                            "C:/workspace/gfs_forecast_data/gfs.t12z.pgrb2.0p25.f001.nc",
                            "C:/workspace/gfs_forecast_data/gfs.t12z.pgrb2.0p25.f002.nc",
                            "C:/workspace/gfs_forecast_data/gfs.t12z.pgrb2.0p25.f003.nc,
                            "C:/workspace/gfs_forecast_data/gfs.t12z.pgrb2.0p25.f004.nc"],
        "geophysicalVar": "TMP_P0_L1_GLL0",
        "interpolMethod": "SCIPY_LINEAR"
   }