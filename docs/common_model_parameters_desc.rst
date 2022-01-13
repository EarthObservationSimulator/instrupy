.. _common_instru_params:

Common Model Parameters
*************************
This section describes the common instrument model parameters expected to be specified by the user.
The expected key/value pairs for initialization of the parameters from json/ dict is described along with examples. 

.. _orientation_json_obj:

:code:`orientation` json object 
================================

Orientation (of a instrument or spacecraft) is parameterized as intrinsic rotations specified by Euler angles and sequence with respect to 
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

Note that the rotation about y-axis corresponds to a roll motion when the reference-frame is the nadir-pointing frame.

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
   eulerSeq1, int, , Axis-number corresponding to the first rotation.
   eulerSeq2, int, , Axis-number corresponding to the second rotation.
   eulerSeq3, int, , Axis-number corresponding to the third rotation.

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

:code:`fieldOfViewGeometry` (Spherical-Geometry) json object 
=============================================================
The :code:`fieldOfViewGeometry` is used to characterize the spherical-geometry of the sensor field-of-view/ scene-field-of-view/ field of regard
in the *SENSOR_BODY_FIXED* frame. The Z-axis is assumed to be the pointing-axis.
Note that the orientation of the sensor is required to complete the field-of-view (/ scene-field-of-view/ field of regard) definition. 

The ``fieldOfViewGeometry`` json object is characterized by the key :code:`shape` definition. 
Three values are allows :code:`"CIRCULAR"`, :code:`RECTANGULAR` and :code:`CUSTOM`.

1. :code:`"shape": "CIRCULAR"`

   Specifies a circular shape about the sensor Z-axis.

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

   Specifies a rectangular shape about the sensor Z-axis.

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

   In this case the field-of-view geometry is specified in terms of list of vertices of a spherical polygon. Each vertex is specified via its clock, cone angles. 
   The definition of the clock, cone angles is the same as used in OrbitPy (propcov) package, i.e.

   Let (:math:`x_P`, :math:`y_P`, :math:`z_P`) be a unit vector describing a point on the unit sphere.
   
   The cone angle for the point is:

   :math:`\pi/2 - \sin^{-1}z_P`.

   Clock angles are angles (right ascensions) measured anti-clockwise from the + X-axis. The clock angle for the point is:

   :math:`atan2(y_P,x_P)`.

   .. figure:: cone_clock_angle.png
      :scale: 100 %
      :align: center   

   The number of values in :code:`customConeAnglesVector` and :code:`customClockAnglesVector` should be the same (except for the case of Circular-shaped geometry in which case see note below).
   The last point of both the vectors should be the same as the first point to ensure polygon closure.

   .. csv-table:: Input parameter description 
      :header: Parameter, Type, Units, Description
      :widths: 10,10,10,40

      customConeAnglesVector, string, degrees, array of cone angle values separated by commas
      customClockAnglesVector, string, degrees, array of clock values separated by commas  
   
   .. note:: In case of circular-shaped spherical geometry, the number of values in :code:`customConeAnglesVector` should be just one (half the circular diameter) and 
             no values in :code:`customClockAnglesVector`.

Example:

.. code-block:: python

               "fieldOfViewGeometry": {
                                          "shape": "CUSTOM",
                                          "customConeAnglesVector": [10,10,10,10,10],
                                          "customClockAnglesVector": [30,120,180,280,30]
                                       }
               "fieldOfViewGeometry": {
                                       "shape": "CUSTOM", 
                                       "customConeAnglesVector": 15, 
                                       "@id": 123}

.. _sceneFieldOfViewGeometry_json_obj:

:code:`sceneFieldOfViewGeometry` json object
==============================================
The scene-field-of-view (sceneFOV) spherical geometry specification characterizes a (approximate) FOV representation of an image 'scene'. 
For example, in the case of stripmap SARs, or pushbroom optical scanners, a scene consists of multiple concatenated narrow strips (in the along-track direction). An 
approximate FOV representation can be specified to represent the observation.  If the sceneFOV geometry is not defined, 
the sceneFOV geometry is assigned to be equal to the instrument FOV geometry. 

The purpose of the sceneFOV is to enable faster coverage calculations. Always the sceneFOV or the FOR is considered for coverage calculations in the
OrbitPy package.

The json structure is identical to the :code:`fieldOfViewGeometry` JSON (see :ref:`fieldOfViewGeometry_json_obj`).

.. _maneuver_json_object:

:code:`maneuver` json object
========================================
This json object specified the total maneuverability of sensor pointing (combining satellite and sensor maneuverability) in the *NADIR_POINTING* reference frame. 
Three types of maneuvers are accepted: `CIRCULAR`, `SINGLE_ROLL_ONLY` and `DOUBLE_ROLL_ONLY`. This should be indicated in the 
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
Such an option is expected for instruments which require a pure-side-looking target geometry such as cross-track scanning radiometers.
The range of possible roll is indicated by the :code:`rollMin` and :code:`rollMax` name, value pairs. Note that these angles are
defined with respect to the *NADIR_POINTING* frame.

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
and the negative roll region) of the nadir-pointing frame (e.g.: synthetic aperture radars). 

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
List of orientations to which the instrument axis can be manuevered. Only the *NADIR_POINTING* reference frame is supported.
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

.. _antenna_json_object:

:code:`antenna` JSON object
==========================================
This json object contains the specifications of the antenna. Two types of antenna-aperture shapes are accepted, which should be indicated in the ``shape``
key/value pair.

1. :code:`"shape":"Circular"`

   This option indicates that the shape of the antenna-aperture is circular.

   .. csv-table:: Expected parameters
      :header: Parameter, Data type, Units, Description
      :widths: 10,10,5,40

      shape, str,, Must be "Circular"
      diameter, float, meters, Diameter of the antenna.
      apertureExcitationProfile, str, , Antenna aperture excitation profile. Accepted values are "UNIFORM" and "COSINE".
      apertureEfficiency, float,, Aperture efficiency (:math:`0 < \eta_{ap} < 1`).
      radiationEfficiency, float,, Radiation efficiency (:math:`0 < \psi < 1`).
      phyTemp, float, Kelvin, Physical temperature of the antenna.


   Example:

   .. code-block:: python
      
      "antenna":{
         "shape":"Circular",
         "diameter": 25,
         "apertureExcitationProfile": "COSINE",
         "apertureEfficiency": 0.6,
         "radiationEfficiency": 0.8,
         "phyTemp": 290
      }      

2. :code:`"shape":"Rectangular"`

   This option indicates that the shape of the antenna-aperture is rectangular.

   .. csv-table:: Expected parameters
      :header: Parameter, Data type, Units, Description
      :widths: 10,10,5,40

      shape, str,, Must be "Circular"
      height, float, meters, Antenna height (along the along-track direction when *SENSOR_BODY_FIXED* is aligned to *NADIR_POINTING* frame).
      width, float, meters, Antenna width (along the cross-track direction when *SENSOR_BODY_FIXED* is aligned to *NADIR_POINTING* frame).
      apertureExcitationProfile, str, , Antenna aperture excitation profile. Accepted values are "UNIFORM" and "COSINE".
      apertureEfficiency, float,, Aperture efficiency.
      radiationEfficiency, float,, Radiation efficiency.
      phyTemp, float, Kelvin, Physical temperature of the antenna.

   Example:

   .. code-block:: python
      
      "antenna":{
         "shape":"rectangular",
         "height": 4.9,
         "width": 0.7,
         "apertureExcitationProfile": "UNIFORM",
         "apertureEfficiency": 0.6,
         "radiationEfficiency": 0.8,
         "phyTemp": 290
      }

.. todo:: The operating frequency is not made as a specification of the antenna. Change behavior in the future?

.. _syntheticDataConfig_json_obj:

:code:`syntheticDataConfig` JSON object
================================================
This JSON object is used to describe the configuration of the synthetic data to be produced by the instrument models. A source data file containing gridded geophysical data
in netCDF format and the name of the geophysical variable appropriate to the instrument is required as input. 
*SCIPY_LINEAR* and *METPY_LINEAR* are inbuilt interpolation methods which can be invoked for interpolation. 

.. csv-table:: Expected parameters
   :header: Parameter, Data type, Units, Description
   :widths: 10,10,5,40

   sourceFilePaths, list str,, List of absolute filepaths of the science-data files in NetCDF format. Each file corresponds to a specific (forecast/analysis) time.
   geophysicalVar, str,, Geophysical variable (name as present in the source NetCDF file) to be used for the synthetic data.
   interpolMethod, str,, Interpolation method to be employed while interpolating the source data onto the pixel-positions. Allowed values are: *SCIPY_LINEAR* and *METPY_LINEAR*. 

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


.. _mode_json_obj:

:code:`mode` JSON object format
================================
The ``mode`` json object is used when initializing an instrument with several modes using the :class:`instrupy.base.Instrument` class. 
Several modes (in a list) maybe specified within a single instrument. Each mode corresponds to a specific operating point. For example, 
consider a *Basic Sensor* instrument which operates at two look-angles: (1) nadir-look (2) side-look at 30 deg. 
Such an instrument is considered to be made up of two modes with one mode specifying the nadir-look and the other mode specifying the side-look.
A mode-identifier can be specified by the user with which the corresponding mode can be referenced.

.. csv-table:: Input parameter description 
   :header: Parameter, Type, Units, Description
   :widths: 10,10,10,40

   @id, string,, Unique identifier of mode.

The parameters outside the mode block are used as the common parameters for all the modes, while the parameters specified
within a mode list entry are specific to the particular mode.

Example:

.. code-block:: python

               specs = '{  "@type": "Basic Sensor",
                           "name": "Atom",
                           "@id": "senX",  
                           "mass": 28, 
                           "volume": 0.12, 
                           "power": 32, 
                           "bitsPerPixel": 8, 
                           "fieldOfViewGeometry": {
                                       "shape": "CIRCULAR",
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
