.. _radiometer_model_desc:

Radiometer Model
*****************
The radiometer model is based on the reference listed below. The following system types are supported: total-power, 
unbalanced-Dicke, balanced-Dicke and noise-adding. The following scan types are supported: fixed (no scan), cross-track 
and conical.

The FOV of the instrument is calculated from the antenna specifications (beamwidth), *scan-type* and the instrument orientation. 
A sceneFOV can be specified separately. The FOR is built based on the sceneFOV and the maneuver specifications. 
The FOV/ sceneFOV/ FOR is used in the coverage calculations (using the OrbitPy package) to find the locations accessed on the ground.

.. todo:: Field-of-view for conical-scan radiometers.

.. note:: See :ref:`radiometer_glossary` for names of the variables used in any discussion below.

**References:**

1. Chapter 6,7 in "Microwave Radar and Radiometric Remote Sensing," David Gardner Long , Fawwaz T. Ulaby 2014 

Model parameters
=================
A ``RadiometerModel`` object can be obtained from a json/ dict by using the ``from_json(.)`` or ``from_dict(.)`` functions. The expected key/value
pairs are described below:

.. csv-table:: Input parameter description 
    :header: Parameter, Data type, Units, Description
    :widths: 10,10,5,40

    @type, string, ,Must be *Basic Sensor*
    @id, string, , Unique identifier for the instrument. If ``None`` a random string is assigned.
    name, string, ,Full name of the instrument 
    mass, float, kilograms, Total mass of this entity.
    volume, float, :math:`m^3`, Total volume of this entity.
    power, float, Watts, Nominal operating power.
    orientation, :ref:`orientation_json_obj`, ,Orientation of the instrument. Default is alignment to the SC_BODY_FIXED frame.
    fieldOfViewGeometry, :ref:`fieldOfViewGeometry_json_obj`, , Field of view spherical geometry specification of the instrument. 
    sceneFieldOfViewGeometry, :ref:`sceneFieldOfViewGeometry_json_obj`, , The SceneFOV spherical geometry specification of the instrument. Default is the field-of-view spherical geometry specification.
    maneuver, :ref:`maneuver_json_object`, , Maneuver specifications (see :ref:`maneuv_desc`).
    pointingOption, :ref:`pointing_opt_json_obj`, , List of orientations to which the instrument axis can be maneuvered.
    dataRate, float, Mega-bits-per-s, Rate of data recorded during nominal operations.
    bitsPerPixel, integer, ,Bits encoded per pixel of image.
    antenna, :ref:`antenna_json_object`, , Antenna specifications. Only rectangular shape and uniform aperture excitation profile is accepted.
    operatingFrequency, float, Hertz, Operating radar center frequency.
    system, :ref:`radiometer_sys_json_object`, , Radiometer system.
    scan, :ref:`radiometer_scan_json_object`, , Scan specifications. Default is a *Fixed* specification (no-scan).
    targetBrightnessTemperature, float, Kelvin, Target brightness temperature. Default value is 290K.

.. _radiometer_sys_json_object:

:code:`system` JSON object
----------------------------------
The radiometer system refers to the electronics configuration from the antenna to the output of the integrator.
Following system-types can be modelled: *TOTAL_POWER*, *UNBALANCED_DICKE*, *BALANCED_DICKE* or *NOISE_ADDING*. 

The ``antenna``, ``operatingFrequency`` and the ``targetBrightnessTemperature`` specifications of the system are obtained as external inputs 
when required to compute some system parameters and the radiometric performance. 
In each of the systems, the predetection stage parameters can be specified in two ways: (1) component-level specification or (2) black-box specification.

The key/value pairs of each of the system types is described below:

1. :code:`"@type":"TOTAL_POWER"` 

    The expected key/value pairs for a total-power radiometer system system, **excluding** that of the predetection-stage are given below.

    .. csv-table:: Common parameters
        :header: Parameter, Data type, Units, Description
        :widths: 10,10,5,40

        @type, string, ,Must be *TOTAL_POWER*
        integrationTime, float, seconds, Integration time.
        bandwidth, float, Hertz, Pre-detection bandwidth.
        integratorVoltageGain, float, , Integrator voltage gain.
    
    Below are the expected key/value pairs of the predetection stage (black-box specification).

    .. csv-table:: Predetection stage parameters (black-box specification)
        :header: Parameter, Data type, Units, Description
        :widths: 10,10,5,40

        predetectionGain, float, decibels, Pre-detection stage gain.
        predetectionInpNoiseTemp, float, Kelvin, Pre-detection *input noise* temperature.
        predetectionGainVariation, float, , Pre-detection stage gain variation. Linear units.
    
    Below are the expected key/value pairs of the predetection stage (component-level specification).

    .. csv-table:: Predetection stage parameters (component-level specification)
        :header: Parameter, Data type, Units, Description
        :widths: 10,10,5,40

        tlLoss, float, decibels, Transmission line loss.
        tlPhyTemp, float, Kelvin, Transmission line *physical* temperature.
        rfAmpGain, float, decibels, RF amplifier gain.
        rfAmpInpNoiseTemp, float, Kelvin, RF amplifier *input noise* temperature.
        rfAmpGainVariation, float, , RF amplifier gain variation. Linear units.
        mixerGain, float, decibels, Mixer gain.
        mixerInpNoiseAmp, float, Kelvin, Mixer *input noise* temperature.
        mixerGainVariation, float, , Mixer gain variation. Linear units.
        ifAmpGain, float, decibels, Intermediate frequency amplifier gain.
        ifAmpInpNoiseTemp, float, Kelvin, Intermediate frequency amplifier *input noise* temperature.
        ifAmpGainVariation, float, , IF amplifier gain variation. Linear units.       

    Example:

    *Total-power System with component-level predetection-stage specification.*

    .. code-block:: python

        tpr_sys1_json = {"tlLoss": 0.5,
                         "tlPhyTemp": 290,
                         "rfAmpGain": 30,
                         "rfAmpInpNoiseTemp": 200,
                         "rfAmpGainVariation": 10,
                         "mixerGain": 23,
                         "mixerInpNoiseTemp": 1200,
                         "mixerGainVariation": 2,
                         "ifAmpGain": 30,
                         "ifAmpInputNoiseTemp": 100,
                         "ifAmpGainVariation": 10,
                         "integratorVoltageGain": 1,
                         "integrationTime": 100e-3,
                         "bandwidth": 10e6,
                        }

    
    *Total-power System with block-box predetection-stage specification.*

    .. code-block:: python

        tpr_sys2_json = {"predetectionGain": 83,
                         "predetectionInpNoiseTemp": 200,
                         "predetectionGainVariation": 2000000,
                         "integrationTime": 100e-3,
                         "bandwidth": 10e6,
                         "integratorVoltageGain": 1 
                        }

2. :code:`"@type":"UNBALANCED_DICKE"` 

    The expected key/value pairs for a unbalanced-Dicke radiometer system system is similar to the *TOTAL_POWER* system. 
    The ``referenceTemperature`` is an additional key/value pairs required as compared to the total-power radiometer system.

    .. csv-table:: Common parameters
        :header: Parameter, Data type, Units, Description
        :widths: 10,10,5,40

        @type, string, ,Must be *UNBALANCED_DICKE*
        integrationTime, float, seconds, Integration time.
        bandwidth, float, Hertz, Pre-detection bandwidth.
        integratorVoltageGain, float, , Integrator voltage gain.
        referenceTemperature, float, Kelvin, Reference source noise temperature.
    
    Below are the expected key/value pairs of the predetection stage (black-box specification). 

    .. csv-table:: Predetection stage parameters (black-box specification)
        :header: Parameter, Data type, Units, Description
        :widths: 10,10,5,40

        predetectionGain, float, decibels, Pre-detection stage gain.
        predetectionInpNoiseTemp, float, Kelvin, Pre-detection *input noise* temperature.
        predetectionGainVariation, float, , Pre-detection stage gain variation. Linear units.
    
    Below are the expected key/value pairs of the predetection stage (component-level specification).
    The ``dickeSwitchOutputNoiseTemperature`` is an additional key/value pairs required as compared to the total-power radiometer system.

    .. csv-table:: Predetection stage parameters (component-level specification)
        :header: Parameter, Data type, Units, Description
        :widths: 10,10,5,40

        tlLoss, float, decibels, Transmission line loss.
        tlPhyTemp, float, Kelvin, Transmission line *physical* temperature.
        rfAmpGain, float, decibels, RF amplifier gain.
        rfAmpInpNoiseTemp, float, Kelvin, RF amplifier *input noise* temperature.
        rfAmpGainVariation, float, , RF amplifier gain variation. Linear units.
        mixerGain, float, decibels, Mixer gain.
        mixerInpNoiseAmp, float, Kelvin, Mixer *input noise* temperature.
        mixerGainVariation, float, , Mixer gain variation. Linear units.
        ifAmpGain, float, decibels, Intermediate frequency amplifier gain.
        ifAmpInpNoiseTemp, float, Kelvin, Intermediate frequency amplifier *input noise* temperature.
        ifAmpGainVariation, float, , IF amplifier gain variation. Linear units.
        dickeSwitchOutputNoiseTemperature, float, Kelvin, Dicke switch noise temperature *referenced to the output port.*

3. :code:`"@type":"BALANCED_DICKE"` 

    The expected key/value pairs for a balanced-Dicke radiometer system system is similar to the *TOTAL_POWER* system. 

    .. csv-table:: Common parameters
        :header: Parameter, Data type, Units, Description
        :widths: 10,10,5,40

        @type, string, ,Must be *BALANCED_DICKE*
        integrationTime, float, seconds, Integration time.
        bandwidth, float, Hertz, Pre-detection bandwidth.
        integratorVoltageGain, float, , Integrator voltage gain.
    
    Below are the expected key/value pairs of the predetection stage (black-box specification). 

    .. csv-table:: Predetection stage parameters (black-box specification)
        :header: Parameter, Data type, Units, Description
        :widths: 10,10,5,40

        predetectionGain, float, decibels, Pre-detection stage gain.
        predetectionInpNoiseTemp, float, Kelvin, Pre-detection *input noise* temperature.
        predetectionGainVariation, float, , Pre-detection stage gain variation. Linear units.
    
    Below are the expected key/value pairs of the predetection stage (component-level specification).
    The ``dickeSwitchOutputNoiseTemperature`` is an additional key/value pairs required as compared to the total-power radiometer system.

    .. csv-table:: Predetection stage parameters (component-level specification)
        :header: Parameter, Data type, Units, Description
        :widths: 10,10,5,40

        tlLoss, float, decibels, Transmission line loss.
        tlPhyTemp, float, Kelvin, Transmission line *physical* temperature.
        rfAmpGain, float, decibels, RF amplifier gain.
        rfAmpInpNoiseTemp, float, Kelvin, RF amplifier *input noise* temperature.
        rfAmpGainVariation, float, , RF amplifier gain variation. Linear units.
        mixerGain, float, decibels, Mixer gain.
        mixerInpNoiseAmp, float, Kelvin, Mixer *input noise* temperature.
        mixerGainVariation, float, , Mixer gain variation. Linear units.
        ifAmpGain, float, decibels, Intermediate frequency amplifier gain.
        ifAmpInpNoiseTemp, float, Kelvin, Intermediate frequency amplifier *input noise* temperature.
        ifAmpGainVariation, float, , IF amplifier gain variation. Linear units.
        dickeSwitchOutputNoiseTemperature, float, Kelvin, Dicke switch noise temperature *referenced to the output port.*

4. :code:`"@type":"NOISE_ADDING"` 

    The expected key/value pairs for a noise-adding radiometer system system is similar to the *TOTAL_POWER* system. 
   
    The ``excessNoiseTemperature`` is an additional key/value pair required as compared to the total-power radiometer system.
    
    .. csv-table:: Common parameters
        :header: Parameter, Data type, Units, Description
        :widths: 10,10,5,40

        @type, string, ,Must be *NOISE_ADDING*
        integrationTime, float, seconds, Integration time.
        bandwidth, float, Hertz, Pre-detection bandwidth.
        integratorVoltageGain, float, , Integrator voltage gain.
        excessNoiseTemperature, float, Kelvin, Excess noise temperature (added noise to the receiver input during the diode ON half-cycle) in Kelvin *referenced to the output port.*
    
    Below are the expected key/value pairs of the predetection stage (black-box specification). 

    .. csv-table:: Predetection stage parameters (black-box specification)
        :header: Parameter, Data type, Units, Description
        :widths: 10,10,5,40

        predetectionGain, float, decibels, Pre-detection stage gain.
        predetectionInpNoiseTemp, float, Kelvin, Pre-detection *input noise* temperature.
        predetectionGainVariation, float, , Pre-detection stage gain variation. Linear units.
    
    Below are the expected key/value pairs of the predetection stage (component-level specification).

    .. csv-table:: Predetection stage parameters (component-level specification)
        :header: Parameter, Data type, Units, Description
        :widths: 10,10,5,40

        tlLoss, float, decibels, Transmission line loss.
        tlPhyTemp, float, Kelvin, Transmission line *physical* temperature.
        rfAmpGain, float, decibels, RF amplifier gain.
        rfAmpInpNoiseTemp, float, Kelvin, RF amplifier *input noise* temperature.
        rfAmpGainVariation, float, , RF amplifier gain variation. Linear units.
        mixerGain, float, decibels, Mixer gain.
        mixerInpNoiseAmp, float, Kelvin, Mixer *input noise* temperature.
        mixerGainVariation, float, , Mixer gain variation. Linear units.
        ifAmpGain, float, decibels, Intermediate frequency amplifier gain.
        ifAmpInpNoiseTemp, float, Kelvin, Intermediate frequency amplifier *input noise* temperature.
        ifAmpGainVariation, float, , IF amplifier gain variation. Linear units.

.. _radiometer_scan_json_object:

:code:`scan` JSON object
----------------------------------
Three scan-techniques are supported: *FIXED* (no-scan), *CROSS_TRACK* and *CONICAL*. The scan-technique determines the instrument field-of-view (and hence the swath-width), 
dwell-time (and hence the maximum integration-time).

1. :code:`"@type":"FIXED"`
   
   This scan-technique specifies that there is no scan. The antenna (or the feeder) is held fixed with respect to the spacecraft. No parameters are required.

   Example:

   .. code-block:: python
      
      "scan":{
            "@type": "FIXED"
      }
 
2. :code:`"@type":"CROSS_TRACK"`
   
    In this scan-technique the antenna foot-print is scanned in the cross-track direction. The ``scanWidth`` parameter specifies the angular width
    about the instrument orientation (which in general is *SIDE_LOOK*), while the ``interScanOverheadTime`` specifies the time taken to go from scan of 
    one strip (in the cross-track direction) to the next. 

    .. csv-table:: 
        :header: Parameter, Data type, Units, Description
        :widths: 10,10,5,40

        scanWidth, float, degrees, Angular scan-width.
        interScanOverheadTime, float, seconds, Time taken from ending current scan to starting next scan. Significant in case of mechanical scanning. Default value is 0.
    
    Example:

    .. code-block:: python

       "scan":{
                "@type": "CROSS_TRACK, 
                "scanWidth": 120, 
                "interScanOverheadTime": 1e-3
              }

3. :code:`"@type":"CONICAL"`

    In this scan-technique the antenna footprint is scanned along the cone-perimeter. The ``offNadirAngle``specifies the (half) cone angle while 
    the ``clockAngleRange`` parameter specifies the azimuth extent of the scan (symmetrically about the along-track direction). 
    The ``interScanOverheadTime`` specifies the time taken to go from scan of one strip to the next. 

    For illustration of off-nadir angle and clock angles see Fig.7 in T. Kawanishi et al., "The Advanced Microwave Scanning Radiometer for the Earth Observing System (AMSR-E), NASDA's contribution to the EOS for global energy and water cycle studies," in IEEE Transactions on Geoscience and Remote Sensing, vol. 41, no. 2, pp. 184-194, Feb. 2003, doi: 10.1109/TGRS.2002.808331.

    .. csv-table:: 
        :header: Parameter, Data type, Units, Description
        :widths: 10,10,5,40

        offNadirAngle, float, degrees, Off-nadir angle (i.e. the half-cone angle of the conical scan).
        clockAngleRange, float, degrees, Scan clock angle range in degrees.
        interScanOverheadTime, float, seconds, Time taken from ending current scan to starting next scan. Significant in case of mechanical scanning. Default value is 0.
    
    Example:

    .. code-block:: python

        "scan": {
                 "offNadirAngle": 30, 
                 "clockAngleRange": 60, 
                 "interScanOverheadTime": 1e-3
                 }


Model results
------------------
Using the radiometer model, coverage calculations (using the OrbitPy package) can be carried out over a region of interest. Coverage calculations which involve 
a grid (list of grid-points) evaluate to see if the grid-points fall within the instrument sceneFOV (sceneFOV = FOV in most cases) or the FOR. The pointing-options feature further 
allows to automate coverage calculations for numerous instrument orientations. 

Once the coverage has been evaluated, the observable locations and the observer (satellite) locations is known. The following data metrics at the observable location 
on the surface of Earth can be calculated:

.. csv-table:: Observation data metrics table
    :widths: 8,4,4,20
    :header: Metric/Aux data,Data Type,Units,Description

    radiometric res [K], float, Kelvin, Radiometric resolution/ sensitivity.
    ground pixel along-track resolution [m], float, meters, Along-track resolution of an ground-pixel centered about observation point.
    ground pixel cross-track resolution [m], float, meters, Cross-track resolution of an ground-pixel centered about observation point.
    swath-width [m], float, meters, Swath-width of the strip of which the imaged pixel is part off.
    beam efficiency, float, ,Beam efficiency of the antenna.
    incidence angle [deg], float, degrees, Observation incidence angle at the ground-pixel.

.. note:: Coverage calculations for radiometers with conical-scan is currently not supported unless a sceneFOV has been explicitly specified.

.. todo:: The along-track and cross-track pixel resolutions are accurate only pixels imaged at strictly sidelooking geometry (roll-only, no pitch). Needs revision.

Model description
------------------

Definition of the predetection stage:

From Pg 273, Fig.7-13 in [1] , the predetection stage includes all subsystems between the antenna and the input terminals of the square-law detector.
The specifications of the radiometric system can be made by either defining the specification of the entire predetection stage or of their individual components.

. note:: The swath-width is calculated more precisely as compared to the pixel-resolution calculations. This leads to a slight difference
             as can be seen in the below results. 



.. _radiometer_glossary:

Glossary
---------

Examples
---------
Please see the ``examples`` folder.