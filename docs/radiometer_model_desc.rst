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
------------------
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
The radiometer system refers to the electronics configuration from the antenna to the output of the integrator. It includes the pre-detection stage, detector and the integrator. 
Following system-types can be modelled: *TOTAL_POWER*, *UNBALANCED_DICKE*, *BALANCED_DICKE* or *NOISE_ADDING*. The expected key/value pairs for each system type is described below:

1. :code:`"@type":"TOTAL_POWER"` 

   The expected key/value pairs for a *TOTAL_POWER* system are given below:

   .. csv-table:: Input parameter description 
    :header: Parameter, Data type, Units, Description
    :widths: 10,10,5,40

        @type, string, ,Must be *Basic Sensor*
        @id, string, , Unique identifier for the instrument. If ``None`` a random string is assigned.
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
        integratorVoltageGain, float, , Integrator voltage gain.
        predetectionGain, float, decibels, Pre-detection stage gain.
        predetectionInpNoiseTemp, float, Kelvin, Pre-detection *input noise* temperature.
        predetectionGainVariation, float, , Pre-detection stage gain variation. Linear units.
        integrationTime, float, seconds, Integration time.
        bandwidth, float, Hertz, Pre-detection bandwidth.



   Example:

   .. code-block:: python
      
      "swathConfig":{
            "@type": "full"
      }

.. _radiometer_scan_json_object:

:code:`scan` JSON object
----------------------------------
Two configurations (types) are accepted: *FULL* and *FIXED*.


Model results
------------------

Model description
------------------

Definition of the predetection stage:

From Pg 273, Fig.7-13 in [1] , the predetection stage includes all subsystems between the antenna and the input terminals of the square-law detector.
The specifications of the radiometric system can be made by either defining the specification of the entire predetection stage or of their individual components.

.. _radiometer_glossary:

Glossary
---------

Examples
---------
Please see the ``examples`` folder.