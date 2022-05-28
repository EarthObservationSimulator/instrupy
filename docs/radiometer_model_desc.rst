.. _radiometer_model_desc:

Radiometer Model
*****************
The radiometer model is based on the reference listed below. The following system types are supported: total-power, 
unbalanced-Dicke, balanced-Dicke and noise-adding. The following scan types are supported: fixed (no scan), cross-track 
and conical.

The FOV of the instrument is calculated from the antenna specifications (beamwidth), *scan-type* and the instrument orientation. 
A sceneFOV can be specified separately. The FOR is built based on the sceneFOV and the maneuver specifications. 
The sceneFOV/ FOR is used in the coverage calculations (using the OrbitPy package) to find the locations accessed on the ground.

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

    @type, string, ,Must be *Radiometer*
    @id, string, , Unique identifier for the instrument. If ``None`` a random string is assigned.
    name, string, ,Full name of the instrument 
    mass, float, kilograms, Total mass of this entity.
    volume, float, :math:`m^3`, Total volume of this entity.
    power, float, Watts, Nominal operating power.
    orientation, :ref:`orientation_json_obj`, ,Orientation of the instrument. Default is alignment to the SC_BODY_FIXED frame.
    fieldOfViewGeometry, :ref:`fieldOfViewGeometry_json_obj`, , Field of view spherical geometry specification of the instrument. 
    sceneFieldOfViewGeometry, :ref:`sceneFieldOfViewGeometry_json_obj`, , The SceneFOV spherical geometry specification of the instrument. Default is the field-of-view spherical geometry.
    maneuver, :ref:`maneuver_json_object`, , Maneuver specifications (see :ref:`maneuv_desc`).
    pointingOption, :ref:`pointing_opt_json_obj`, , List of orientations to which the instrument axis can be maneuvered.
    dataRate, float, Mega-bits-per-s, Rate of data recorded during nominal operations.
    bitsPerPixel, integer, ,Bits encoded per pixel of image.
    antenna, :ref:`antenna_json_object`, , Antenna specifications.
    operatingFrequency, float, Hertz, Operating center frequency.
    system, :ref:`radiometer_sys_json_object`, , Radiometer system.
    scan, :ref:`radiometer_scan_json_object`, , Scan specifications. Default is a *FIXED* specification (no-scan).
    targetBrightnessTemperature, float, Kelvin, Target brightness temperature. Default value is 290K.

.. _radiometer_sys_json_object:

:code:`system` JSON object
----------------------------------
The radiometer-system refers to the electronics configuration from the antenna to the output of the integrator.
Following system-types can be modelled: *TOTAL_POWER*, *UNBALANCED_DICKE*, *BALANCED_DICKE* or *NOISE_ADDING*. 

The ``antenna``, ``operatingFrequency`` and the ``targetBrightnessTemperature`` specifications of the system are obtained as external inputs 
when required to compute some system parameters and the radiometric performance. 
In each of the systems, the predetection stage parameters can be specified in two ways: (1) component-level specification or (2) black-box specification.

The key/value pairs of each of the system types is described below:

1. :code:`"@type":"TOTAL_POWER"` 

    The expected key/value pairs for a total-power radiometer system system **excluding** that of the predetection-stage are given below.

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

    The expected key/value pairs for a unbalanced-Dicke radiometer system system **excluding** that of the predetection-stage consists of all the kep/value pairs
    of the *TOTAL_POWER* system **and** the ``referenceTemperature`` key/value pair. The ``@type`` key must have "UNBALANCED_DICKE" as the value.

    .. csv-table:: Common parameters
        :header: Parameter, Data type, Units, Description
        :widths: 10,10,5,40

        @type, string, ,Must be *UNBALANCED_DICKE*
        referenceTemperature, float, Kelvin, Reference source noise temperature.
    
    The expected key/value pairs of the predetection stage (black-box specification) is the same as that of the *TOTAL_POWER* system.
    
    The expected key/value pairs of the predetection stage (component-level specification) consists of all the key/value pairs
    of the *TOTAL_POWER* system **and** the ``dickeSwitchOutputNoiseTemperature`` key/value pair.

    .. csv-table:: Predetection stage parameters (component-level specification)
        :header: Parameter, Data type, Units, Description
        :widths: 10,10,5,40

        dickeSwitchOutputNoiseTemperature, float, Kelvin, Dicke switch noise temperature *referenced to the output port.*

3. :code:`"@type":"BALANCED_DICKE"` 

    The expected key/value pairs for a balanced-Dicke radiometer system system is similar to the *TOTAL_POWER* system. 

    The expected key/value pairs for a unbalanced-Dicke radiometer system system **excluding** that of the predetection-stage consists of all the kep/value pairs
    of the *TOTAL_POWER* system. The ``@type`` key must have "BALANCED_DICKE" as the value.

    .. csv-table:: Common parameters
        :header: Parameter, Data type, Units, Description
        :widths: 10,10,5,40

        @type, string, ,Must be *BALANCED_DICKE*
    
    The expected key/value pairs of the predetection stage (black-box specification) is the same as that of the *TOTAL_POWER* system. 
    
    The expected key/value pairs of the predetection stage (component-level specification) consists of all the key/value pairs
    of the *TOTAL_POWER* system **and** the ``dickeSwitchOutputNoiseTemperature`` key/value pair.

    .. csv-table:: Predetection stage parameters (component-level specification)
        :header: Parameter, Data type, Units, Description
        :widths: 10,10,5,40

        dickeSwitchOutputNoiseTemperature, float, Kelvin, Dicke switch noise temperature *referenced to the output port.*

4. :code:`"@type":"NOISE_ADDING"` 

    The expected key/value pairs for a noise-adding radiometer system system is similar to the *TOTAL_POWER* system. 
   
    The expected key/value pairs for a unbalanced-Dicke radiometer system system **excluding** that of the predetection-stage consists of all the kep/value pairs
    of the *TOTAL_POWER* system **and** the ``excessNoiseTemperature`` key/value pair. The ``@type`` key must have "NOISE_ADDING" as the value.
    
    .. csv-table:: Common parameters
        :header: Parameter, Data type, Units, Description
        :widths: 10,10,5,40

        @type, string, ,Must be *NOISE_ADDING*
        excessNoiseTemperature, float, Kelvin, Excess noise temperature (added noise to the receiver input during the diode ON half-cycle) in Kelvin *referenced to the output port.*
    
    The expected key/value pairs of the predetection stage (black-box specification) is the same as that of the *TOTAL_POWER* system. 
    
    The expected key/value pairs of the predetection stage (component-level specification) is the same as that of the *TOTAL_POWER* system. 

.. _radiometer_scan_json_object:

:code:`scan` JSON object
----------------------------------
Three scan-techniques are supported: *FIXED* (no-scan), *CROSS_TRACK* and *CONICAL*. The scan-technique determines the instrument field-of-view (and hence the swath-width), 
dwell-time (and hence the maximum integration-time).

1. :code:`"@type":"FIXED"`
   
   This scan-technique specifies that there is no scan. The antenna (and the feeder) is held fixed with respect to the spacecraft. No parameters are required.

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

    In this scan-technique the antenna footprint is scanned along the cone-perimeter. The ``offNadirAngle`` specifies the (half) cone angle while 
    the ``clockAngleRange`` parameter specifies the azimuth extent of the scan (symmetrically about the along-track direction). 
    The ``interScanOverheadTime`` specifies the time taken to go from scan of one strip to the next. 

    For illustration of off-nadir angle and clock angles see Fig.7 in T. Kawanishi et al., "The Advanced Microwave Scanning Radiometer for the Earth Observing System (AMSR-E), NASDA's contribution to the EOS for global energy and water cycle studies," in IEEE Transactions on Geoscience and Remote Sensing, vol. 41, no. 2, pp. 184-194, Feb. 2003.

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
===============
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

.. todo:: The along-track and cross-track pixel resolutions are accurate only for pixels imaged at strictly sidelooking geometry (roll-only, no pitch). Needs revision.

Model description
===================
Below text lays down the formulae coded into the model based on reference [1].

Viewing geometry
-----------------------
The viewing geometry parameters, i.e. :math:`\mathbf{S}`, :math:`\mathbf{T}`, :math:`\mathbf{R}`, :math:`\theta_i` and :math:`\gamma` are determined using the setup 
described in :ref:`basic sensor model description<basic_sensor_model_desc>`.

Pixel-resolutions
-------------------
Note that the current formulation is accurate only when ground-pixel is being imaged at the nadir or is at purely side-looking geometry.

:math:`\rho_{at} = R \mu_{at}`

:math:`\rho_{ct} = R \mu_{ct}/ \cos(\theta_i)`

.. todo:: Update for the general target geometry. 

Radiometric resolution
-----------------------

Integration time calculation
...................................
The dwell-time :math:`t_d` of the antenna over a pixel gives the maximum possible integration time. It depends on the scan technique:

1. *FIXED* 
   
   :math:`t_d = \rho_{at}/v_{g}`

2. *CROSS_TRACK* 
   
    :math:`n_{pps} = \Delta_{asw} / \rho_{ct}`
    
    :math:`t_d = \dfrac{\rho_{at}/v_{g} - \Delta_{is}}{n_{pps}}`

3. *CONICAL* 

    :math:`n_{pps} = \Delta_{car} / \rho_{ct}`

    :math:`t_d = \dfrac{\rho_{at}/v_{g} - \Delta_{is}}{n_{pps}}`

**FInally,** if the calculated dwell time is lesser than the user-defined integration-time, the integration-time is set to the calculated dwell time,
else the integration time is set to the user-specified integration-time.

:math:`if \hspace{2mm} \tau_{spec} > t_d, \hspace{2mm} \tau =  t_d` else :math:`\tau = \tau_{spec}`

Predetection section parameters
.................................
The predetection stage includes all subsystems between the antenna and the input terminals of the square-law detector (Pg 273, Fig.7-13 in [1]).
The specifications of the radiometric system can be made by either defining the specification of the entire predetection stage (as a black-box)
or of their individual components. 

*If the black-box specifications are provided:*

:math:`G_{PD}^- = G_{PD} - 0.5 \Delta G_{PD}`

:math:`G_{PD}^+ = G_{PD} + 0.5 \Delta G_{PD}`

*If the component-level specifications are provided:*

(Fig.7-9 in [1] describes the gain of the transmission line as 1/L, where L is the transmission line loss.)

:math:`G_{TL} = 1/L`  (transmission line "gain")

:math:`G_{PD} = G_{TL} G_{RF} G_{MIX} G_{IF}`
    
:math:`G_{PD}^- = G_{TL} * (G_{RF} - 0.5 \Delta G_{RF}) (G_{MIX} - 0.5 \Delta G_{MIX}) (G_{IF} - 0.5 \Delta G_{IF})` 

:math:`G_{PD}^+ = G_{TL} * (G_{RF} + 0.5 \Delta G_{RF}) (G_{MIX} + 0.5 \Delta G_{MIX})  (G_{IF} + 0.5 \Delta G_{IF})`

(See Section 7-3.1 in [1] for example calculation of noise temperature from cascaded stages.)

:math:`T_{REC} = T_{RF} + T_{MIX}/ G_{RF} + T_{IF}/ (G_{RF} G_{IF})`  (Eqn 7.29 in [1])

:math:`T'_{REC} = (L-1) T_{TL}^{P} + L T_{REC}`

In case of *UNBALANCED_DICKE* and *BALANCED_DICKE* radiometer-system:

:math:`T'_{REC} = T'_{REC} + T_{DSW}^o`

System parameters
....................

Calculate system gain factor (eqn 7.43 in [1]):

:math:`G_s = 2 G_{INT} G_{PD} k_B B`        

Calculate the system gain variation:

:math:`G_s^- = 2 G_{INT} G_{PD}^- k_B B`

:math:`G_s^+ = 2 G_{INT} G_{PD}^+ k_B B`

:math:`\Delta G_s = G_s^+ - G_s^-`

:math:`\bar{G_s} = G_s`  (average system power gain, TODO: check)

Calculate system temperature:

(antenna radiation efficiency (:math:`\psi`) = 1/ antenna loss)

:math:`T_A = \psi T'_A + (1-psi) T_A^p`

:math:`T_{SYS} = T_A + T'_{REC}`  (eqn 7.31 in [1])

Resolution calculation
........................
*TOTAL_POWER* radiometer system:

:math:`\Delta T = T_{SYS} \sqrt{\dfrac{1}{B \tau} + (\dfrac{\Delta G_{SYS}}{\bar{G_{SYS}}})^2}`

*UNBALANCED_DICKE* radiometer system:

:math:`\Delta T = \sqrt{\dfrac{2 T_{SYS}^2 + 2 (T_{REF} + T'_{REC})^2}{B \tau} + (\dfrac{\Delta G_{SYS}}{\bar{G_{SYS}}})^2 (T_A - T_{REF})^2}`

*BALANCED_DICKE* radiometer system:

:math:`\Delta T = 2 \dfrac{T_{SYS}}{\sqrt{B \tau}}`

*NOISE_ADDING* radiometer system:

:math:`\Delta T = 2 \dfrac{T_{SYS}}{\sqrt{B \tau}}  (1 + \dfrac{2 T_{SYS}}{T_{N}''})`

Instrument field-of-View spherical-geometry calculations
-----------------------------------------------------------
The instrument field-of-view depends on the chosen scan technique and antenna specifications.

*FIXED* scan:

The FOV spherical-geometry shape is determined by the antenna shape (*CIRCULAR* or *RECTANGULAR*).

:math:`\theta_{AT} = \mu_{at}`

:math:`\theta_{CT} = \mu_{ct}`

Note that for circular antenna shape :math`\mu_{at} = \mu_{ct}`.

*CROSS_TRACK* scan:

The FOV spherical-geometry shape is always *RECTANGULAR*.

:math:`\theta_{AT} = \mu_{at}`

:math:`\theta_{CT} = \mu_{ct} + \Delta_{asw}` 

*CONICAL* scan: 

TBD. The instrument orientation has to be nadir-pointing.


Swath-width
-------------

THe swath-width is calculated from the instrument look-angle and not the look-angle to the target ground-point.
The swath-width depends on the scan technique. 

*FIXED* and *CROSS_TRACK* scan:

In case of fixed-scan mode, there is only 1 imaged ground-pixel per swath. Swath-width is computed to be equal to the antenna-footprint cross-track size. 
See Fig.5.1.3.1 in Spaceborne SAR Study: LDRD 92 Final Report SANDIA Report March 1993.

:math:`R_S = R_E + h`   

:math:`\gamma_n = \gamma_I - 0.5 \hspace{1mm} \theta_{CT}`

:math:`\gamma_f = \gamma_I  + 0.5 \hspace{1mm} \theta_{CT}`

:math:`\theta_{in} = \sin^{-1}(\sin(\gamma_n) R_S/R_E)`

:math:`\theta_{if} = \sin^{-1}(\sin(\gamma_f) R_S/R_E)`

:math:`\alpha_n = \theta_{in} - \gamma_n`

:math:`\alpha_f = \theta_{if} - \gamma_f`

if :math:`\gamma_n` <= 0, the radiometer footprint falls in the nadir-direction, and we have:

:math:`\alpha_s = \alpha_f + \alpha_n`

if :math:`\gamma_n` > 0 we have:

:math:`\alpha_s = |\alpha_f - \alpha_n|`

:math:`W_{gr} = R_E \alpha_s`   

(:math:`\theta_{CT} = \mu_{ct}` for the vase of *FIXED* scan.)

.. note:: The swath-width is calculated more precisely as compared to the pixel-resolution calculations. This leads to a small but noticeable
    difference while examining the results of a *FIXED* scan radiometer, in which the swath-width should be equal to the pixel-size in the cross-track direction.

*CONICAL* scan:

Calculate the radius of the small-circle on the Earth surface on which the imaged arc lies.

:math:`\theta_i^{cs} = \sin^{-1}(\sin(\gamma^{cs}) \dfrac{R_S}{R_E})`

:math:`\alpha^{cs} = \theta_i^{cs} - \gamma^{cs}`

:math:`r^{cs} = R_E \sin{\alpha^{cs}}`

:math:`A^{cs} = \Delta_{car} r^{cs}` 

Beam-efficiency
----------------

Please refer to the antenna description.


Examples
=========
Please see the ``examples`` folder.

.. _radiometer_glossary:

Glossary
=========
* :math:`\mathbf{R}`: Range vector from satellite to target ground point.
* :math:`\theta_i`: Incidence angle at the target ground point.
* :math:`R_E`: Nominal equatorial radius of Earth.
* :math:`c`: Speed of light.
* :math:`h`: Altitude of the satellite.
* :math:`\lambda`: Operating center wavelength of the radiometer.
* :math:`\rho_{at}`: Along-track pixel resolution.
* :math:`\rho_{ct}`: Cross-track pixel resolution.
* :math:`\mu_{AT}`: Along-track *antenna* FOV.
* :math:`\mu_{CT}`: Cross-track *antenna* FOV.
* :math:`\theta_{AT}`: Along-track *instrument* FOV.
* :math:`\theta_{CT}`: Cross-track *instrument* FOV.
* :math:`v_g`: Ground speed of satellite footprint.
* :math:`t_d`: Dwell time available over the ground-pixel.
* :math:`\tau`: Integration time.
* :math:`\tau_{spec}`: Integration time specification from user.
* :math:`\Delta_{asw}`: Angular scan width in case of *CROSS_TRACK* scan.
* :math:`\Delta_{car}`: Clock angle range (of scan) in case of *CONICAL* scan.
* :math:`n_{pps}`: Number of pixels per strip in case of *CROSS_TRACK* and *CONICAL* scans.
* :math:`\Delta_{is}``: Overhead time to go switch scan from one strip to another in case of *CROSS_TRACK* and *CONICAL* scans.
* :math:`G_{PD}`: Predetection gain (linear units).
* :math:`G_{PD}^+`: Predetection gain + (linear units).
* :math:`G_{PD}^-`: Predetection gain - (linear units).
* :math:`\Delta G_{PD}`: Predetection Gain variation (linear units).
* :math:`L`: Transmission line loss (linear units).
* :math:`G_{TL}`: Transmission line gain.
* :math:`G_{RF}`: RF amplifier gain.
* :math:`G_{MIX}`: Mixer gain.
* :math:`G_{IF}`: IF (Intermediate frequency) amplifier gain.
* :math:`\Delta G_{RF}`: RF amplifier gain variation.
* :math:`\Delta G_{MIX}`: Mixer gain variation.
* :math:`\Delta G_{IF}`: IF amplifier gain variation.
* :math:`T_{REC}`: Predetection stage (*excluding* the transmission line from antenna to the RF amplifier) input noise temperature. (Receiver noise temperature.)
* :math:`T'_{REC}`: Predetection stage (*including* the transmission line from antenna to the RF amplifier) input noise temperature. (Receiver noise temperature referred to the antenna terminals.)
* :math:`T_{RF}`: RF amplifier input noise temperature.
* :math:`T_{MIX}`: Mixer input noise temperature.
* :math:`T_{IF}`: IF amplifier input noise temperature.
* :math:`T_{TL}^{P}`: Transmission line physical temperature.
* :math:`T_{DSW}^o` : Dicke switch output noise temperature.
* :math:`G_s`: System Gain (linear units).
* :math:`G_{INT}`: Integrator *voltage* gain.
* :math:`B`: Predetection bandwidth.
* :math:`k_B`: Boltzmann constant.
* :math:`G_s^-`: System Gain - (linear units).
* :math:`G_s^+`: System Gain + (linear units).
* :math:`\Delta G_s`: System Gain variation (linear units).
* :math:`\psi`: Antenna radiation efficiency (= 1/ antenna loss).
* :math:`T_A^p`: Antenna physical temperature.
* :math:`T'_A`: Scene brightness temperature :math:`T_B(\theta,\phi)`, weighted with the antenna pattern.
* :math:`T_A`: Antenna (radiometric) temperature referred at the output terminal of the antenna.
* :math:`\bar{G_s}`: Average system gain.
* :math:`T_{SYS}`: System noise temperature.
* :math:`\Delta T`: Radiometric resolution of the radiometer.
* :math:`T_{REF}``: Reference noise temperature for Dicke radiometer systems.
* :math:`T_{N}''`: Excess noise temperature for *NOISE_ADDING* radiometer system.
* :math:`R_S`: Distance of satellite from center of Earth.
* :math:`\gamma_I`: Instrument look angle. 
* :math:`R_n`: Slant-range to near edge of swath.
* :math:`R_f`: Slant-range to far edge of swath.
* :math:`\gamma_n`: Look angle to nearest (to the satellite) part of swath.
* :math:`\gamma_f`: Look angle to farthest (to the satellite) part of swath.
* :math:`\theta_{in}`: Incidence angle to nearest (to the satellite) part of swath.
* :math:`\theta_{if}`: Incidence angle to farthest (to the satellite) part of swath.
* :math:`\theta_{im}`: Incidence angle at ground corresponding to the instrument look-angle (~middle of swath).
* :math:`\alpha_n`: Core angle of nearest part of swath.
* :math:`\alpha_f`: Core angle of farthest part of swath.
* :math:`\alpha_m`: Core angle corresponding to the instrument look-angle (~middle of swath).
* :math:`W_{gr}`: Swath-width in case of *FIXED* and *CROSS_TRACK* scans.
* :math:`\gamma^{cs}`: *CONICAL* scan off-nadir angle (= look angle to the scanned strip).
* :math:`\theta_i^{cs}`: Incidence angle to the *CONICAL* scan strip.
* :math:`\alpha^{cs}`: Earth centric angle (angle b/w the nadir position to the scanned strip about center of Earth) in *CONICAL* scan swath calculations.
* :math:`r^{cs}`: Small circle (on Earth) radius in *CONICAL* scan swath calculations.
* :math:`A^{cs}`: Scanned arc length (*CONICAL* scan).
