Synthetic Aperture Radar Description
*************************************

References:

1. *Performance Limits for Synthetic Aperture Radar - second edition SANDIA Report 2006.* ----> Main reference.
2. *Spaceborne SAR Study: LDRD 92 Final Report SANDIA Report March 1993.* ----> Reference for PRF validity calculations, corrections for spaceborne radar.

 .. note:: See :ref:`synthetic_aperture_radar_glossary` for names of the variables used in any discussion below.


Input JSON format specifications description
===============================================

.. csv-table:: Input parameter description 
   :header: Parameter, Data type,Units,Description
   :widths: 10,10,8,40

   name, string, ,Full name of the instrument 
   acronym, string, ,Acronym or initialism or abbreviation.
   mass, float, kilograms,Total mass of this entity.
   volume, float, :code:`m^3`,Total volume of this entity.
   power, float, Watts, Nominal operating power.
   orientation, :ref:`orientation_json_string`, ,Orientation of the instrument with respect to Nadir-frame. Only orientation of :code:`"convention": "SIDE_LOOK"` is accepted.
   dataRate, float, Megabits per s,Rate of data recorded during nominal operations.
   bitsPerPixel, integer, ,Bits encoded per pixel of image.
   pulseWidth, float, seconds, Actual pulse width.
   antennaAlongTrackDim, float, meters, Antenna size in the along-track direction.
   antennaCrossTrackDim, float, meters, Antenna size in the cross-track direction.
   antennaApertureEfficiency, float, ,Aperture efficiency of antenna (:math:`0 < \eta_{ap} < 1`).
   operatingFrequency, float, Hertz, Operating radar center frequency.
   peakTransmitPower, float, Watts, Peak transmit power.
   chirpBandwidth, float, Hertz, Bandwidth of radar operation.
   minimumPRF, float, Hertz, The minimum pulse-repetition-frequency of operation.
   maximumPRF, float,  Hertz, The maximum pulse-repetition-frequency of operation.
   sceneNoiseTemp, float, Kelvin, Nominal scene noise temperature.
   systemNoiseFigure, float, decibels, System noise figure for the receiver. 
   radarLosses, float, decibels, These include a variety of losses primarily over the microwave signal path but doesn't include the atmosphere.
   sigmaNEZ0threshold, float, decibels, The :math:`\sigma_{NEZ0}` threshold for classification as a valid observation.

.. _synthetic_aperture_radar_csv_output:

Synthetic Aperture Radar Level-0 CSV output file description
=============================================================

Description of the header elements:

.. csv-table:: Level-0 output data-metrics description
    :widths: 8,4,4,20
    :header: Metric/Aux data,Data Type,Units,Description
                                                                                                                                                                                                                         
    :code:`Access From [JDUT1]`                      , float   , Julian Date UT1, Access from time
    :code:`Access Duration [s]`                      , float   , seconds , Duration of access
    :code:`POI index`                                , integer ,         , Index of ground-point.                                                                                                                                   
    :code:`Coverage [T/F]`                           , string    ,         , Indicates if observation was  possible during the access event  (True/ False).                                                                           
    :code:`Incidence Angle [deg]`                    , float   , degrees , Incidence angle at target point calculated assuming spherical Earth.                                                                                                                       
    :code:`Swath-Width [m]`                          , float   , meters  , Swath-width of the strip of which  the imaged pixel is part off.                                                                                         
    :code:`Sigma NEZ Nought [dB]`                    , float   , decibels, The backscatter coefficient of a  target for which the signal power level in final image is equal to the noise power level.  **Lesser is better.**       
    :code:`Ground Pixel Along-Track  Resolution [m]` , float   , meters  , Along-track pixel resolution                                                                                                                             
    :code:`Ground Pixel Cross-Track Resolution [m]`  , float   , meters  , Cross-track pixel resolution                                                                                                                             


Example 
-------

.. csv-table:: Basic Sensor typical data metrics example CSV output file
   :header: Access From [JDUT1],Access Duration [s],POI index,Ground Pixel Along-Track Resolution [m],Ground Pixel Cross-Track Resolution [m],Sigma NEZ Nought [dB],Incidence angle [deg],Swath-width [km],Coverage [T/F]
   :widths: 10,10,10,10,10,10,10,10,10

    2458636.099819839,0.430776178837,0,2.2119895193094723,5191.424108924452,-62.1930195781922,0.006617391394361553,44.160809452464996,True
    2458636.168546318,0.430776178837,0,2.2119895204877777,5324.762373152094,-62.30315683467477,0.006451684190133409,44.16080932581713,True
    2458636.237272797,0.430776178837,0,2.2119895216361742,5465.130457250957,-62.41615989259542,0.006285977157196085,44.16080920238835,True
    2458636.305999276,0.430776178837,0,2.2119895227546684,5613.099672872992,-62.53218213714776,0.006120269942440346,44.16080908217726,True

.. _synthetic_aperture_radar_calc:

Typical observation metrics calculation
=========================================

Viewing geometry
-----------------

See :ref:`satellite_to_target_viewing_geometry` for the calculation of the viewing geometry parameters.

Swath-width
------------
.. warning:: While calculating swath width the instrument nominal look angle (not look angle to the target ground-pixel) 
             must be used!!!!          

*See [2] Pg 23 and 24 (Fig. 5.1.3.1)*

:math:`R_S = R_E + h`

:math:`\gamma_I = \theta_{roll}`       

:math:`\gamma_n = \gamma_I - 0.5 \hspace{1mm} \theta_{elv}`

:math:`\gamma_f = \gamma_I  + 0.5 \hspace{1mm} \theta_{elv}`

:math:`\theta_{in} = \sin^{-1}(\sin(\gamma_n) R_S/R_E)`

:math:`\theta_{if} = \sin^{-1}(\sin(\gamma_f) R_S/R_E)`

:math:`\alpha_n = \theta_{in} - \gamma_n`

:math:`\alpha_f = \theta_{if} - \gamma_f`

:math:`\alpha_s = \alpha_f - \alpha_n`

:math:`W_{gr} = R_E \alpha_s`   

Ground pixel resolution calculations
-------------------------------------

From *[1] equations 36, 23* we can get the target ground-pixel range resolution :math:`\rho_y`

:math:`\rho_y = \dfrac{a_{wr} c}{2 B_T \cos\psi_g}`

From *[2] equation (5.3.6.3)* we get the minimum possible azimuth resolution (for strip mapping) of the ground-pixel resolution.

:math:`\rho_a = \dfrac{D_{az}}{2} \dfrac{v_g}{v_s}`

:math:`\sigma_{NEZ0}` calculations
-----------------------------------

:math:`\psi_g = \dfrac{\pi}{2} - \theta_i` 

Use *[1] equation (17)* to find average transmit power :math:`P_{avg}`

:math:`T_{eff} = \tau_p` (approximate effective pulse duration to be actual pulse duration, as in case of matched filter processing)

:math:`d = T_{eff} \hspace{1mm} f_P` 

:math:`P_{avg} = d \hspace{1mm} P_T`

Use *[1] equation 8*, find :math:`G_A`

:math:`A_A = D_{elv} \hspace{1mm} D_{az}`

:math:`G_A = 4 \pi \dfrac{\eta_{ap} A_A}{\lambda^2}`                

*[1] equation 37* we can get the :math:`\sigma_{NEZ0}`

:math:`\sigma_{NEZ0} = \dfrac{265 \pi^3 k T}{c} (R^3  v_s  \cos\psi_g) \dfrac{ B_T F_N L_{radar} L_{atmos}}{P_{avg} G_A^2 \lambda^3} \dfrac{L_r L_a}{a_{wr} a_{wa}}`

:math:`\sigma_{NEZ0},_{dB} = 10 log_{10}\sigma_{NEZ0}`

.. note:: :math:`v_s` is to be used here. See [2] for more explanation.

.. todo:: Write documentation about calculation of image-footprint velocity

Auxillary calculations
=========================================

Field-of-View calculations
---------------------------
The antenna is assumed to be planar with dimensions :math:`D_{az} \hspace{1mm} \times \hspace{1mm} D_{elv}`. The along-track and cross-track 
beamwidth is calculated as: 

:math:`\theta_{az} = \lambda / D_{az}`,     *[1] (eqn 41)*  

:math:`\theta_{elv} = \lambda / D_{elv}`

The along-track and cross-track antenna beamwidths are set to be the along-track and cross-track (full) field-of-view angles,
hence a rectangular field-of-view geometry.

Checking validity of pulse repetition frequency (PRF)
------------------------------------------------------
The user supplies a range of operable PRFs of the SAR instrument. Depending on the orbit conditions (the altitude of satellite
in our case) a usable/ valid PRF has to be selected for target observation. [2] is the primary reference for this formulation, although some errors have been found (and corrected for the current
implementation) in the text. 
The below conditions need to be satisfied:

1. The length of the echo from 3-dB antenna beam illuminated swath is less than inter-pulse period. See [2] Pg 22, 23 and 24.

    :math:`R_n = \sqrt(R_E^2 + R_S^2 - 2 R_E R_S \cos\alpha_n)` 

    :math:`R_f = \sqrt(R_E^2 + R_S^2 - 2 R_E R_S \cos\alpha_f))` 
            
    :math:`\tau_{near} = 2*Rn/c`

    :math:`\tau_{far} = 2*Rf/c` 

    :math:`PRF_{MAX} = 1.0/(2.0*\tau_p + \tau_{far} - \tau_{near})` 

2. The PRFs are high enough to allow for unambiguous detection of doppler shifts.

    :math:`PRF_{MIN} = \dfrac{v_s}{\rho_{a}}` *[2] equation 5.4.4.2*

3. The echos from target doesn't overlap with a transmit pulse (in the future).

    :math:`N = int(f_P \dfrac{2 R_n}{c}) + 1`

    :math:`\dfrac{N-1}{\tau_{near}-\tau_p} < f_P  < \dfrac{N}{\tau_{far} + \tau_p}` *[2] inequality 5.1.4.1*

4. The echo from Nadir (or a previous transmit pulse) doesn't overlap with the desired echo. Nadir echo is very strong
   (even though the antenna gain in the Nadir direction maybe small) since the range to Nadir is small.

    .. warning:: [2] inequality 5.1.5.2 which gives the Nadir interference condition seems wrong. 
                     Refer my notes for the nadir interference condition.             

    :math:`\tau_{nadir} = \dfrac{2 h}{c}`

    :math:`M = int(f_P \dfrac{2 R_f}{c}) + 1`

    :math:`1 <= m <= M`

    :math:`\dfrac{m}{\tau_{near} - \tau_p - \tau_{nadir}} < f_P` (or)
    :math:`f_P< \dfrac{m}{\tau_{far} + \tau_p - \tau_{nadir}}`     
     

Of all the available valid PRFs, the highest PRF is chosen since it improves the :math:`\sigma_{NEZ0}` observation data-metric.
The reason is that the average transmit power increases (since we keep the transmit pulse length constant), and hence the received 
image signal-to-noise-ratio increases.

.. _synthetic_aperture_radar_glossary:

Glossary
==========

.. note:: The same variable names as in the references are followed as much as possible. However it becomes difficult when merging the formulation in
          case of multiple references. 

* :math:`\mathbf{S}`: Position vector of the satellite in the Earth-Centered-Inertial frame (equatorial-plane)
* :math:`\mathbf{T}`: Position vector of the Target ground-point in the Earth-Centered-Inertial  (equatorial-plane)
* :math:`\mathbf{R}`: Range vector from satellite to target ground pixel
* :math:`\gamma`:  Look-angle to target ground pixel from satellite
* :math:`\theta_i`: Incidence angle at the target ground pixel
* :math:`R_E`: Nominal radius of Earth
* :math:`c`: speed of light
* :math:`h`: altitude of satellite
* :math:`D_{az}`: Dimension of antenna in along-track direction
* :math:`D_{elv}`: Dimension of antenna in cross-track direction
* :math:`\lambda`: Operating center wavelength of the radar
* :math:`\theta_{az}`: Beamwidth of antenna in along-track direction
* :math:`\theta_{elv}`: Beamwidth of antenna in cross-track direction
* :math:`\gamma_I`: Instrument look angle 
* :math:`\theta_{roll}`: Roll angle of the instrument, assuming the instrument is aligned to spacecraft body frame, which in turn is aligned to the *nadir-frame*
* :math:`\gamma_n`: Look angle to nearest part of swath
* :math:`\gamma_f`: Look angle to farthest part of swath
* :math:`\theta_{in}`: Incidence angle to nearest part of swath
* :math:`\theta_{if}`: Incidence angle to farthest part of swath
* :math:`\alpha_n`: Core angle of nearest part of swath
* :math:`\alpha_f`: Core angle of farthest part of swath
* :math:`W_{gr}`: Swath-width 
* :math:`\theta_{im}`: Incidence angle to middle of swath
* :math:`\gamma_m`: Look angle to middle of swath
* :math:`\rho_a`: Azimuth resolution
* :math:`\rho_y`: Ground (projected) cross-range resolution
* :math:`\psi_g`: Grazing angle to target ground pixel
* :math:`T_{eff}`: Effective pulse width 
* :math:`f_P`: pulse-repetition-frequency
* :math:`d`: Duty-cycle
* :math:`P_T`: Peak transmit power 
* :math:`P_{avg}`: Average transmit power
* :math:`A_A`: Area of antenna
* :math:`\eta_{ap}`: aperture efficiency of antenna
* :math:`G_A`: Gain of antenna
* :math:`v_s`: Velocity of satellite
* :math:`v_g`: Ground velocity of satellite footprint
* :math:`R_n`: Slant-range to near edge of swath
* :math:`R_f`: Slant-range to far edge of swath
* :math:`\tau_{near}`: Time of return of echo (from transmit time) from the near end of swath
* :math:`\tau_{far}`:  Time of return of echo (from transmit time) from the far end of swath
* :math:`PRF_{MAX}`: Maximum allowable PRF
* :math:`PRF_{MIN}`: Maximum allowable PRF
* :math:`N`: The number of transmit pulses after which echo from desired swath is received
* :math:`\tau_{nadir}`: Time of return of pulse from Nadir
* :math:`M`: Maximum number of transmit pulses after which echo from desired region completes