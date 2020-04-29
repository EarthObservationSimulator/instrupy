Passive Optical Scanner Description
************************************

References:

1. James R. Wertz and  Wiley J. Larson  (editors), *Space Mission Analysis and Design*, 3rd edition, chapter 9. 


Input JSON format specifications description
===============================================

.. csv-table:: Input parameter description 
   :header: Parameter, Data type, Units, Description
   :widths: 10,10,5,40

   @type, string, ,Must be *Passive Optical Scanner*
   @id, string, , Unique identifier for the instrument.
   name, string, ,Full name of the instrument 
   acronym, string, ,Acronym or initialism or abbreviation.
   mass, float, kilograms,Total mass of this entity.
   volume, float, :code:`m^3`,Total volume of this entity.
   power, float, Watts, Nominal operating power.
   orientation, :ref:`orientation_json_obj`, ,Orientation of the instrument with respect to Nadir-frame.
   fieldOfView, :ref:`fieldOfView_json_obj`, ,Field of view specification of instrument. Only field of view of :code:`"sensorGeometry": "RECTANGULAR"` is accepted.
   numStripsInScene, int, , Number of consequetive scanned strips in a scene See :ref:`ifov_fov_scenefov_for_desc`.
   dataRate, float, Mega-bits per s, Rate of data recorded during nominal operations.
   scanTechnique, string, ,Accepted values are ":code:`PUSHBROOM`" or ":code:`WHISKBROOM`" or ":code:`MATRIX_IMAGER`".
   numberDetectorRowsAT, integer, ,Number of detector rows in along-track direction.
   numberDetectorColsCT, integer, ,Number of detector columns in cross-track direction.
   Fnum, float, ,F-number/ F# of lens.
   focalLength, float, meters, Focal length of lens.
   operatingWavelength, float, meters, Center operating wavelength.
   bandwidth, float, meters, Bandwidth of operation.
   quantumEff, float, , Quantum efficiency of the detector element (:math:`0 < QE < 1`).
   numOfReadOutE, float, , Number of read out electrons of detector.
   targetBlackBodyTemp, float, Kelvin, Target equivalent black-body temperature.
   bitsPerPixel, integer, ,Bits encoded per pixel of image.
   detectorWidth, float, meters,Width of detector element.
   apertureDia, float, meters, Telescope aperture diameter.
   maxDetectorExposureTime, float, seconds, maximum exposure time on the detector elements (optional parameter).
   snrThreshold, float,, Threshold value of SNR for observation to be classified as 'Valid'
   considerAtmosLoss, bool,, True/False flag to specify if atmospheric losses should be taken into account using LOWTRAN 3rd party package. Default is `False`.
   minRequiredAccessTime, float, seconds, Minimum required access time over a ground-point for observation to be possible. Required for matrix imagers.
   maneuverability, :ref:`maneuverability_json_object`, ,Payload maneuverability (see :ref:`manuv_desc`)

.. figure:: passive_scanner_aperture_figure.png
   :scale: 75 %
   :align: center

   Diagram of rectangular aperture illustrating the input parameters :code:`numberDetectorRowsAT`, :code:`numberDetectorColsCT` and :code:`detectorWidth`.

.. warning:: Some of the inputs are interdependent. The dependency **must** be satisfied by the values input by the user.
             The present version of the instrupy package does **not** check for the consistency of the values.

             Following relations between the inputs must be satisfied:

             *  Only square detectors are supported. Hence the IFOV of the detectors must be equal for the along-track 
                and cross-track directions. This results in following relationship: 

                :math:`IFOV = \dfrac{\theta_{AT}}{N_{pix}^{AT}} = \dfrac{\theta_{CT}}{N_{pix}^{CT}} = \dfrac{d}{f}`

                where,
                :math:`IFOV` is the instantaneous FOV or FOV per detector, 
                :math:`\theta_{AT}` is the along-track (angular) FOV,
                :math:`\theta_{CT}` is the cross-track (angular) FOV,
                :math:`N_{pix}^{AT}` is the number of ground-pixels in along-track direction,
                :math:`N_{pix}^{CT}` is the number of ground-pixels in cross-track direction,
                :math:`d` is detector element length,
                :math:`f` is the focal length.

             *  :math:`F\# = \dfrac{f}{D}`

                where,
                :math:`F\#` is the F-number and :math:`D` is the aperture diameter.

.. warning:: Note there is difference between **"ground-pixel"** and **"detectors"**. Detectors refer to the actual physical discrete sensing elements on the scanner aperture. While ground-pixels refer 
             to the imaged pixels on the ground. The number of detectors in the cross-track direction will be less than the number of ground-pixels in the cross-track direction in case of Whiskbroom scanners.

.. _passive_optical_scanner_data_metrics_calc:

Output observation metrics calculation
========================================================

 .. note:: See :ref:`passive_optical_scanner_glossary` for names of the variables used in any discussion below.

.. csv-table:: Observation data metrics table
    :widths: 8,4,4,20
    :header: Metric/Aux data,Data Type,Units,Description 
                                                                                                                                                                                                  
    Coverage [T/F], string,, Indicates if observation was  possible during the access event  (True/ False).                                                                        
    Noise-Equivalent delta T [K], float, Kelvin  , Noise Equivalent delta temperature. Characterizes the instrument in its ability to resolve temperature variations for a given background temperature. 
    DR, float,, Dynamic Range. Is the quotient of the signal and read-out noise electrons the sensor sees between dark and bright scenes.                            
    SNR, float,, Signal-to-Noise ratio                                                                                                                                 
    Ground Pixel Along-Track  Resolution [m], float, meters, Along-track pixel resolution                                                                                                                          
    Ground Pixel Cross-Track Resolution [m] , float, meters, Cross-track pixel resolution 

Viewing geometry
-----------------

See :ref:`satellite_to_target_viewing_geometry` for the calculation of the viewing sensorGeometry parameters.

Ground-pixel resolution calculations
--------------------------------------
Accurate only when ground-pixel is being imaged at Nadir or at strictly sidelooking geometry to the ground track.

:math:`\xi = \dfrac{d}{f}`

:math:`\rho_{CT} = \xi \dfrac{R}{\cos\theta_i}`

:math:`\rho_{AT} = \xi R`

.. todo:: Update for the general target geometry. 

Integration time calculation
----------------------------- 

Let :math:`t_{acc}` be the total access time of the instrument over a ground-point. It can be calculated analytically as:
      
:math:`t_{acc} = \theta_{AT} \hspace{2mm} h/ v_g`

.. todo:: Update access time calculation for general target geometry. Above formulation is valid only for the Nadir case or for strictly 
          sidelooking geometry.

PUSHBROOM
^^^^^^^^^^^^^^^^^^

.. note:: Only one detector array (in cross-track) supported.

:math:`T_i =  t_{acc}`

WHISKBROOM
^^^^^^^^^^^^^^^^^^

.. note:: Only one detector array (in along-track) supported

:math:`T_i =  \dfrac{t_{acc}  N_{pix}^{AT}}{N_{pix}^{CT}}`

MATRIX_IMAGER
^^^^^^^^^^^^^^^^^^

:math:`T_i =  t_{acc}`

If the calculated integration time is greater than the user-defined maximum detector exposure time, it is set to the user-defined maximum detector exposure
time.

:math:`if \hspace{2mm} T_i > T^{exp}_{max}, T_i =  T^{exp}_{max}`


Calculation of signal electrons
-----------------------------------

.. note:: The units of radiance used is [:math:`photons \hspace{1mm} s^{-1} \hspace{1mm} m^{-2} \hspace{1mm} sr^{-1}`]

Radiance with Earth as blackbody radiator
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Assume Earth (target under observation) is a black-body and a Lambertian surface, i.e. the radiance
is independent of the angle. 

:math:`L_{E} = \int_{\lambda_1}^{\lambda_2} L_{\lambda} \tau_{\lambda}^{atm} \cos\theta_i`

where the spectral radiance is given from Planks blackbody radiation equation,

:math:`L_{\lambda} = \dfrac{2 \Upsilon c^2}{\lambda^5} \dfrac{1}{\exp{\dfrac{\Upsilon c}{\lambda k_B T} - 1}}`


Radiance with Earth as reflector of Solar energy
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Assume Sun is a blackbody with temperature 6000K. Also assumed is that the reflectivity of the Earths surface is unity over all wavelength.

:math:`L_S =  \int_{\lambda_1}^{\lambda_2} L_{\lambda} \tau_{\lambda}^{atm}`

.. note:: :math:`\tau_{\lambda}^{atm}` here considers the two-way atmospheric losses, i.e. Sun to Ground and Ground to Satellite. 
          Strictly speaking the Ground to Satellite atmospheric loss appears later, but mathematically either way the result
          is the same. In the present implementation framework it is easier to consider the term here since after this stage
          of calculation, the spectral information (energy per unit wavelength/frequency) is lost.

:math:`{\bf V_{Sun2T}} = {\bf T} - {\bf P_{Sun}}`

:math:`\theta_i^{Solar} = \cos^{-1}(\dfrac{{\bf T} \cdot -{\bf V_{Sun2T}}}{|{\bf T}||\bf V_{Sun2T}|})`

:math:`L^{dw}_S = L_S  \cos\theta_i^{Solar}`

:math:`A_{gp} = \rho_{CT} \rho_{AT}`

:math:`R^{dw}_S|_{ph} = L^{dw}_S A_{gp} \dfrac{\pi r_{Solar}^2}{|{\bf V_{Sun2T}}|^2}`
        
:math:`R^{uw}_S|_{ph} = R^{dw}_S|_{ph} \cos\theta_i` 

:math:`L^{uw}_S = \dfrac{R^{uw}_S|_{ph}}{4 \pi A_{gp}}`
 
Radiance to Signal electrons calculation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

:math:`L_T = L_{E} + L^{uw}_S`

:math:`R^{rad}_T|_{ph} = L_T A_{gp}`

:math:`R^{sen}_T|_{ph} = \dfrac{R^{rad}_T|_{ph}}{|{\bf R}|^2} (\dfrac{D_{ap}}{2})^2 \pi`

:math:`R^{det}_T|_{ph} = R^{sen}_T|_{ph} \tau_{op}`

:math:`N_{ph} = R^{det}_T|_{ph} T_i`

:math:`N_e = N_{ph} Q_E`



Calculation of signal-to-noise-ratio
---------------------------------------

:math:`N_{sh} = \sqrt{N_e}`

:math:`N_t = \sqrt{N_n^2 + N_r^2}`

:math:`SNR = \dfrac{N_e}{N_t}`

Calculation of dynamic range
-----------------------------------

:math:`DR = \dfrac{N_e}{N_r}`

Calculation of Noise-Equivalent Delta T
----------------------------------------

Calculate number of signal electrons for a 1K raise in the temperature of observation pixel.

:math:`\Delta N = N_{e,new} - N_e`

:math:`NE\Delta T = \dfrac{N_e}{\Delta N}`


.. _passive_optical_scanner_glossary:


Glossary
==========

* :math:`\mathbf{S}`: Position vector of the satellite in the ECI frame (equatorial-plane)
* :math:`\mathbf{T}`: Position vector of the Target ground-point in the ECI frame (equatorial-plane)
* :math:`\mathbf{R}`: Range vector from satellite to target ground point
* :math:`\gamma`:  Look-angle to target ground point from satellite
* :math:`\theta_i`: Incidence angle at the target ground point
* :math:`h`: altitude of satellite
* :math:`v_g`: Ground speed of satellite
* :math:`\xi`: The instantaneous field-of-view / field-of-view of detector
* :math:`d`: Detector width/ length (only square detectors allowed)
* :math:`f`: Focal-length of lens
* :math:`\rho_{CT}`: Cross-track ground-pixel resolution
* :math:`\rho_{AT}`: Along-track ground-pixel resolution
* :math:`T_i`: Integration time of ground-pixel
* :math:`T^{exp}_{max}`: Maximum exposure time on detector
* :math:`t_{acc}`: Access time over the ground-point
* :math:`\theta_{AT}`: Along-track FOV
* :math:`\theta_{CT}`: Cross-track FOV
* :math:`N_{pix}^{AT}`: Number of ground-pixels in along-track direction
* :math:`N_{pix}^{CT}`: Number of ground-pixels in cross-track direction
* :math:`L_{\lambda}`: Plancks spectral blackbody radiance
* :math:`\tau_{\lambda}^{atm}`: Wavelength dependent atmospheric loss (Target to Space) as computed by the software `LowTran-7`
* :math:`L_{E}`: Radiance from Earth in the direction of target ground-pixel.
* :math:`\lambda_{op}`: Operating center wavelength
* :math:`\lambda_1`: Lower end wavelength of operating band
* :math:`\lambda_2`: Upper end wavelength of operating band
* :math:`\Upsilon`: Planks constant
* :math:`T`: Target equivalent blackbody temperature
* :math:`k_B`: Boltzmann constant
* :math:`\lambda`: wavelength
* :math:`{\bf P_{Sun}}`: Position vector of Sun in ECI frame (equatorial-plane)
* :math:`L_S`: The radiance from the Sun
* :math:`{\bf V_{Sun2T}}`: Vector from Sun to Target in ECI frame (equatorial-plane)
* :math:`\theta_i^{Solar}`: Solar incidence angle at ground-pixel
* :math:`A_{gp}`: Observation ground pixel area
* :math:`L^{dw}_S`: Downwelling radiance at target observation ground-pixel
* :math:`R^{dw}_S|_{ph}`: Downwelling photon rate at observation ground-pixel
* :math:`R^{uw}_S|_{ph}`: Upwelling photon rate from the ground-pixel to the observing satellite
* :math:`L^{uw}_S`: Upwelling reflected Solar radiance from the ground-pixel
* :math:`L_T`: Total radiance from the target area
* :math:`R^{rad}_T|_{ph}`: Rate of photons radiated, reflected
* :math:`R^{sen}_T|_{ph}`: Rate of photons at sensor aperture
* :math:`R^{det}_T|_{ph}`: Rate of photons at detector
* :math:`N_{ph}`: Number of photons at the detector
* :math:`N_e`: Number of electrons at the detector
* :math:`Q_E`: Quantum efficiency of detector
* :math:`N_{sh}`: Number of Shott noise electrons
* :math:`N_r`: Number of read out noise electrons 
* :math:`N_{t}`: Total number of noise electrons
* :math:`N_{e,new}`: Number of signal electrons for 1K raise in temperature of observation ground pixel 
* :math:`\Delta N`: Change in number of charge carriers for 1K temperature change
* :math:`NE\Delta T`: Noise equivalent delta temperature difference
* :math:`r_{Solar}`: Solar radius

