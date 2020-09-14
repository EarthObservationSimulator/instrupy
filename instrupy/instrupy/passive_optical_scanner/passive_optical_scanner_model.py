""" 
.. module:: passive_optical_scanner_model

:synopsis: *Module to handle passive optical scanners with detectors operating at 
            Visible and near-Visible (IR and UV) wavelengths.*
            
            **Formulation based on Chapter 9 in Space Mission Analysis and Design, 3rd edition.**

"""

import json
import numpy
import copy
import pandas, csv
import sys
from instrupy.util import Entity, EnumEntity, OrientationConvention, Orientation, SensorGeometry, FieldOfView, MathUtilityFunctions, Constants, FileUtilityFunctions

class ScanTech(EnumEntity):
    """Enumeration of recognized passive optical scanner scanning techniques."""
    PUSHBROOM = "PUSHBROOM",
    WHISKBROOM = "WHISKBROOM",
    MATRIX_IMAGER = "MATRIX_IMAGER"

class PassiveOpticalScannerModel(Entity):

    """A passive optical scanner class. Supports following sub-types of passive optical scanners: 
       
       * Rectangular aperture imager with square detector elements. 
       * WHISKBROOM, PUSHBROOM and MATRIX_IMAGER scanning technique supported.       
      
       :ivar name: Full name of the instrument.
       :vartype name: str
        
       :ivar acronym: Acronym, initialism, or abbreviation.
       :vartype acronym: str

       :ivar mass: Total mass (kg) of this entity.
       :vartype mass: float

       :ivar volume: Total volume (m3) of this entity.
       :vartype volume: float
        
       :ivar power: Nominal operating power (W).
       :vartype power: float

       :ivar orientation: Orientation of the instrument with respect to Nadir-frame. 
       :vartype orientation: :class:`instrupy.util.Orientation`
                                
       :ivar fieldOfView: Field of view specification of instrument. 
       :vartype fieldOfView: :class:`instrupy.util.FieldOfView` 

       :ivar sceneFieldOfView: Field of view corresponding to a "scene" captured by the instrument. A scene is made of multiple concatenated strips.
       :vartype sceneFieldOfView: :class:`instrupy.util.FieldOfView`  

       :ivar fieldOfRegard: Field of view calculated taking into account manuverability of the payload.
       :vartype fieldOfRegard: :class:`instrupy.util.FieldOfView`  
       
       :ivar dataRate: Rate of data recorded (Mbps) during nominal operations.
       :vartype dataRate: float  

       :ivar scanTechnique: Accepted values are "PUSHBROOM" or "WHISKBROOM" or "MATRIX_IMAGER".
       :vartype scanTechnique: str
        
       :ivar numberDetectorRowsAT: Number of detector rows in along-track direction
       :vartype numberDetectorRowsAT: int

       :ivar numberDetectorColsCT: Number of detector columns in cross-track direction
       :vartype numberDetectorColsCT: int

       :ivar Fnum: F-number/ F# of lens
       :vartype Fnum: float

       :ivar focalLength: Focal length of lens (m)
       :vartype focalLength: float

       :ivar operatingWavelength: Center operating wavelength (m)
       :vartype operatingWavelength: float

       :ivar bandwidth: Bandwidth of operation (m)

                        .. note:: It is assumed that the detector element supports the entire bandwidth with same quantum efficiency for all wavelengths.
                                  Assumption maybe reasonable for narrow-bandwidths.
       
       :vartype bandwidth: float

       :ivar quantumEff: Quantum efficiency of the detector element (:math:`0 < QE < 1`)
       :vartype quantumEff: float

       :ivar opticsSysEff: Optical systems efficiency (:math:`0 < \\tau_0 < 1`)
       :vartype opticsSysEff: float

       :ivar numOfReadOutE: Number of read out electrons of detector
       :vartype numOfReadOutE: float

       :ivar targetBlackBodyTemp: Target equivalent black-body temperature
       :vartype targetBlackBodyTemp: float

       :ivar bitsPerPixel: Bits encoded per pixel of image
       :vartype bitsPerPixel: int

       :ivar detectorWidth: width of detector element (m)
       :vartype detectorWidth: float

       :ivar apertureDia: telescope aperture diameter (m)
       :vartype apertureDia: float

       :ivar maxDetectorExposureTime: maximum exposure time of detector (s)
       :vartype maxDetectorExposureTime: float

       :ivar snrThreshold: Threshold value of SNR for observation to be classified as 'Valid'
       :vartype snrThreshold: float
       
       :ivar considerAtmosLoss: Flag to specify that atmos loss is to be consdered. 
       :vartype considerAtmosLoss: bool

       .. todo:: Accommodate input of specific detectivity D* for infrared detectors
   
    """

    def __init__(self, name=None, acronym=None, mass=None,
            volume=None, power=None,  orientation = None,
            fieldOfView=None, sceneFieldOfView = None, fieldOfRegard = None, dataRate=None, scanTechnique = None,
            numberDetectorRowsAT=None, numberDetectorColsCT=None, apertureDia = None,
            Fnum = None, focalLength = None, 
            operatingWavelength = None, bandwidth = None, quantumEff = None, 
            opticsSysEff = None, numOfReadOutE = None, targetBlackBodyTemp = None,
            bitsPerPixel = None, detectorWidth = None, maxDetectorExposureTime= None, snrThreshold = None,
            considerAtmosLoss= None, _id=None):
        """Initialization

        """
        self.name = str(name) if name is not None else None
        self.acronym = str(acronym) if acronym is not None else self.name
        self.mass = float(mass) if mass is not None else None
        self.volume = float(volume) if volume is not None else None
        self.power = float(power) if power is not None else power
        self.orientation = copy.deepcopy(orientation) if orientation is not None else Orientation(0,0,0,1,2,3)
        self.fieldOfView = copy.deepcopy(fieldOfView) if fieldOfView is not None else None
        self.sceneFieldOfView = copy.deepcopy(sceneFieldOfView) if sceneFieldOfView is not None else None
        self.fieldOfRegard = copy.deepcopy(fieldOfRegard) if fieldOfRegard is not None else None
        self.dataRate = float(dataRate) if dataRate is not None else None    
        self.scanTechnique = ScanTech.get(scanTechnique) if scanTechnique is not None else None   
        self.numberDetectorRowsAT = int(numberDetectorRowsAT) if numberDetectorRowsAT is not None else None
        self.numberDetectorColsCT = int(numberDetectorColsCT) if numberDetectorColsCT is not None else None
        self.apertureDia = float(apertureDia) if apertureDia is not None else None
        self.Fnum = float(Fnum) if Fnum is not None else None
        self.focalLength = float(focalLength) if focalLength is not None else None 
        self.operatingWavelength = float(operatingWavelength) if operatingWavelength is not None else None 
        self.bandwidth = float(bandwidth) if bandwidth is not None else None 
        self.quantumEff = float(quantumEff) if quantumEff is not None else None 
        self.opticsSysEff = float(opticsSysEff) if opticsSysEff is not None else None 
        self.numOfReadOutE = float(numOfReadOutE) if numOfReadOutE is not None else None 
        self.targetBlackBodyTemp = float(targetBlackBodyTemp) if targetBlackBodyTemp is not None else None 
        self.bitsPerPixel = int(bitsPerPixel) if bitsPerPixel is not None else None 
        self.detectorWidth = float(detectorWidth) if detectorWidth is not None else None
        self.snrThreshold = float(snrThreshold) if snrThreshold is not None else None
        self.maxDetectorExposureTime = float(maxDetectorExposureTime) if maxDetectorExposureTime is not None else None    
        self.considerAtmosLoss = bool(considerAtmosLoss) if considerAtmosLoss is not None else False # Set to False by default

        super(PassiveOpticalScannerModel,self).__init__(_id, "Passive Optical Scanner")

    @staticmethod
    def from_dict(d):
        """ Parses an instrument from a normalized JSON dictionary.
        
            .. warning:: Some of the inputs are interdependent. The dependency **must** be satisfied by the values input by the user.
                         The present version of the instrupy package does **not** check for the consistency of the values.

                         Following relations between the inputs must be satisfied:

                         *   Only square detectors are supported. Hence the IFOV of the detectors must be equal for the along-track 
                             and cross-track directions. This results in following relationship: 

                             :math:`IFOV = \\dfrac{\\theta_{AT}}{N_{pix}^{AT}} = \\dfrac{\\theta_{CT}}{N_{pix}^{CT}} = \\dfrac{d}{f}`

                             where,
                             :math:`IFOV` is the instantaneous FOV or FOV per detector, 
                             :math:`\\theta_{AT}` is the along-track (angular) FOV,
                             :math:`\\theta_{CT}` is the cross-track (angular) FOV,
                             :math:`N_{pix}^{AT}` is the number of ground-pixels in along-track direction,
                             :math:`N_{pix}^{CT}` is the number of ground-pixels in cross-track direction,
                             :math:`d` is detector element length,
                             :math:`f` is the focal length.

                         *   :math:`F\# = \\dfrac{f}{D}`

                             where,
                             :math:`F\#` is the F-number and :math:`D` is the aperture diameter.

            .. todo:: Change behavior to put less burden on user. 
            
        """
        _scan = ScanTech.get(d.get("scanTechnique", None))
        # Only whiskbroom, pushbroom and step-and-stare scan techniques supported.
        if(_scan == "WHISKBROOM") or (_scan == "PUSHBROOM") or (_scan == "MATRIX_IMAGER"):

            # Only rectangular FOV specs supported
            fov_json_str = d.get("fieldOfView", None)
            fov_geometry = SensorGeometry.get(FileUtilityFunctions.from_json(fov_json_str).get("sensorGeometry",None))
            if(fov_geometry == "RECTANGULAR"):

                if(_scan == "PUSHBROOM"):
                    if(d.get("numberDetectorRowsAT", None) != 1):
                        raise Exception("For PUSHBROOM scanning, only 1 detector-row of along-track detectors allowed.")

                if(_scan == "WHISKBROOM"):
                    if(d.get("numberDetectorColsCT", None) != 1):
                        raise Exception("For whiskbroom scanning only one detector-column in cross-track direction is allowed.")  

                # initialize "Scene FOV" if required        
                numStripsInScene = d.get("numStripsInScene", None)
                if(numStripsInScene):                
                    sc_AT_fov_deg = numStripsInScene * float(fov_json_str["alongTrackFieldOfView"]) 
                    sc_CT_fov_deg = float(fov_json_str["crossTrackFieldOfView"]) 
                    sc_fov_json_str = '{ "sensorGeometry": "RECTANGULAR", "alongTrackFieldOfView":' + str(sc_AT_fov_deg)+ ',"crossTrackFieldOfView":' + str(sc_CT_fov_deg) + '}' 
                    scene_fov = FieldOfView.from_json(sc_fov_json_str)
                else:
                    sc_fov_json_str = None
                    scene_fov = None

                # initialize field-of-regard
                if(sc_fov_json_str):
                    fldofreg_str = {**json.loads(sc_fov_json_str),  **{"maneuverability": d.get("maneuverability", None)}}
                else:
                    fldofreg_str = {**fov_json_str, **{"maneuverability": d.get("maneuverability", None)}}

                return PassiveOpticalScannerModel(
                        name = d.get("name", None),
                        acronym = d.get("acronym", None),
                        mass = d.get("mass", None),
                        volume = d.get("volume", None),
                        power = d.get("power", None),
                        orientation = Orientation.from_json(d.get("orientation", None)),
                        fieldOfView = FieldOfView.from_json(d.get("fieldOfView", None)),
                        fieldOfRegard= FieldOfView.from_json(fldofreg_str),
                        sceneFieldOfView = scene_fov,
                        dataRate = d.get("dataRate", None),
                        operatingWavelength = d.get("operatingWavelength", None),
                        bandwidth = d.get("bandwidth", None),
                        quantumEff = d.get("quantumEff", None),
                        targetBlackBodyTemp = d.get("targetBlackBodyTemp", None),
                        bitsPerPixel = d.get("bitsPerPixel", None),
                        scanTechnique = _scan,
                        opticsSysEff = d.get("opticsSysEff", None),
                        numOfReadOutE = d.get("numOfReadOutE", None),
                        numberDetectorRowsAT = d.get("numberDetectorRowsAT", None),
                        numberDetectorColsCT = d.get("numberDetectorColsCT", None),
                        detectorWidth = d.get("detectorWidth", None),
                        focalLength = d.get("focalLength", None),
                        apertureDia = d.get("apertureDia", None),
                        Fnum = d.get("Fnum", None),
                        maxDetectorExposureTime = d.get("maxDetectorExposureTime", None),
                        snrThreshold = d.get("snrThreshold", None),
                        considerAtmosLoss = d.get("considerAtmosLoss", None),
                        _id = d.get("@id", None)
                        )

            else:
                raise Exception("Error message: Only Rectangular FOVs supported for WHISKBROOM and PUSHBROOM scanning techniques.")
        else:
            raise Exception("Error message: Please specify scanning technique as one of WHISKBROOM/ PUSHBROOM/ MATRIX_IMAGER scanning.")

    def calc_typ_data_metrics(self, SpacecraftOrbitState, TargetCoords):
        ''' Calculate typical observation data metrics.

            :param SpacecraftOrbitState: Spacecraft state at the time of observation. This is approximately taken to be the middle (or as close as possible to the middle) of the access interval.

                               Dictionary keys are: 
                               
                               * :code:`Time[JDUT1]` (:class:`float`), Time in Julian Day UT1. Corresponds to the time of observation. 
                               * :code:`x[km]` (:class:`float`), :code:`y[km]` (:class:`float`), :code:`z[km]` (:class:`float`), Cartesian spatial coordinates of satellite in Earth Centered Inertial frame with equatorial plane at the time of observation.
                               * :code:`vx[km/s]` (:class:`float`), :code:`vy[km/s]` (:class:`float`), :code:`vz[km/s]` (:class:`float`), velocity of spacecraft in Earth Centered Inertial frame with equatorial plane at the time of observation.
            :paramtype SpacecraftOrbitState: dict
            
            :param TargetCoords: Location of the observation.

                               Dictionary keys are: 
                                
                               * :code:`Lat [deg]` (:class:`float`), :code:`Lon [deg]` (:class:`float`), indicating the corresponding ground-point accessed (latitude, longitude) in degrees.
            :paramtype TargetCoords: dict

            :returns: Typical calculated observation data metrics.
            
                      Dictionary keys are:  
                    
                      * :code:`Coverage [T/F]` (:class:`bool`) indicating if observation was possible during the access event.
                      * :code:`Noise-Equivalent Delta T [K]` (:class:`float`) Noise-equivalent delta temperature in Kelvin
                      * :code:`DR` (:class:`float`) Dynamic Range
                      * :code:`SNR` (:class:`float`) Signal-to_noise Ratio
                      * :code:`Ground Pixel Along-Track Resolution [m]` (:class:`float`) Spatial resolution of a hypothetical ground-pixel centered about observation point in meters
                      * :code:`Ground Pixel Cross-Track Resolution [m]` (:class:`float`) Spatial resolution of a hypothetical ground-pixel centered about observation point in meters

            :rtype: dict
                     
            .. note:: We differentiate between **access** and **coverage**. **Access** is when the target location
                      falls under the sensor FOV. **Coverage** is when the target location falls under sensor FOV *and* 
                      can be observed.
            
            .. todo:: revise pixel resolution calculations. Current formula works only when sensor has a pure-sidelooking geometry with the ground-point.
                      revise the analytical calculation of the access duration
                        
        '''        
        # Observation time in Julian Day UT1
        tObs_JDUT1 = SpacecraftOrbitState["Time[JDUT1]"]

        # Target position in ECI frame
        TargetPosition_km = MathUtilityFunctions.geo2eci([TargetCoords["Lat [deg]"], TargetCoords["Lon [deg]"], 0.0], tObs_JDUT1)

        # Spacecraft position, velocity in ECI frame
        SpacecraftPosition_km = numpy.array([SpacecraftOrbitState["x[km]"], SpacecraftOrbitState["y[km]"], SpacecraftOrbitState["z[km]"]])  
        SpacecraftVelocity_kmps = [SpacecraftOrbitState["vx[km/s]"], SpacecraftOrbitState["vy[km/s]"], SpacecraftOrbitState["vz[km/s]"]] 
        
        #  Calculate range vector between spacecraft and POI (Target)
        range_vector_km = TargetPosition_km - SpacecraftPosition_km

        alt_km = numpy.linalg.norm(SpacecraftPosition_km) - Constants.radiusOfEarthInKM
        look_angle = numpy.arccos(numpy.dot(MathUtilityFunctions.normalize(range_vector_km), -1*MathUtilityFunctions.normalize(SpacecraftPosition_km)))
        incidence_angle_rad = numpy.arcsin(numpy.sin(look_angle)*(Constants.radiusOfEarthInKM + alt_km)/Constants.radiusOfEarthInKM)
        
        range_vec_norm_km = numpy.linalg.norm(range_vector_km)

        # Calculate FOV of a single detector, i.e. the IFOV
        iFOV_deg = numpy.rad2deg(self.detectorWidth / self.focalLength)
        # Calculate the cross track spatial resolution of the ground-pixel
        pixelSpatialRes_CT_m = numpy.deg2rad(iFOV_deg)*range_vec_norm_km*1.0e3/numpy.cos(incidence_angle_rad)
        # Calculate along-track spatial resolution of the ground-pixel
        pixelSpatialRes_AT_m = numpy.deg2rad(iFOV_deg)*range_vec_norm_km*1.0e3
        pixelArea_m2 = pixelSpatialRes_AT_m * pixelSpatialRes_CT_m                    
        
        # Analytical calculation of the access duration from the satellite altitude (accurate as long as the observation is considered to be made at zero-pitch, and yaw is zero).
        nadir_accessDuration_s = numpy.deg2rad(self.fieldOfView.get_ATCT_fov()[0])*alt_km/ (MathUtilityFunctions.compute_satellite_footprint_speed(SpacecraftPosition_km,SpacecraftVelocity_kmps) *1e-3) # analytical calculation of the access duration
        
        # The analytically caclulated access duration at the nadir is given as input to the integration time caclulations
        Ti_s = PassiveOpticalScannerModel.calculate_integration_time(self.scanTechnique, self.numberDetectorRowsAT, self.numberDetectorColsCT, nadir_accessDuration_s, iFOV_deg, self.maxDetectorExposureTime, self.fieldOfView.get_rectangular_fov_specs_from_custom_fov_specs()[1])
        
        Ne = PassiveOpticalScannerModel.calculate_number_of_signal_electrons(self.operatingWavelength, self.bandwidth, self.targetBlackBodyTemp, 
                                                                       self.apertureDia, self.opticsSysEff, self.quantumEff,  
                                                                       tObs_JDUT1, SpacecraftPosition_km, TargetPosition_km, pixelArea_m2, Ti_s,
                                                                       self.considerAtmosLoss)
   
        # number of noise electrons
        Nn = numpy.sqrt(Ne)
        # total number of noise electrons
        Nt = numpy.sqrt(Nn**2 + self.numOfReadOutE**2)
       
        # signal to noise ratio
        try:
            SNR = Ne/ Nt
        except ZeroDivisionError:
            SNR = numpy.NaN

        # dynamic range of sensor
        try:
            DR = Ne/self.numOfReadOutE
        except ZeroDivisionError:
            DR = numpy.NaN

        # find noise equivalent temperature difference
        Ne_new = PassiveOpticalScannerModel.calculate_number_of_signal_electrons(self.operatingWavelength, self.bandwidth, self.targetBlackBodyTemp + 1,
                                                                           self.apertureDia, self.opticsSysEff, self.quantumEff, 
                                                                           tObs_JDUT1, SpacecraftPosition_km, TargetPosition_km, pixelArea_m2, Ti_s,
                                                                           self.considerAtmosLoss)

        deltaN = Ne_new - Ne


        if deltaN == 0:
            NEdeltaT = numpy.NaN
        else:
            NEdeltaT = Nn/ deltaN        
            

        if(SNR > self.snrThreshold):
            isCovered = True
        else:
            isCovered = False
    
        obsv_metrics = {}
        obsv_metrics["Ground Pixel Along-Track Resolution [m]"] = pixelSpatialRes_AT_m
        obsv_metrics["Ground Pixel Cross-Track Resolution [m]"] = pixelSpatialRes_CT_m
        obsv_metrics["SNR"] = SNR
        obsv_metrics["DR"] = DR
        obsv_metrics["Noise-Equivalent Delta T [K]"] = NEdeltaT
        obsv_metrics["Coverage [T/F]"] = isCovered

        return obsv_metrics

    @staticmethod
    def calculate_integration_time(scanTechnique,numDetRowsAT, numDetColsCT, accessDuration_s, iFOV_deg,  maxDetectorExposureTime = None, crossTrack_fov_deg = None):
        """ Calculate integration time based on scanning method. 

            :param scanTechnique: Scanning technique, can be *PUSHBROOM*, *WHISKBROOM* or  *MATRIX_IMAGER*.
            :paramtype scanTechnique: Str

            :param numDetRowsAT: number of detector rows in along-track direction, laid on the sensor image-plane.
            :paramtype numDetRowsAT: int

            :param numDetColsCT: number of detector columns in cross-track direction, laid on the sensor image-plane.
            :paramtype numDetColsCT: int

            :param accessDuration_s: access duration in seconds during which the image is built
            :paramtype accessDuration_s: float

            :param iFOV_deg: instantaneous (or) the field-of-view per detector element in degrees.
            :paramtype iFOV_deg: float

            :ivar maxDetectorExposureTime: maximum exposure time of detector (s)
            :vartype maxDetectorExposureTime: float

            :param crossTrack_fov_deg: Cross-track field-of-view in degrees (input needed in case *WHISKBROOM* scanning technique is specified).
            :paramtype crossTrack_fov_deg: float

            :returns: integration time in seconds
            :rtype: float

        """
        scanTechnique = ScanTech.get(scanTechnique)
        if(scanTechnique == "PUSHBROOM"):
            if(numDetRowsAT != 1):
                raise Exception("For pushbroom scanning only one detector-row in along-track direction is allowed.")    
            # integration time for pushbroom is same as that of the time taken by the 
            # satellite to go over one groundpixel in the along-track direction 
            # For the case of single-row array of detectors, this is same as the access time.            
            Ti_s = accessDuration_s

        elif(scanTechnique == "WHISKBROOM"):
            if(numDetColsCT != 1):
                raise Exception("For whiskbroom scanning only one detector-column in cross-track direction is allowed.")  
            # For whiskbroom scanning, a single detector is used to scan across many cross-track ground-pixels
            # hence the integration time is a fraction of the access time over a single ground-pixel.
            # Note that the total access time by a row of detectors corresponds to the total along-track FOV of the instrument (not just the AT-FOV of the single row of detectors)
            N_pixel_CT = crossTrack_fov_deg/ iFOV_deg
            Ti_s = accessDuration_s/ N_pixel_CT 

        elif(scanTechnique == "MATRIX_IMAGER"):
            # For MATRIX_IMAGER scanning the integration time at each detector is the total access time over the entire 2D scene imaged.
            # Note that the total access time by a row of detectors corresponds to the total along-track FOV of the instrument (not just the AT-FOV of the single row of detectors)
            Ti_s = accessDuration_s  

        else:
            raise Exception("Unknown scanning technique specified.")

        if(maxDetectorExposureTime): # is a maximum detector exposure time is specified
            if Ti_s > maxDetectorExposureTime:
                Ti_s = maxDetectorExposureTime

        return Ti_s

    @staticmethod
    def calculate_number_of_signal_electrons(opWav_m, bw_m, bbT_K, apDia_m, opTrns, QE, tObs_JDUT1,
                                            obs_pos_km, tar_pos_km, pixelArea_m2, Ti_s, considerAtmosLoss):
        """Calculate the typical number of signal electrons at the detector due to incoming light energy from a hypothetical pixel centered at the target position.
           The light energy from the observed pixel is due to the following two sources:

           1. Radiant energy from the target by considering the target as a black-body with the supplied black-body temperature. 
           2. Reflected solar energy off the target surface. The Sun is considered as a black body at 6000 K, whose energy is reflected
              off the observed pixel area which has an unit reflectivity and Lambertian properties.
             
           See SMAD 3rd edition, Table 9-15 for the forumlation.
           
           :param opWav_m: operating wavelength [m]
           :paramtype opWav_m: float

           :param bw_m: bandwidth [m]
           :paramtype bw_m: float

           :param bbT_K: black-body temperature [K]
           :paramtype bbT_K: float

           :param apDia_m: Lens aperture diameter [m]
           :paramtype apDia_m: float

           :param opTrns: optical-transmission factor (:math:`0 < \\tau_0  < 1`)
           :paramtype opTrns: float

           :param QE: Quantum efficiency (:math:`0 < QE  < 1`)
           :paramtype QE: float

           :param tObs_JDUT1: observation time [Julian Day UT1]
           :paramtype tObs_JDUT1: float

           :param obs_pos_km: observer position [km]. Must be referenced to the center of Earth.
           :paramtype obs_pos_km: list, float

           :param tar_pos_km: target position [km]. Must be referenced to the center of Earth.
           :paramtype tar_pos_km: list, float

           :param pixelArea_m2: area of observed pixel [m2]
           :paramtype pixelArea_m2: float

           :param Ti_s: integration time [s]
           :paramtype Ti_s: float

           :param considerAtmosLoss: Flag to specify that atmos loss is to be consdered. 
           :paramtype considerAtmosLoss: bool

           :return: Number of electrons
           :rtype: float

        """
        obs2tar_vec_km = numpy.array(tar_pos_km) - numpy.array(obs_pos_km)
        distance_km = numpy.linalg.norm(obs2tar_vec_km)
        alt_km = numpy.linalg.norm(obs_pos_km) - Constants.radiusOfEarthInKM
        lookAngle_rad = numpy.arccos(numpy.dot(MathUtilityFunctions.normalize(obs2tar_vec_km), -1*MathUtilityFunctions.normalize(obs_pos_km)))
        obsIncAng_rad = numpy.arcsin(numpy.sin(lookAngle_rad)*(Constants.radiusOfEarthInKM + alt_km)/Constants.radiusOfEarthInKM)
 
        # estimate total radiance from the surface [photons/s/m2/sr] to the direction of the observer                     
        Lint = PassiveOpticalScannerModel.radianceWithEarthAsBlackBodyRadiator(opWav_m, bw_m, bbT_K, obsIncAng_rad, considerAtmosLoss) + PassiveOpticalScannerModel.radianceWithEarthAsReflector(opWav_m, bw_m, tObs_JDUT1, obs_pos_km, tar_pos_km, pixelArea_m2, considerAtmosLoss)

        # total radiated, reflected photon rate
        L = Lint * pixelArea_m2
        
        # input photons rate at the sensor         
        Pin = (L/((distance_km*1e3)**2)) * ((apDia_m / 2)**2) *numpy.pi
 
        # input photon rate at the detector
        PD = Pin * opTrns
        # input number of photons at the detector
        E = PD * Ti_s
        Np = E # since the unit of E are number of photons

        # number of electrons available
        Ne = Np *QE
        return Ne

    @staticmethod
    def radianceWithEarthAsBlackBodyRadiator(opWav_m, bw_m, bbT_K, obsIncAng_rad, considerAtmosLoss):
        """ Determine radiance from Earth as blackbody radiator in bandwidth of interest.
            Assumes the Earth as blackbody with user-input equivalent black-body temperature.
            Also assumes Earth as Lambertian surface obeying Lambert's coside law.

            :param opWav_m: operating wavelength [m]
            :paramtype opWav_m: float

            :param bw_m: bandwidth [m]
            :paramtype bw_m: float

            :param bbT_K: black-body temperature [K]
            :paramtype bbT_K: float

            :param obsIncAng_rad: observation incidence angle in radians
            :paramtype obsIncAng_rad: float

            :param considerAtmosLoss: Flag to specify if atmospheric loss is to be considered. 
            :paramtype considerAtmosLoss: bool
    
            :return: radiance from Earth [photons/s/m2/sr]
            :rtype: float            
          
        """
        # Ensure that observation incidence angle is input in an positive angular domain
        while(obsIncAng_rad < 0):
            obsIncAng_rad = obsIncAng_rad + (2 * numpy.pi )
        obsIncAng_rad = obsIncAng_rad % (numpy.pi*2)
        
        # Observation incidence angle must be in range 270 deg to 90 deg deg, else line-of-sight doesn't exist.
        if((obsIncAng_rad > numpy.pi/2) and (obsIncAng_rad < (3/2) * numpy.pi)):
            raise Exception("Observation incidence angle should be in range -90 deg to 90 deg.")
        
        Lint_ER = PassiveOpticalScannerModel.planck_photon_integral_with_wavlen_dependent_atmos_loss_1((opWav_m - bw_m*0.5), (opWav_m + bw_m*0.5), bbT_K, obsIncAng_rad, considerAtmosLoss)
        
        # Assume Lambertian surface obeying Lambert's cosine law.
        Lint_ER = Lint_ER * numpy.cos(obsIncAng_rad)  

        return Lint_ER 

    @staticmethod
    def radianceWithEarthAsReflector(opWav_m, bw_m, tObs_JDUT1, obs_pos_km, tar_pos_km, obs_area_m2, considerAtmosLoss):
        """ Determine radiance due to reflection of solar radiation off the Earths surface. Reflectance from the Earth surface is taken as **1**.        

            :param opWav_m: operating wavelength [m]
            :paramtype opWav_m: float

            :param bw_m: bandwidth [m]
            :paramtype bw_m: float

            :param tObs_JDUT1: observation time [Julian Day UT1]
            :paramtype tObs_JDUT1: float

            :param obs_pos_km: observer (satellite) position vector [km]. Must be referenced to ECI (equatorial-plane) frame.
            :paramtype obs_pos_km: list, float

            :param tar_pos_km: target (observed ground pixel) position vector [km]. Must be referenced to ECI (equatorial-plane) frame.
            :paramtype tar_pos_km: list, float

            :param obs_area_m2: area of observed target (m2)
            :paramtype obs_area_m2: float
            
            :param considerAtmosLoss: Flag to specify that atmos loss is to be consdered. 
            :paramtype considerAtmosLoss: bool

            .. todo:: account for wavelength dependent, surface dependent reflectance of Earth
        
        """
        obs2tar_vec_km = numpy.array(tar_pos_km) - numpy.array(obs_pos_km)
        alt_km = numpy.linalg.norm(obs_pos_km) - Constants.radiusOfEarthInKM
        lookAngle_rad = numpy.arccos(numpy.dot(MathUtilityFunctions.normalize(obs2tar_vec_km), -1*MathUtilityFunctions.normalize(obs_pos_km)))
        obsIncAng_rad = numpy.arcsin(numpy.sin(lookAngle_rad)*(Constants.radiusOfEarthInKM + alt_km)/Constants.radiusOfEarthInKM)

        # calculate solar incidence angle
        [solar_inc_angle_rad, solar_distance] = MathUtilityFunctions.compute_sun_zenith(tObs_JDUT1, tar_pos_km)
        if(solar_inc_angle_rad is None) or (solar_inc_angle_rad >= numpy.pi/2): # check if Sun is below horizon
            return 1e-19 # return small number
        ''' 
        Calculate radiance from Sun [photons/s/m2/sr]. Note that here the two-way atmospheric losses (Sun to Ground to Observer) is taken into account, not just one way (Sun to Ground).
        Strictly speaking just the Sun to Ground atmospheric loss should be taken into account at this stage and a later stage the Ground to Sun atmos
        loss should be taken into account. Mathmatically the below implementation is correct, and is favored since it is easier to code. 
        '''
        Lint = PassiveOpticalScannerModel.planck_photon_integral_with_wavlen_dependent_atmos_loss_2((opWav_m - bw_m*0.5), (opWav_m + bw_m*0.5), Constants.SunBlackBodyTemperature, solar_inc_angle_rad, obsIncAng_rad, considerAtmosLoss)

        # Calculate downwelling radiance assuming Lambertian surface obeying Lambert's cosine law.
        Lint_dwnwell = Lint * numpy.cos(solar_inc_angle_rad)  

        # calculate downwelling photon rate 
        Pin_dwnwell = Lint_dwnwell * obs_area_m2 * (numpy.pi*(Constants.SolarRadius /(solar_distance*1e3))**2) # [photons/s]      
   
        # upwelling photon rate (surface reflectivity = 1) 
        Pin_upwell = Pin_dwnwell

        # calculate the reflected upwelling radiance assuming Lambertian surface with unit reflectivity
        Lint_upWell = Pin_upwell * numpy.cos(obsIncAng_rad)  * 1/ (4*numpy.pi * obs_area_m2) 

        return Lint_upWell # [photons/s/m2/sr]

    @staticmethod
    def planck_photon_integral_with_wavlen_dependent_atmos_loss_1(wav_low_m, wav_high_m, bbT_K, obs_zen_rad, considerAtmosLoss):
        """ Blackbody integral of spectral photon radiance in given bandwidth of interest, taking into account
            atmospheric losses. 

            :param wav_low_m: lower wavelength in bandwidth of interest [m]
            :paramtype wav_low_m: float

            :param wav_high_m: higher wavelength in bandwidth of interest  [m]
            :paramtype wav_high_m: float

            :param bbT_K: black-body temperature [K]
            :paramtype bbT_K: float

            :param obs_zen_rad: zenith angle of observer
            :paramtype obs_zen_rad: float

            :param considerAtmosLoss: Flag to specify that atmos loss is to be consdered. 
            :paramtype considerAtmosLoss: bool

            :return: radiance [photons/s/m2/sr] after consideration of atmospheric losses
            :rtype: float
        
        """    
        Lint = 0 # integrated radiance
        if(considerAtmosLoss is True):
            TR = MathUtilityFunctions.get_transmission_Obs2Space(wav_low_m, wav_high_m, obs_zen_rad)
            trEff = TR['transmission'][0].values[:,0]
            wav_m = TR['wavelength_nm'].values * 1e-9
            for indx in range(0,len(wav_m)-1):
                Lint = Lint + 0.5*(trEff[indx]+trEff[indx+1])*(PassiveOpticalScannerModel.planck_photon_integral(wav_m[indx], bbT_K) - PassiveOpticalScannerModel.planck_photon_integral(wav_m[indx+1], bbT_K))
        else:
            Lint = (PassiveOpticalScannerModel.planck_photon_integral(wav_high_m, bbT_K) - PassiveOpticalScannerModel.planck_photon_integral(wav_low_m, bbT_K))
        return Lint

    @staticmethod
    def planck_photon_integral_with_wavlen_dependent_atmos_loss_2(wav_low_m, wav_high_m, bbT_K, sun_zen_rad, obs_zen_rad, considerAtmosLoss):
        """ Blackbody integral of reflected spectral photon radiance in given bandwidth of interest, taking into account
            atmospheric losses.  The two-way atmospheric losses (Sun to Ground to Observer) is taken into account, not just one way 
            (Sun to Ground).  Strictly speaking just the Sun to Ground atmospheric loss should be taken into account at this stage
            and a later stage the Ground to Sun atmos loss should be taken into account. 
            Mathmatically the below implementation is correct, and is favored since it is easier to code. 

            :param wav_low_m: lower wavelength in bandwidth of interest [m]
            :paramtype wav_low_m: float

            :param wav_high_m: higher wavelength in bandwidth of interest  [m]
            :paramtype wav_high_m: float

            :param bbT_K: black-body temperature [K]
            :paramtype bbT_K: float

            :param sun_zen_rad: zenith angle of Sun
            :paramtype sun_zen_rad: float

            :param obs_zen_rad: zenith angle of observer
            :paramtype obs_zen_rad: float

            :param considerAtmosLoss: Flag to specify that atmos loss is to be consdered. 
            :paramtype considerAtmosLoss: bool

            :return: radiance [photons/s/m2/sr] after consideration of atmospheric losses
            :rtype: float
        
        """          
        Lint = 0 # integrated radiance
        if(considerAtmosLoss is True):
            TR_SunPath = MathUtilityFunctions.get_transmission_Obs2Space(wav_low_m, wav_high_m, sun_zen_rad)
            TR_ObsPath = MathUtilityFunctions.get_transmission_Obs2Space(wav_low_m, wav_high_m, obs_zen_rad)
            
            trEff_SunPath = TR_SunPath['transmission'][0].values[:,0]
            trEff_ObsPath = TR_ObsPath['transmission'][0].values[:,0]
            wav_m = TR_SunPath['wavelength_nm'].values * 1e-9
            for indx in range(0,len(wav_m)-1):
                Lint = Lint + (0.5*(trEff_SunPath[indx]+trEff_SunPath[indx+1]))*(0.5*(trEff_ObsPath[indx]+trEff_ObsPath[indx+1]))*(PassiveOpticalScannerModel.planck_photon_integral(wav_m[indx], bbT_K) - PassiveOpticalScannerModel.planck_photon_integral(wav_m[indx+1], bbT_K))
        else:
            Lint = (PassiveOpticalScannerModel.planck_photon_integral(wav_high_m, bbT_K) - PassiveOpticalScannerModel.planck_photon_integral(wav_low_m, bbT_K))

        return Lint


    
    @staticmethod
    def planck_photon_integral(lambda_m, T_K):
        """ Blackbody integral of spectral photon radiance from  
            Algorithm from `Spectral Calc <http://www.spectralcalc.com/blackbody/CalculatingBlackbodyRadianceV2.pdf>`_
            Follows Widger and Woodall, Bulletin of the American Meteorological Society, Vol. 57, No. 10, pp. 1217

            :param lambda_m: wavelength (m)
            :paramtype lambda_m: float

            :param T_K: Equivalent blackbody absolute temperature of body (Kelvin).

            :returns: integrated radiance from supplied wavelength to zero-wavelegnth (photons/s/m2/sr)
                      
                      .. warning:: note that it is **not** from supplied wavelength to infinity wavelength.

            :rtype: float
        
        """
        try:
            sigma_cm = 1/(lambda_m*100)
        except ZeroDivisionError:
            sigma_cm = 1e25 # very high number

        Planck = Constants.Planck
        Boltzmann = Constants.Boltzmann
        Speed_of_light = Constants.speedOfLight
        # compute powers of x, the dimensionless spectral coordinate
        c1 = Planck*Speed_of_light/Boltzmann
        x = c1*100*sigma_cm/T_K 
        x2 = x * x 

        # decide how many terms of sum are needed
        iterations = 2.0 + 20.0/x 
        if (iterations<512):
              pass 
        else:
             iterations = 512
        iter = int(iterations) 
        # add up terms of sum
        sum = 0 
        for n in range(1,iter): 
            dn = 1.0/n 
            sum += numpy.exp(-n*x) * (x2 + 2.0*(x + dn)*dn)*dn 
       
        # return result, in units of photons/s/m2/sr
        kTohc = Boltzmann*T_K/(Planck*Speed_of_light) 
        c2 = 2.0* pow(kTohc,3)*Speed_of_light 
        radiance = c2 *sum 

        return radiance
    
    

    