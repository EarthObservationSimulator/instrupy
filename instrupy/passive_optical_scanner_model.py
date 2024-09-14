""" 
.. module:: passive_optical_scanner_model

:synopsis: *Module to handle passive optical scanners with detectors operating at 
            Visible and near-Visible (IR and UV) wavelengths.*
            
            Formulation based on Chapter 9 in Space Mission Analysis and Design, 3rd edition. 
            
            Also see (for optical instruments): V. Ravindra and S. Nag, "Instrument Data Metrics Evaluator for Tradespace Analysis of Earth Observing Constellations," 
            2020 IEEE Aerospace Conference, Big Sky, MT, USA, 2020.

"""
import numpy as np
import copy
import uuid
from instrupy.util import Entity, EnumEntity, Orientation, SphericalGeometry, ViewGeometry, Maneuver, MathUtilityFunctions, GeoUtilityFunctions, Constants, FileUtilityFunctions

class ScanTech(EnumEntity):
    """Enumeration of recognized scanning techniques.
    
    :cvar PUSHBROOM: Indicates pushbroom scanners such as the Landsat-8 TIRS, OLI instruments, which consist of a linear array of square detectors along the cross-track direction.
    :vartype PUSHBROOM: str

    :cvar WHISKBROOM: Indicates whiskbroom (multielement) scanners such as the MODIS radiometer, which consists of linear array of square detectors in along-track direction with movable mirror assembly.
    :vartype WHISKBROOM: str

    :cvar MATRIX_IMAGER: Indicates matrix imager (also called as the step-and-stare imager) with 2D array of square detectors.
    :vartype MATRIX_IMAGER: str
    
    """
    PUSHBROOM = "PUSHBROOM",
    WHISKBROOM = "WHISKBROOM",
    MATRIX_IMAGER = "MATRIX_IMAGER"

class AtmosphericLossModel(EnumEntity):
    """Enumeration of recognized atmospheric loss models.
    
    :cvar LOWTRAN7: Indicates the Lowtran atmospheric loss model as implemented in the python package available here: https://pypi.org/project/lowtran/2.2.2a0/.
    :vartype LOWTRAN7: str
    
    """
    LOWTRAN7 = "LOWTRAN7"

class PassiveOpticalScannerModel(Entity):
    """A passive optical scanner class. Supports following types of passive optical scanners: 
       
       * Rectangular aperture imager with square detector elements. 
       * WHISKBROOM, PUSHBROOM and MATRIX_IMAGER scanning techniques.       
      
       :ivar name: Full name of the instrument.
       :vartype name: str

       :ivar mass: Total mass (kg) of this entity.
       :vartype mass: float

       :ivar volume: Total volume (cm3) of this entity.
       :vartype volume: float
        
       :ivar power: Nominal operating power (W).
       :vartype power: float

       :ivar orientation: Orientation of the instrument.
       :vartype orientation: :class:`instrupy.util.Orientation`

       :ivar fieldOfView: Field of view of instrument specification (SphericalGeometry and Orientation).
       :vartype fieldOfView: :class:`instrupy.util.ViewGeometry`

       :ivar sceneFieldOfView: Scene field of view specification (SphericalGeometry and Orientation).
       :vartype fieldOfView: :class:`instrupy.util.ViewGeometry`

       :ivar maneuver: Maneuver specification of the instrument. 
       :vartype maneuver: :class:`instrupy.util.Maneuver`  

       :ivar fieldOfRegard: Field of regard of the instrument taking into account the sceneFOV and manueverability of the satellite-sensor system. 
                            Note that this shall be a list in order to accommodate non-intersecting view-geometries. 
                            TODO: Modify behavior to have FOR = FOV when no maneuver is specified (hence fixed pointing).
       :vartype fieldOfRegard: list, :class:`instrupy.util.ViewGeometry`
       
       :ivar pointingOption: List of ``Orientation`` objects which specify the orientationsinto which the instrument-axis can be maneuvered. 
                             The orientation must be specified in the NADIR_POINTING frame.
       :vartype pointingOption: list, :class:`orbitpy.util.Orientation`    
       
       :ivar dataRate: Rate of data recorded (Mbps) during nominal operations.
       :vartype dataRate: float  

       :ivar scanTechnique: Scan technique.
       :vartype scanTechnique: :class:`instrupy.passive_optical_scanner_model.ScanTech`
        
       :ivar numberDetectorRows: Number of detector rows (along the Y-axis of the SENOR_BODY_FIXED frame). If the SENSOR_BODY_FIXED frame is aligned to the NADIR_POINTING frame, this direction corresponds to the along-track direction.
       :vartype numberDetectorRows: int

       :ivar numberDetectorCols: Number of detector columns (along the X-axis of the SENOR_BODY_FIXED frame). If the SENSOR_BODY_FIXED frame is aligned to the NADIR_POINTING frame, this direction corresponds to the cross-track direction.
       :vartype numberDetectorCols: int

       :ivar Fnum: F-number/ F# of lens.
       :vartype Fnum: float

       :ivar focalLength: Focal length of lens in meters.
       :vartype focalLength: float

       :ivar operatingWavelength: Center operating wavelength in meters.
       :vartype operatingWavelength: float

       :ivar bandwidth: Bandwidth of operation in meters.

                    .. note:: It is assumed that the detector element supports the entire bandwidth with same quantum efficiency for all wavelengths.
                                Assumption maybe reasonable for narrow-bandwidths.
    
       :vartype bandwidth: float

       :ivar quantumEff: Quantum efficiency of the detector element (:math:`0 < qe < 1`)
       :vartype quantumEff: float

       :ivar opticsSysEff: Optical systems efficiency (:math:`0 < \\tau_0 < 1`)
       :vartype opticsSysEff: float

       :ivar numOfReadOutE: Number of read out electrons of the detector.
       :vartype numOfReadOutE: float

       :ivar targetBlackBodyTemp: Target body's equivalent black-body temperature in Kelvin.
       :vartype targetBlackBodyTemp: float

       :ivar bitsPerPixel: Bits encoded per pixel of image.
       :vartype bitsPerPixel: int

       :ivar detectorWidth: Width of detector element in meters.
       :vartype detectorWidth: float

       :ivar apertureDia: Telescope aperture diameter in meters.
       :vartype apertureDia: float

       :ivar maxDetectorExposureTime: Maximum exposure time of detector in seconds.
       :vartype maxDetectorExposureTime: float
       
       :ivar atmosLossModel: The atmospheric loss model. 'LOWTRAN7' model is supported. If ``None`` the atmospheric losses are not be considered.
       :vartype atmosLossModel: :class:`instrupy.passive_optical_scanner_model.AtmosphericLossModel` or None

       :ivar _id: Unique instrument identifier.
       :vartype _id: str or int

       .. todo:: Accommodate input of specific detectivity D* for infrared detectors
   
    """

    def __init__(self, name=None, mass=None, volume=None, power=None, orientation=None,
            fieldOfViewGeometry=None, sceneFieldOfViewGeometry=None, maneuver=None, pointingOption=None, dataRate=None, 
            scanTechnique=None, numberDetectorRows=None, numberDetectorCols=None, 
            apertureDia=None, Fnum=None, focalLength=None, 
            operatingWavelength=None, bandwidth=None, quantumEff=None, 
            opticsSysEff=None, numOfReadOutE=None, targetBlackBodyTemp=None,
            bitsPerPixel=None, detectorWidth=None, maxDetectorExposureTime=None,
            atmosLossModel=None, _id=None):
        """ Initialization. All except the below two parameters have identical description as that of the corresponding class instance variables.

        :param fieldOfViewGeometry: Instrument field-of-view spherical geometry.
        :paramtype fieldOfViewGeometry: :class:`instrupy.util.SphericalGeometry`

        :param sceneFieldOfViewGeometry: Scene field-of-view spherical geometry.
        :paramtype sceneFieldOfViewGeometry: :class:`instrupy.util.SphericalGeometry`

        :return: None
        :rtype: None

        """
        self.name = str(name) if name is not None else None
        self.mass = float(mass) if mass is not None else None
        self.volume = float(volume) if volume is not None else None
        self.power = float(power) if power is not None else power
        self.orientation = copy.deepcopy(orientation) if orientation is not None else None
        self.fieldOfView = ViewGeometry(orien=self.orientation, sph_geom=fieldOfViewGeometry) if self.orientation is not None and fieldOfViewGeometry is not None and isinstance(fieldOfViewGeometry, SphericalGeometry) else None
        self.sceneFieldOfView = ViewGeometry(orien=self.orientation, sph_geom=sceneFieldOfViewGeometry) if self.orientation is not None and sceneFieldOfViewGeometry is not None and isinstance(sceneFieldOfViewGeometry, SphericalGeometry) else None
        self.maneuver = copy.deepcopy(maneuver) if maneuver is not None and isinstance(maneuver, Maneuver) else None
        self.fieldOfRegard = self.maneuver.calc_field_of_regard(self.sceneFieldOfView.sph_geom) if (self.maneuver is not None and self.sceneFieldOfView is not None) else None
        self.pointingOption = None
        if isinstance(pointingOption, list):
            if all(isinstance(x, Orientation) for x in pointingOption):
                self.pointingOption = pointingOption 
        elif isinstance(pointingOption, Orientation):
            self.pointingOption = [pointingOption] # make into single-element lis
        self.dataRate = float(dataRate) if dataRate is not None else None    
        self.scanTechnique = ScanTech.get(scanTechnique) if scanTechnique is not None else None   
        self.numberDetectorRows = int(numberDetectorRows) if numberDetectorRows is not None else None
        self.numberDetectorCols = int(numberDetectorCols) if numberDetectorCols is not None else None
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
        self.maxDetectorExposureTime = float(maxDetectorExposureTime) if maxDetectorExposureTime is not None else None    
        self.atmosLossModel = AtmosphericLossModel.get(atmosLossModel) if atmosLossModel is not None else None # Set to None by default

        super(PassiveOpticalScannerModel,self).__init__(_id, "Passive Optical Scanner")

    @staticmethod
    def from_dict(d):
        """ Parses a ``PassiveOpticalScannerModel`` object from a normalized JSON dictionary. 
            Refer to :ref:`passive_optical_scanner_model_desc` for description of the accepted key/value pairs.
        
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

        The following default values are assigned to the object instance parameters in case of 
        :class:`None` values or missing key/value pairs in the input dictionary.

        .. csv-table:: Default values
            :header: Parameter, Default Value
            :widths: 10,40

            orientation, Referenced and aligned to the SC_BODY_FIXED frame.
            sceneFieldOfViewGeometry, fieldOfViewGeometry
            atmosLossModel, None (i.e. do not consider atmospheric loss).
            targetBlackBodyTemp, 290 Kelvin
            _id, random string

        :param d: Normalized JSON dictionary with the corresponding model specifications. 
        :paramtype d: dict

        :returns: ``PassiveOpticalScannerModel`` object initialized with the input specifications.
        :rtype: :class:`instrupy.PassiveOpticalScannerModel`
            
        """
        # parse the pointing options as a list of Orientation objects.
        pnt_opt_dict = d.get("pointingOption", None)
        _pointing_option = None
        if pnt_opt_dict:
            # translate to a list of Orientation objects
            if isinstance(pnt_opt_dict, list):
                _pointing_option = [Orientation.from_dict(x) for x in pnt_opt_dict]
            else:
                _pointing_option = [Orientation.from_dict(pnt_opt_dict)]

        _scan = ScanTech.get(d.get("scanTechnique", None))
        # Only whiskbroom, pushbroom and step-and-stare scan techniques supported.
        if(_scan == ScanTech.WHISKBROOM) or (_scan == ScanTech.PUSHBROOM) or (_scan == ScanTech.MATRIX_IMAGER):

            # Only rectangular FOV specs supported
            fov_dict = d.get("fieldOfViewGeometry", None)
            if fov_dict:
                fov = SphericalGeometry.from_dict(fov_dict)
            else:
                raise RuntimeError('Please specify a valid field-of-view geometry specification.')
            if(fov.shape ==SphericalGeometry.Shape.RECTANGULAR):
                if(_scan == "PUSHBROOM"):
                    if(d.get("numberDetectorRows", None) != 1):
                        raise Exception("For PUSHBROOM scanning, only 1 detector-row is allowed.")

                if(_scan == "WHISKBROOM"):
                    if(d.get("numberDetectorCols", None) != 1):
                        raise Exception("For whiskbroom scanning only one detector-column is allowed.")  

                scene_fov_geom = d.get("sceneFieldOfViewGeometry", fov_dict)  # default sceneFOV geometry is the instrument FOV geometry

                return PassiveOpticalScannerModel(
                        name = d.get("name", None),
                        mass = d.get("mass", None),
                        volume = d.get("volume", None),
                        power = d.get("power", None),
                        orientation = Orientation.from_json(d.get("orientation", {"referenceFrame": "SC_BODY_FIXED", "convention": "REF_FRAME_ALIGNED"})), # default orientation is referenced and aligned to to the SENSOR_BODY_FIXED frame.
                        fieldOfViewGeometry = fov,
                        sceneFieldOfViewGeometry = SphericalGeometry.from_json(scene_fov_geom), 
                        maneuver= Maneuver.from_json(d.get("maneuver", None)),
                        pointingOption = _pointing_option,
                        dataRate = d.get("dataRate", None),
                        operatingWavelength = d.get("operatingWavelength", None),
                        bandwidth = d.get("bandwidth", None),
                        quantumEff = d.get("quantumEff", None),
                        targetBlackBodyTemp = d.get("targetBlackBodyTemp", 290), # default is 290 K
                        bitsPerPixel = d.get("bitsPerPixel", None),
                        scanTechnique = _scan,
                        opticsSysEff = d.get("opticsSysEff", None),
                        numOfReadOutE = d.get("numOfReadOutE", None),
                        numberDetectorRows = d.get("numberDetectorRows", None),
                        numberDetectorCols = d.get("numberDetectorCols", None),
                        detectorWidth = d.get("detectorWidth", None),
                        focalLength = d.get("focalLength", None),
                        apertureDia = d.get("apertureDia", None),
                        Fnum = d.get("Fnum", None),
                        maxDetectorExposureTime = d.get("maxDetectorExposureTime", None),
                        atmosLossModel = d.get("atmosLossModel", None),
                        _id = d.get("@id", str(uuid.uuid4()))
                        )

            else:
                raise Exception("Only Rectangular FOV geometries supported for passive-optical scanner instrument type.")
        else:
            raise Exception("Please specify scanning technique as one of WHISKBROOM/ PUSHBROOM/ MATRIX_IMAGER.")

    def to_dict(self):
        """ Translate the ``PassiveOpticalScannerModel`` object to a python dictionary such that it can be uniquely reconstructed back from the dictionary.

        :returns: ``PassiveOpticalScannerModel`` object as a python dictionary.
        :rtype: dict

        """
        fieldOfViewGeometry_dict = self.fieldOfView.sph_geom.to_dict() if self.fieldOfView is not None and isinstance(self.fieldOfView, ViewGeometry) else None
        sceneFieldOfViewGeometry_dict = self.sceneFieldOfView.sph_geom.to_dict() if self.sceneFieldOfView is not None and isinstance(self.sceneFieldOfView, ViewGeometry) else None
        orientation_dict = self.orientation.to_dict() if self.orientation is not None and isinstance(self.orientation, Orientation) else None
        maneuver_dict = self.maneuver.to_dict() if self.maneuver is not None and isinstance(self.maneuver, Maneuver) else None
        pointing_opt_dict = [Orientation.to_dict(x) for x in self.pointingOption] if self.pointingOption is not None else None
        return dict({
                "@type": "Passive Optical Scanner",
                "name":self.name,
                "mass":self.mass,
                "volume":self.volume,
                "power":self.power,
                "fieldOfViewGeometry":fieldOfViewGeometry_dict,
                "sceneFieldOfViewGeometry": sceneFieldOfViewGeometry_dict,
                "orientation":orientation_dict,
                "maneuver":maneuver_dict,
                "pointingOption": pointing_opt_dict,
                "dataRate":self.dataRate,
                "bitsPerPixel": self.bitsPerPixel,
                "operatingWavelength": self.operatingWavelength,
                "bandwidth": self.bandwidth,
                "quantumEff": self.quantumEff,
                "targetBlackBodyTemp": self.targetBlackBodyTemp,
                "scanTechnique": self.scanTechnique,
                "opticsSysEff": self.opticsSysEff,
                "numOfReadOutE": self.numOfReadOutE,
                "numberDetectorRows": self.numberDetectorRows,
                "numberDetectorCols": self.numberDetectorCols,
                "detectorWidth": self.detectorWidth,
                "focalLength": self.focalLength,
                "apertureDia": self.apertureDia,
                "Fnum": self.Fnum,
                "maxDetectorExposureTime": self.maxDetectorExposureTime,
                "atmosLossModel": self.atmosLossModel.value,
                "@id": self._id
                })

    def __repr__(self):
        return "PassiveOpticalScannerModel.from_dict({})".format(self.to_dict())

    def calc_data_metrics(self, sc_orbit_state, target_coords):
        """ Calculate typical observation data metrics.

        :param sc_orbit_state: Spacecraft state at the time of observation.

        Dictionary keys are: 
        
        * :code:`time [JDUT1]` (:class:`float`), Time in Julian Day UT1. Corresponds to the time of observation. 
        * :code:`x [km]` (:class:`float`), :code:`y [km]` (:class:`float`), :code:`z [km]` (:class:`float`), Cartesian coordinates of spacecraft in EARTH_CENTERED_INERTIAL frame at the time of observation.
        * :code:`vx [km/s]` (:class:`float`), :code:`vy [km/s]` (:class:`float`), :code:`vz [km/s]` (:class:`float`), Velocity of spacecraft in EARTH_CENTERED_INERTIAL frame at the time of observation.
        
        :paramtype sc_orbit_state: dict
        
        :param target_coords: Location of the observation.

        Dictionary keys are: 
    
        * :code:`lat [deg]` (:class:`float`), :code:`lon [deg]` (:class:`float`) indicating the corresponding ground-point accessed (latitude, longitude) in degrees.

        :paramtype target_coords: dict

        :returns: Typical calculated observation data metrics.
        
        Dictionary keys are:  
    
        * :code:`noise-equivalent delta T [K]` (:class:`float`) Noise Equivalent delta temperature in Kelvin. Characterizes the instrument in its ability to resolve temperature variations for a given background temperature. 
        * :code:`dynamic range` (:class:`float`) Dynamic Range. Is the quotient of the signal and read-out noise electrons the sensor sees between dark and bright scenes.
        * :code:`SNR` (:class:`float`) Signal-to-noise Ratio.
        * :code:`ground pixel along-track resolution [m]` (:class:`float`) Spatial resolution [meters] of a hypothetical ground-pixel centered about observation point in the along-track direction.
        * :code:`ground pixel cross-track resolution [m]` (:class:`float`) Spatial resolution [meters] of a hypothetical ground-pixel centered about observation point in the cross-track direction.

        :rtype: dict
        
        .. todo:: Revise pixel-resolution and access-time calculations. Current formula is accurate only for nadir-looking or (purely) side-looking instruments.
                    
        """        
        # Observation time in Julian Day UT1
        tObs_JDUT1 = sc_orbit_state["time [JDUT1]"]

        # Target position in ECI frame
        target_pos = GeoUtilityFunctions.geo2eci([target_coords["lat [deg]"], target_coords["lon [deg]"], 0.0], tObs_JDUT1)

        # Spacecraft position, velocity in ECI frame
        sc_pos = np.array([sc_orbit_state["x [km]"], sc_orbit_state["y [km]"], sc_orbit_state["z [km]"]])  
        sc_vel = [sc_orbit_state["vx [km/s]"], sc_orbit_state["vy [km/s]"], sc_orbit_state["vz [km/s]"]] 
        
        #  Calculate range vector between spacecraft and POI (Target)
        range_vector_km = target_pos - sc_pos

        alt_km = np.linalg.norm(sc_pos) - Constants.radiusOfEarthInKM

        # below snippet required to prevent invalid 'x' due to numerical errors
        x = np.dot(MathUtilityFunctions.normalize(range_vector_km), -1*MathUtilityFunctions.normalize(sc_pos))
        if abs(x-1) < 1e-7:
            x=1
        look_angle = np.arccos(x)
        
        incidence_angle_rad = np.arcsin(np.sin(look_angle)*(Constants.radiusOfEarthInKM + alt_km)/Constants.radiusOfEarthInKM)
        
        range_vec_norm_km = np.linalg.norm(range_vector_km)

        # Calculate FOV of a single detector, i.e. the IFOV (instantaneous FOV) as referred in texts.
        iFOV_deg = np.rad2deg(self.detectorWidth / self.focalLength)
        # Calculate the cross track spatial resolution of the ground-pixel
        res_CT_m = np.deg2rad(iFOV_deg)*range_vec_norm_km*1.0e3/np.cos(incidence_angle_rad)
        # Calculate along-track spatial resolution of the ground-pixel
        res_AT_m = np.deg2rad(iFOV_deg)*range_vec_norm_km*1.0e3
        pixel_area_m2 = res_AT_m * res_CT_m                    
        
        # Analytical calculation of the access duration from the satellite altitude (accurate as long as the observation is considered to be made at Nadir or in side-looking geometry).
        nadir_access_duration_s = np.deg2rad(self.fieldOfView.sph_geom.angle_height)*alt_km/ (GeoUtilityFunctions.compute_satellite_footprint_speed(sc_pos,sc_vel) *1e-3) # analytical calculation of the access duration
        
        # The analytically calculated access duration is given as input to the integration time calculations
        Ti_s = PassiveOpticalScannerModel.calculate_integration_time(self.scanTechnique, self.numberDetectorRows, self.numberDetectorCols, nadir_access_duration_s, iFOV_deg, self.maxDetectorExposureTime, self.fieldOfView.sph_geom.angle_width)
        
        Ne = PassiveOpticalScannerModel.calculate_number_of_signal_electrons(self.operatingWavelength, self.bandwidth, self.targetBlackBodyTemp, 
                                                                       self.apertureDia, self.opticsSysEff, self.quantumEff,  
                                                                       tObs_JDUT1, sc_pos, target_pos, pixel_area_m2, Ti_s,
                                                                       self.atmosLossModel)
   
        # number of noise electrons
        Nn = np.sqrt(Ne)
        # total number of noise electrons
        Nt = np.sqrt(Nn**2 + self.numOfReadOutE**2)
       
        # signal to noise ratio
        if Nt==0:
            SNR = np.NaN
        else:
            SNR = Ne/ Nt

        # dynamic range of sensor
        if self.numOfReadOutE==0:
            DR = np.NaN
        else:
            DR = Ne/self.numOfReadOutE

        # find noise equivalent temperature difference
        Ne_new = PassiveOpticalScannerModel.calculate_number_of_signal_electrons(self.operatingWavelength, self.bandwidth, self.targetBlackBodyTemp + 1,
                                                                           self.apertureDia, self.opticsSysEff, self.quantumEff, 
                                                                           tObs_JDUT1, sc_pos, target_pos, pixel_area_m2, Ti_s,
                                                                           self.atmosLossModel)

        deltaN = Ne_new - Ne

        # noise-equivalent delta temperature
        if deltaN==0:
            NEdeltaT = np.inf
        else:    
            NEdeltaT = Nn/ deltaN     
    
        obsv_metrics = {}
        obsv_metrics["ground pixel along-track resolution [m]"] = round(res_AT_m, 2)
        obsv_metrics["ground pixel cross-track resolution [m]"] = round(res_CT_m, 2)
        obsv_metrics["SNR"] = round(SNR, 2)
        obsv_metrics["dynamic range"] = round(DR, 2)
        obsv_metrics["noise-equivalent delta T [K]"] = round(NEdeltaT, 5)

        return obsv_metrics

    @staticmethod
    def calculate_integration_time(scan_tech, num_det_rows, num_det_cols, access_dur_s, iFOV_deg,  max_det_exp_time=None, angle_width_deg=None):
        """ Calculate integration time based on the scanning method. 

        :param scan_tech: Scanning technique.
        :paramtype scan_tech: :class:`instrupy.passive_optical_scanner_model.ScanTech`

        :param num_det_rows: Number of detector rows on the focal-plane-array. If the SENSOR_BODY_FIXED frame is aligned to the NADIR_POINTING frame, this direction corresponds to the along-track direction.
        :paramtype num_det_rows: int

        :param num_det_cols: Number of detector columns on the focal-plane-array. If the SENSOR_BODY_FIXED frame is aligned to the NADIR_POINTING frame, this direction corresponds to the cross-track direction.
        :paramtype num_det_cols: int

        :param access_dur_s: Access duration in seconds during which the image is built.
        :paramtype access_dur_s: float

        :param iFOV_deg: Instantaneous (or) the field-of-view per detector element in degrees.
        :paramtype iFOV_deg: float

        :ivar max_det_exp_time: Maximum exposure time of detector in seconds.
        :vartype max_det_exp_time: float

        :param angle_width_deg: Angular width of the field-of-view in degrees (input needed in case *WHISKBROOM* scanning technique is specified).
        :paramtype angle_width_deg: float

        :returns: Integration time in seconds.
        :rtype: float

        """
        if(scan_tech == ScanTech.PUSHBROOM):
            if(num_det_rows != 1):
                raise Exception("For pushbroom scanning only one detector-row is allowed.")    
            # integration time for pushbroom is same as that of the time taken by the 
            # satellite to go over one groundpixel in the along-track direction 
            # For the case of single-row array of detectors, this is same as the access time.            
            Ti_s = access_dur_s

        elif(scan_tech == ScanTech.WHISKBROOM):
            if(num_det_cols != 1):
                raise Exception("For whiskbroom scanning only one detector-column is allowed.")  
            # For whiskbroom scanning, a single detector is used to scan across many cross-track ground-pixels
            # hence the integration time is a fraction of the access time over a single ground-pixel.
            # Note that the total access time by the column of detectors ( = access_dur_s) corresponds to the total along-track FOV of the instrument.
            N_pixel_CT = angle_width_deg/ iFOV_deg 
            Ti_s = access_dur_s/ N_pixel_CT 

        elif(scan_tech == ScanTech.MATRIX_IMAGER):
            # For MATRIX_IMAGER scanning the integration time at each detector is the total access time over the entire 2D scene imaged.
            # Note that the total access time by the column of detectors ( = access_dur_s) corresponds to the total along-track FOV of the instrument.
            Ti_s = access_dur_s  

        else:
            raise Exception("Unknown scanning technique specified.")

        if(max_det_exp_time is not None): # if a maximum detector exposure time is specified
            if Ti_s > max_det_exp_time:
                Ti_s = max_det_exp_time

        return Ti_s

    @staticmethod
    def calculate_number_of_signal_electrons(op_wav_m, bw_m, bb_T_K, ap_dia_m, op_tx_fac, qe, tObs_JDUT1,
                                            obs_pos_km, tar_pos_km, pixel_area_m2, Ti_s, atmos_loss_model):
        """Calculate the number of signal electrons at the detector due to incoming light energy from a hypothetical pixel centered at the target position.
           The light energy from the observed pixel is due to the following two sources:

           1. Radiant energy from the target by considering the target as a black-body with the supplied black-body temperature. 
           2. Reflected solar energy off the target surface. The Sun is considered as a black body at 6000 K, whose energy is reflected
              off the observed pixel area which has an unit reflectivity and Lambertian properties.
             
           See SMAD 3rd edition, Table 9-15 for the forumlation.
           
           .. note:: The simultaneous assumption of unit reflectance and blackbody to calculate the SNR is contradictory and can be troubling in spectral regions between the thermal and optical. 
                    If the spectral region is in the IR or optical (where the radiative and reflective characteristics of Earth do not
                    intersect), either of the radiated or the reflected energy is dominant, and the contradictory assumption is not troubling.

           :param op_wav_m: Operating wavelength in meters.
           :paramtype op_wav_m: float

           :param bw_m: Bandwidth in meters.
           :paramtype bw_m: float

           :param bb_T_K: Black-body temperature in Kelvin.
           :paramtype bb_T_K: float

           :param ap_dia_m: Lens aperture diameter in meters.
           :paramtype ap_dia_m: float

           :param op_tx_fac: Optical-transmission factor (:math:`0 < \\tau_0  < 1`)
           :paramtype op_tx_fac: float

           :param qe: Quantum efficiency (:math:`0 < qe  < 1`) of the detector.
           :paramtype qe: float

           :param tObs_JDUT1: Observation time [Julian Day UT1].
           :paramtype tObs_JDUT1: float

           :param obs_pos_km: Observer position vector [km]. Must be referenced to the center of Earth.
           :paramtype obs_pos_km: list, float

           :param tar_pos_km: Target position vector [km]. Must be referenced to the center of Earth.
           :paramtype tar_pos_km: list, float

           :param pixel_area_m2: Area of observed pixel [m2].
           :paramtype pixel_area_m2: float

           :param Ti_s: Integration time in seconds.
           :paramtype Ti_s: float

           :param atmos_loss_model: The atmospheric loss model. 'LOWTRAN7' model is supported. If ``None`` the atmospheric losses are not be considered.
           :paramtype atmos_loss_model: :class:`instrupy.passive_optical_scanner_model.AtmosphericLossModel` or None

           :return: Number of signal electrons.
           :rtype: float

        """
        obs2tar_vec_km = np.array(tar_pos_km) - np.array(obs_pos_km)
        distance_km = np.linalg.norm(obs2tar_vec_km)
        alt_km = np.linalg.norm(obs_pos_km) - Constants.radiusOfEarthInKM

        # below snippet required to prevent invalid 'x' due to numerical errors
        x = np.dot(MathUtilityFunctions.normalize(obs2tar_vec_km), -1*MathUtilityFunctions.normalize(obs_pos_km))
        if abs(x-1) < 1e-7:
            x=1
        look_angle = np.arccos(x)

        obs_inc_angle = np.arcsin(np.sin(look_angle)*(Constants.radiusOfEarthInKM + alt_km)/Constants.radiusOfEarthInKM)
 
        # estimate total radiance from the surface [photons/s/m2/sr] to the direction of the observer                     
        Lint = PassiveOpticalScannerModel.radiance_with_earth_as_bb_radiator(op_wav_m, bw_m, bb_T_K, obs_inc_angle, atmos_loss_model) + \
               PassiveOpticalScannerModel.radiance_with_earth_as_reflector(op_wav_m, bw_m, tObs_JDUT1, obs_pos_km, tar_pos_km, pixel_area_m2, atmos_loss_model)

        # total radiated and reflected photon rate
        L = Lint * pixel_area_m2
        
        # input photons rate at the sensor         
        Pin = (L/((distance_km*1e3)**2)) * ((ap_dia_m / 2)**2) *np.pi
 
        # input photon rate at the detector
        PD = Pin * op_tx_fac
        # input number of photons at the detector
        E = PD * Ti_s
        Np = E # since the unit of E are number of photons

        # number of electrons available
        Ne = Np *qe
        return Ne

    @staticmethod
    def radiance_with_earth_as_bb_radiator(wav_m, bw_m, bb_T_K, obs_inc_angle, atmos_loss_model):
        """ Determine radiance from Earth as blackbody radiator in bandwidth of interest.
            Assumes the Earth as blackbody with user-input equivalent black-body temperature.
            Also assumes Earth as Lambertian surface obeying Lambert's cosine law.

        :param wav_m: Wavelength of interest in meters.
        :paramtype wav_m: float

        :param bw_m: Bandwidth in meters.
        :paramtype bw_m: float

        :param bb_T_K: Black-body temperature in Kelvin.
        :paramtype bb_T_K: float

        :param obs_inc_angle: Observation incidence angle in radians.
        :paramtype obs_inc_angle: float

        :param atmos_loss_model: The atmospheric loss model. 'LOWTRAN7' model is supported. If ``None`` the atmospheric losses are not be considered.
        :paramtype atmos_loss_model: :class:`instrupy.passive_optical_scanner_model.AtmosphericLossModel` or ``None``

        :return: radiance from Earth [photons/s/m2/sr]
        :rtype: float            
          
        """
        # Ensure that observation incidence angle is input in an positive angular domain
        while(obs_inc_angle < 0):
            obs_inc_angle = obs_inc_angle + (2 * np.pi )
        obs_inc_angle = obs_inc_angle % (np.pi*2)
        
        # Observation incidence angle must be in range 270 deg to 90 deg deg, else line-of-sight doesn't exist.
        if((obs_inc_angle > np.pi/2) and (obs_inc_angle < (3/2) * np.pi)):
            raise Exception("Observation incidence angle should be in range -90 deg to 90 deg.")
        
        Lint_ER = PassiveOpticalScannerModel.planck_photon_integral_with_wavelen_dependent_atmos_loss_1((wav_m - bw_m*0.5), (wav_m + bw_m*0.5), bb_T_K, obs_inc_angle, atmos_loss_model)
        
        # Assume Lambertian surface obeying Lambert's cosine law.
        Lint_ER = Lint_ER * np.cos(obs_inc_angle)  

        return Lint_ER 

    @staticmethod
    def radiance_with_earth_as_reflector(wav_m, bw_m, tObs_JDUT1, obs_pos_km, tar_pos_km, obs_area_m2, atmos_loss_model):
        """ Determine radiance due to reflection of solar radiation off the Earths surface. Reflectance from the Earth surface is taken to be unity.        

        :param wav_m: Wavelength of interest in meters.
        :paramtype wav_m: float

        :param bw_m: Bandwidth in meters.
        :paramtype bw_m: float

        :param tObs_JDUT1: observation time [Julian Day UT1]
        :paramtype tObs_JDUT1: float

        :param obs_pos_km: observer (satellite) position vector [km]. Must be referenced to EARTH_CENTERED_INERTIAL frame.
        :paramtype obs_pos_km: list, float

        :param tar_pos_km: target (observed ground pixel) position vector [km]. Must be referenced to EARTH_CENTERED_INERTIAL frame.
        :paramtype tar_pos_km: list, float

        :param obs_area_m2: area of observed target (m2)
        :paramtype obs_area_m2: float
        
        :param atmos_loss_model: The atmospheric loss model. 'LOWTRAN7' model is supported. If ``None`` the atmospheric losses are not be considered.
        :paramtype atmos_loss_model: :class:`instrupy.passive_optical_scanner_model.AtmosphericLossModel` or None

        :return: upwelling radiance with Earth as reflector of Solar energy [photons/s/m2/sr]
        :rtype: float
        
        """
        obs2tar_vec_km = np.array(tar_pos_km) - np.array(obs_pos_km)
        alt_km = np.linalg.norm(obs_pos_km) - Constants.radiusOfEarthInKM
        
        # below snippet required to prevent invalid 'x' due to numerical errors
        x = np.dot(MathUtilityFunctions.normalize(obs2tar_vec_km), -1*MathUtilityFunctions.normalize(obs_pos_km))
        if abs(x-1) < 1e-7:
            x=1
        look_angle = np.arccos(x)

        obs_inc_angle = np.arcsin(np.sin(look_angle)*(Constants.radiusOfEarthInKM + alt_km)/Constants.radiusOfEarthInKM)

        # calculate solar incidence angle
        [solar_inc_angle_rad, solar_distance] = GeoUtilityFunctions.compute_sun_zenith(tObs_JDUT1, tar_pos_km)
        if(solar_inc_angle_rad is None) or (solar_inc_angle_rad >= np.pi/2): # check if Sun is below horizon
            return 1e-19 # return small number
        """
        Calculate radiance from Sun [photons/s/m2/sr]. Note that here the two-way atmospheric losses (Sun to Ground to Observer) is taken into account, not just one way (Sun to Ground).
        Strictly speaking just the Sun to Ground atmospheric loss should be taken into account at this stage and a later stage the Ground to Sun atmos
        loss should be taken into account. Mathematically the below implementation is correct, and is favored since it is easier to code. 
        """
        Lint = PassiveOpticalScannerModel.planck_photon_integral_with_wavelen_dependent_atmos_loss_2((wav_m - bw_m*0.5), (wav_m + bw_m*0.5), Constants.SunBlackBodyTemperature, solar_inc_angle_rad, obs_inc_angle, atmos_loss_model)

        # Calculate downwelling radiance assuming Lambertian surface obeying Lambert's cosine law.
        Lint_downwell = Lint * np.cos(solar_inc_angle_rad)  

        # calculate downwelling photon rate 
        Pin_downwell = Lint_downwell * obs_area_m2 * (np.pi*(Constants.SolarRadius /(solar_distance*1e3))**2) # [photons/s]      
   
        # upwelling photon rate (surface reflectivity = 1) 
        Pin_upwell = Pin_downwell

        # calculate the (reflected) upwelling radiance assuming Lambertian surface with unit reflectivity
        Lint_upWell = Pin_upwell * np.cos(obs_inc_angle)  * 1/ (4*np.pi * obs_area_m2) 

        return Lint_upWell # [photons/s/m2/sr]

    @staticmethod
    def planck_photon_integral_with_wavelen_dependent_atmos_loss_1(wav_low_m, wav_high_m, bb_T_K, obs_zen_rad, atmos_loss_model):
        """ Blackbody integral of spectral photon radiance in given bandwidth of interest, taking into account
            atmospheric losses. 
        
        :param wav_low_m: lower wavelength in bandwidth of interest [m]
        :paramtype wav_low_m: float

        :param wav_high_m: higher wavelength in bandwidth of interest  [m]
        :paramtype wav_high_m: float

        :param bb_T_K: black-body temperature [K]
        :paramtype bb_T_K: float

        :param obs_zen_rad: zenith angle of observer
        :paramtype obs_zen_rad: float

        :param atmos_loss_model: The atmospheric loss model. 'LOWTRAN7' model is supported. If ``None`` the atmospheric losses are not be considered.
        :paramtype atmos_loss_model: :class:`instrupy.passive_optical_scanner_model.AtmosphericLossModel` or None

        :return: Radiance [photons/s/m2/sr] after consideration of atmospheric losses.
        :rtype: float
        
        """    
        Lint = 0 # integrated radiance
        if(atmos_loss_model == AtmosphericLossModel.LOWTRAN7):
            TR = GeoUtilityFunctions.get_transmission_Obs2Space(wav_low_m, wav_high_m, obs_zen_rad)
            tr_eff = TR['transmission'][0].values[:,0]
            wav_m = TR['wavelength_nm'].values * 1e-9
            #print(wav_m)
            #print(tr_eff)
            # The array of wavelengths are not necessarily bound within [wav_low_m, wav_high_m]. This could lead to errors. Hence make new modified array of wavelengths. 
            wav_m[0] = wav_high_m
            wav_m[-1] = wav_low_m
            for indx in range(0,len(wav_m)-1):
                Lint = Lint + 0.5*(tr_eff[indx]+tr_eff[indx+1])*(PassiveOpticalScannerModel.planck_photon_integral(wav_m[indx], bb_T_K) - PassiveOpticalScannerModel.planck_photon_integral(wav_m[indx+1], bb_T_K))
        else:
            Lint = (PassiveOpticalScannerModel.planck_photon_integral(wav_high_m, bb_T_K) - PassiveOpticalScannerModel.planck_photon_integral(wav_low_m, bb_T_K))

        return Lint

    @staticmethod
    def planck_photon_integral_with_wavelen_dependent_atmos_loss_2(wav_low_m, wav_high_m, bb_T_K, sun_zen_rad, obs_zen_rad, atmos_loss_model):
        """ Blackbody integral of reflected spectral photon radiance in given bandwidth of interest, taking into account
            atmospheric losses.  The two-way atmospheric losses (Sun to Ground, Ground to Observer) is taken into account, not just one way 
            (Sun to Ground).  Strictly speaking just the Sun to Ground atmospheric loss should be taken into account at this stage
            and at later stage the Ground to Sun atmos loss should be taken into account. 
            Mathematically the below implementation is correct, and is favored since it is easier to code. 

        :param wav_low_m: lower wavelength in bandwidth of interest in meters.
        :paramtype wav_low_m: float

        :param wav_high_m: higher wavelength in bandwidth of interest in meters.
        :paramtype wav_high_m: float

        :param bb_T_K: Black-body temperature in degrees.
        :paramtype bb_T_K: float

        :param sun_zen_rad: Zenith angle of Sun.
        :paramtype sun_zen_rad: float

        :param obs_zen_rad: Zenith angle of observer.
        :paramtype obs_zen_rad: float

        :param atmos_loss_model: The atmospheric loss model. 'LOWTRAN7' model is supported. If ``None`` the atmospheric losses are not be considered.
        :paramtype atmos_loss_model: :class:`instrupy.passive_optical_scanner_model.AtmosphericLossModel` or None

        :return: Radiance [photons/s/m2/sr] after consideration of atmospheric losses.
        :rtype: float
        
        """          
        Lint = 0 # integrated radiance
        if(atmos_loss_model == AtmosphericLossModel.LOWTRAN7):
            TR_SunPath = GeoUtilityFunctions.get_transmission_Obs2Space(wav_low_m, wav_high_m, sun_zen_rad)
            TR_ObsPath = GeoUtilityFunctions.get_transmission_Obs2Space(wav_low_m, wav_high_m, obs_zen_rad)
            
            tr_eff_sun_path = TR_SunPath['transmission'][0].values[:,0]
            tr_eff_obs_path = TR_ObsPath['transmission'][0].values[:,0]
            wav_m = TR_SunPath['wavelength_nm'].values * 1e-9
            for indx in range(0,len(wav_m)-1):
                Lint = Lint + (0.5*(tr_eff_sun_path[indx]+tr_eff_sun_path[indx+1]))*(0.5*(tr_eff_obs_path[indx]+tr_eff_obs_path[indx+1]))*(PassiveOpticalScannerModel.planck_photon_integral(wav_m[indx], bb_T_K) - PassiveOpticalScannerModel.planck_photon_integral(wav_m[indx+1], bb_T_K))
        else:
            Lint = (PassiveOpticalScannerModel.planck_photon_integral(wav_high_m, bb_T_K) - PassiveOpticalScannerModel.planck_photon_integral(wav_low_m, bb_T_K))

        return Lint
    
    @staticmethod
    def planck_photon_integral(lambda_m, T_K):
        """ Blackbody integral of spectral photon radiance from  
            Algorithm from `Spectral Calc <http://www.spectralcalc.com/blackbody/CalculatingBlackbodyRadianceV2.pdf>`_
            Follows Widger and Woodall, Bulletin of the American Meteorological Society, Vol. 57, No. 10, pp. 1217

        :param lambda_m: wavelength (m)
        :paramtype lambda_m: float

        :param T_K: Equivalent blackbody absolute temperature of body (Kelvin).
        :paramtype T_K: float

        :returns: integrated radiance from supplied wavelength to zero-wavelegnth (photons/s/m2/sr)
                    
                    .. warning:: note that integration is **not** from supplied wavelength to infinity wavelength.

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
            sum += np.exp(-n*x) * (x2 + 2.0*(x + dn)*dn)*dn 
       
        # return result, in units of photons/s/m2/sr
        kTohc = Boltzmann*T_K/(Planck*Speed_of_light) 
        c2 = 2.0* pow(kTohc,3)*Speed_of_light 
        radiance = c2 *sum 

        return radiance
    
    def get_id(self):
        """ Get the instrument identifier.

        :returns: instrument identifier.
        :rtype: str

        """
        return self._id
    
    def get_field_of_view(self):
        """ Get the instrument field-of-view.

        :returns: Instrument field-of-view. 
        :rtype: :class:`instrupy.util.ViewGeometry`

        """
        return self.fieldOfView

    def get_scene_field_of_view(self):
        """ Get the scene-field-of-view (sceneFOV).

        :returns: Scene-field-of-view.
        :rtype: :class:`instrupy.util.ViewGeometry`

        """
        return self.sceneFieldOfView

    def get_field_of_regard(self):
        """ Get the instrument field of regard. 

        :returns: Field of regard (list of ``ViewGeometry`` objects). 
        :rtype: list, :class:`instrupy.util.ViewGeometry`

        """
        return self.fieldOfRegard

    def get_orientation(self):
        """ Get the instrument orientation.

        :returns: Instrument orientation.
        :rtype: :class:`instrupy.util.Orientation`

        """
        return self.orientation
    
    def get_pointing_option(self):
        """ Get the list of pointing options.

        :returns: List of pointing options.
        :rtype: list, :class:`instrupy.util.Orientation`

        """
        return self.pointingOption
    
    

    