""" 
.. module:: synthetic_aperture_radar_model

:synopsis: *Module to handle SAR instrument model.*

        Inline code comments contain references to the following articles:

        1. Performance Limits for Synthetic Aperture Radar - second edition SANDIA Report 2006. ----> Main reference.
        2. Spaceborne SAR Study: LDRD 92 Final Report SANDIA Report March 1993. ----> Reference for PRF validity calculations, corrections for spaceborne radar.
        3. V. Ravindra and S. Nag, "Instrument Data Metrics Evaluator for Tradespace Analysis of Earth Observing Constellations," 2020 IEEE Aerospace Conference, Big Sky, MT, USA, 2020, pp. 1-20, doi: 10.1109/AERO47225.2020.9172705.
        
        Polarimetry concepts:
        
        4. *Synthetic Aperture Radar Polarimetry,  Jakob Van Zyl* ----> Reference for compact-pol and AIRSAR implementation for dual-pol.
        5. *SMAP Handbook* ----> Reference for SMAP implementation of dual-pol.
        
        ScanSAR concepts
        
        6. Tomiyasu, Kiyo. "Conceptual performance of a satellite borne, wide swath synthetic aperture radar." IEEE Transactions on Geoscience and Remote Sensing 2 (1981): 108-116.
        7. Moore, Richard K., John P. Claassen, and Y. H. Lin. "Scanning spaceborne synthetic aperture radar with integrated radiometer." IEEE Transactions on Aerospace and Electronic Systems 3 (1981): 410-421.
        8. Currie, A., and Ma A. Brown. "Wide-swath SAR." In IEE Proceedings F (Radar and Signal Processing), vol. 139, no. 2, pp. 122-135. IET Digital Library, 1992.
        9. Chang, Chi-Yung, Michael Y. Jin, Yun-Ling Lou, and Benjamin Holt. "First SIR-C scansar results." IEEE transactions on geoscience and remote sensing 34, no. 5 (1996): 1278-1281.

.. todo:: Include frequency dependent atmospheric losses in :math:`\\sigma_{NESZ}` calculations.

.. todo:: Update antenna specification parameters to the type :class:`instrupy.util.AntennaSpecs`

.. todo:: Make seperate objects for scanning-technique similar to the one in Radiometer model?

"""
import json
import copy
import uuid
import numpy as np
import warnings
from instrupy.util import Entity, EnumEntity, Orientation, SphericalGeometry, ViewGeometry, Maneuver, GeoUtilityFunctions, MathUtilityFunctions, Constants, FileUtilityFunctions

class ScanTech(EnumEntity):
    """Enumeration of recognized SAR scanning techniques.
    
    :cvar STRIPMAP: Stripmap imaging operation.
    :vartype STRIPMAP: str

    :cvar SCANSAR: ScanSAR imaging operation. Multiple strips are scanned in the cross-track direction to increase the overall swath-width 
                   (but resulting in a coarser azimuth resolution due to reduced scan-time per strip). 
    :vartype SCANSAR: str
    
    """
    STRIPMAP = "STRIPMAP",
    SCANSAR = "SCANSAR",

class PolTypeSAR(EnumEntity):
    """Enumeration of recognized SAR polarization types.
    
    :cvar SINGLE: Single transmit and receive polarization.
    :vartype SINGLE: str

    :cvar COMPACT: Single transmit and dual receive polarization.
    :vartype COMPACT: str

    :cvar DUAL: Dual transmit and dual receive polarization.
    :vartype DUAL: str    
    
    """
    SINGLE = "SINGLE",
    COMPACT = "COMPACT",
    DUAL = "DUAL"

class DualPolPulseConfig(EnumEntity):
    """Enumeration of recognized dual-polarization pulse configurations.
    
    :cvar AIRSAR: This pulse configuration is the same as the one implemented by the NASA/JPL AIRSAR systems (see Pg.32, Fig.2-5 in [4]). It consists of transmitting alternating pulses of orthogonal
                  polarization and filtering the received signal into separate orthogonal polarizations.      
    :vartype AIRSAR: str

    :cvar SMAP: This pulse configuration is the same as the one implemented by the SMAP radar (see Pg.41, Fig.26 in [5]). It consists of two slightly separated pulses of 
                orthogonal polarizations at different frequency bands. The received signal is separated into the respective band and the orthogonal 
                polarizations measured. This requires an additional parameter called as the :code:`pulseSeparation` to indicate the separation 
                between the pulses of the two orthogonal polarizations. If not specified a default value of 50% of the pulse-width (:code:`pulseWidth`) is considered.
    :vartype SMAP: str
    
    """
    AIRSAR = "AIRSAR",
    SMAP = "SMAP"

class SwathTypeSAR(EnumEntity):
    """Enumeration of recognized SAR swath imaging configurations.
    
    :cvar FULL: Tne entire illuminated swath by the main-lobe of the antenna is considered.
    :vartype FULL: str

    :cvar FIXED: A fixed swath size (less than the swath illuminated by the main-lobe) is considered. This configuration 
                 could be adapted to meet the Pule Repetition Frequency constraints.
    :vartype FIXED: str

    """
    FULL = "FULL",
    FIXED = "FIXED"

class SyntheticApertureRadarModel(Entity):
    """A synthetic aperture radar class estimating observation data-metrics.      
      
        :cvar L_r: Reduction in SNR gain due to non-ideal range filtering (see [Pg.9, 1]). Default value is 1.2.
        :vartype L_r: float

        :cvar L_a: Reduction in SNR gain due to non-ideal azimuth filtering (see [Pg.10, 1]). Default value is 1.2.
        :vartype L_a: float

        :cvar a_wa:  Azimuth impulse response broadening factor (see [Pg.9, 1]). Default value is 1.2.
        :vartype a_wa: float

        :cvar a_wr: Range impulse response broadening factor (see [Pg.10, 1]). Default value is 1.2.
        :vartype a_wr: float

        :ivar name: Full name of the instrument.
        :vartype name: str

        :ivar mass: Total mass (kg) of this entity.
        :vartype mass: float

        :ivar volume: Total volume (m3) of this entity.
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

        :ivar fieldOfRegard: Field of regard of the instrument taking into account the sensor FOV and manueverability of the satellite-sensor system. 
                             Note that this shall be a list in order to accommodate non-intersecting view-geometries.
        :vartype fieldOfRegard: list, :class:`instrupy.util.ViewGeometry`  
       
        :ivar pointing_option: List of ``Orientation`` objects which specify the orientation of the instrument pointing axis into which the instrument-axis can be maneuvered. 
                               The orientation must be specified in the NADIR_POINTING frame.
        :vartype pointing_option: list, :class:`orbitpy.util.Orientation` 

        :ivar dataRate: Rate of data recorded (Mega-bits-per-sec) during nominal operations.
        :vartype dataRate: float  

        :ivar bitsPerPixel: Number of bits encoded per pixel of image.
        :vartype bitsPerPixel: int    

        :ivar pulseWidth: Actual pulse width in (seconds)  (per channel/polarization).
        :vartype pulseWidth: float

        :ivar antennaHeight: Antenna height (in the along-track direction when SENSOR_BODY_FIXED is aligned to NADIR_POINTING frame).
        :vartype antennaHeight: float

        :ivar antennaWidth: Antenna width (in the cross-track direction when SENSOR_BODY_FIXED is aligned to NADIR_POINTING frame).
        :vartype antennaWidth: float

        :ivar antennaApertureEfficiency: Aperture efficiency of antenna (:math:`0 < \\eta_{ap} < 1`).
        :vartype antennaApertureEfficiency: float

        :ivar operatingFrequency: Operating radar center frequency in (Hertz).
        :vartype operatingFrequency: float

        :ivar peakTransmitPower: Peak transmit power in (Watts).
        :vartype peakTransmitPower: float

        :ivar chirpBandwidth: Chirp bandwidth of radar operation in (Hertz)  (per channel/polarization).
        :vartype chirpBandwidth: float

        :ivar minimumPRF: The minimum allowable pulse-repetition-frequency of operation in (Hertz). 
                          If dual-pol with alternating pol pulses, the PRF specification is considered taking all pulses into account (i.e. is considered as the PRFmaster).
        :vartype minimumPRF: float

        :ivar maximumPRF: The maximum allowable pulse-repetition-frequency of operation in (Hertz).
                          If dual-pol with alternating pol pulses, the PRF specification is considered taking all pulses into account (i.e. is considered as the PRFmaster).
        :vartype maximumPRF: float

        :ivar sceneNoiseTemp: Nominal scene noise temperature in (Kelvin).
        :vartype sceneNoiseTemp: float

        :ivar systemNoiseFigure:  (decibels) System noise figure for the receiver. The system noise figure includes primarily the noise figure of the front-end Low-Noise
                                  Amplifier (LNA) and the losses between the antenna and the LNA. Typical system noise figures for sub-kilowatt radar systems are 3.0 dB to 3.5 dB
                                  at X-band, 3.5 dB to 4.5 dB at Ku-band, and perhaps 6 dB at Ka-band. See [Pg.15, 1].
        :vartype systemNoiseFigure: float

        :ivar radarLoss: (decibels) These include a variety of losses primarily over the microwave signal path, but doesn't include the atmosphere. Included are a power loss from transmitter power amplifier
                           output to the antenna port, and a two-way loss through the radome. Typical numbers might be 0.5 dB to 2 dB from TX amplifier to the
                           antenna port, and perhaps an additional 0.5 dB to 1.5 dB two-way through the radome. See [Pg.15, 1].
        :vartype radarLoss: float

        :ivar atmosLoss: 2-way atmospheric loss of electromagnetic energy (see [Pg.16, 1]).
        :vartype atmosLoss: float     

        :ivar polType: SAR polarization type
        :vartype polType: :class:`instrupy.synthetic_aperture_radar_model.PolTypeSAR`

        :ivar dualPolPulseConfig: In case of dual-pol, this parameter indicates the pulse configuration.
        :vartype dualPolPulseConfig: :class:`DualPolPulseConfig` or None

        :ivar dualPolPulseSep: In case of SMAP dual-pol configuration, this parameter indicates the pulse separation in seconds.
        :vartype dualPolPulseSep: float or None

        :ivar scanTechnique: Scanning technique.
        :vartype scanTechnique: :class:`instrupy.synthetic_aperture_radar_model.ScanTech`

        :ivar swathType: Swath configuration.
        :vartype swathType: :class:`instrupy.synthetic_aperture_radar_model.SwathTypeSAR`       

        :ivar fixedSwathSize: In case of fixed swath configuration this parameter indicates the size of the fixed swath in kilometers.
        :vartype fixedSwathSize: float or None

        :ivar numSubSwaths: Number of sub-swaths (scans) in case of "ScanSAR" operation.
        :vartype numSubSwaths: int or None

        :ivar _id: Unique instrument identifier.
        :vartype _id: str or int

        .. note:: The actual pulse-repetition frequency during the calculation of the observation metrics is taken as the highest PRF within allowed range of PRFs. 
                  The highest PRF is chosen since it maximizes NESZ.
          
    """
    L_r = float(1.2)
    L_a = float(1.2)  
    a_wa = float(1.2)  
    a_wr = float(1.2)    

    def __init__(self, name=None, mass=None, volume=None, power=None,  orientation=None, 
            fieldOfViewGeometry = None, sceneFieldOfViewGeometry = None, maneuver = None, pointingOption=None, 
            dataRate=None, bitsPerPixel = None, pulseWidth = None, antennaHeight= None, 
            antennaWidth = None, antennaApertureEfficiency = None, operatingFrequency = None, 
            peakTransmitPower = None, chirpBandwidth = None, minimumPRF = None, maximumPRF = None, 
            radarLoss = None, atmosLoss=None, sceneNoiseTemp = None, systemNoiseFigure = None,
            polType = None, dualPolPulseConfig = None, dualPolPulseSep = None, swathType = None, 
            scanTechnique = None, fixedSwathSize = None, numSubSwaths = None,  _id=None):
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
        self.power = float(power) if power is not None else None
        self.orientation = copy.deepcopy(orientation) if orientation is not None and isinstance(orientation, Orientation) else None
        self.fieldOfView = ViewGeometry(orien=self.orientation, sph_geom=fieldOfViewGeometry) if self.orientation is not None and fieldOfViewGeometry is not None and isinstance(fieldOfViewGeometry, SphericalGeometry) else None
        self.sceneFieldOfView = ViewGeometry(orien=self.orientation, sph_geom=sceneFieldOfViewGeometry) if self.orientation is not None and sceneFieldOfViewGeometry is not None and isinstance(sceneFieldOfViewGeometry, SphericalGeometry) else None
        self.maneuver = copy.deepcopy(maneuver) if maneuver is not None and isinstance(maneuver, Maneuver) else None
        self.fieldOfRegard = self.maneuver.calc_field_of_regard(self.sceneFieldOfView.sph_geom) if (self.maneuver is not None and self.sceneFieldOfView is not None) else None
        self.pointingOption = None
        if isinstance(pointingOption, list):
            if all(isinstance(x, Orientation) for x in pointingOption):
                self.pointingOption = pointingOption 
        elif isinstance(pointingOption, Orientation):
            self.pointingOption = [pointingOption] # make into single-element list
        self.dataRate = float(dataRate) if dataRate is not None else None          
        self.bitsPerPixel = int(bitsPerPixel) if bitsPerPixel is not None else None 
        self.pulseWidth = float(pulseWidth) if pulseWidth is not None else None
        self.antennaHeight = float(antennaHeight) if antennaHeight is not None else None
        self.antennaWidth = float(antennaWidth) if antennaWidth is not None else None
        self.antennaApertureEfficiency = float(antennaApertureEfficiency) if antennaApertureEfficiency is not None else None
        self.operatingFrequency = float(operatingFrequency) if operatingFrequency is not None else None
        self.peakTransmitPower = float(peakTransmitPower) if peakTransmitPower is not None else None
        self.chirpBandwidth = float(chirpBandwidth) if chirpBandwidth is not None else None        
        self.minimumPRF = float(minimumPRF) if minimumPRF is not None else None   
        self.maximumPRF = float(maximumPRF) if maximumPRF is not None else None 
        self.sceneNoiseTemp = float(sceneNoiseTemp) if sceneNoiseTemp is not None else float(290) # 290 K is default  
        self.radarLoss = float(radarLoss) if radarLoss is not None else None  
        self.atmosLoss = float(atmosLoss) if atmosLoss is not None else None         
        self.systemNoiseFigure = float(systemNoiseFigure) if systemNoiseFigure is not None else None 
        self.polType = PolTypeSAR.get(polType) if polType is not None else None  
        self.dualPolPulseConfig = DualPolPulseConfig.get(dualPolPulseConfig) if dualPolPulseConfig is not None else None  
        self.dualPolPulseSep = float(dualPolPulseSep) if dualPolPulseSep is not None else None  
        self.scanTechnique = ScanTech.get(scanTechnique) if scanTechnique is not None else None
        self.swathType = SwathTypeSAR.get(swathType) if swathType is not None else None  
        self.fixedSwathSize = float(fixedSwathSize) if fixedSwathSize is not None else None  
        self.numSubSwaths = int(numSubSwaths) if numSubSwaths is not None else None

        super(SyntheticApertureRadarModel,self).__init__(_id, "Synthetic Aperture Radar")
        
    @staticmethod
    def from_dict(d):
        """ Parses an SAR instrument from a normalized JSON dictionary.

        The following default values are assigned to the object instance parameters in case of 
        :class:`None` values or missing key/value pairs in the input dictionary.

        .. csv-table:: Default values
            :header: Parameter, Default Value
            :widths: 10,40

            scanTech, ScanTech.STRIPMAP
            orientation, Orientation.Convention.SIDE_LOOK at 25 deg
            sceneFieldOfViewGeometry, (Instrument) fieldOfViewGeometry
            polType, PolTypeSAR.SINGLE (single transmit and single receive)
            pulseSeparation, 50% of pulse length 
            swathType, SwathTypeSAR.FULL  
            fixedSwathSize, 10 km
            numSubSwaths, 1
            atmosLoss, 2 dB
            _id, random string
        
        :param d: Normalized JSON dictionary with the corresponding model specifications. 
        :paramtype d: dict

        :returns: SyntheticApertureRadarModel object initialized with the input specifications.
        :rtype: :class:`instrupy.SyntheticApertureRadarModel`

        """
        # Only side-looking orientation of instrument supported for synthetic aperture radar 
        orien_dict = d.get("orientation", {"convention": "SIDE_LOOK", "sideLookAngle":25})
        orien_conv= orien_dict.get("convention",None)
        if(orien_conv is None):
            raise Exception("Please specify valid Orientation convention.")
        if(Orientation.Convention.get(orien_conv) != Orientation.Convention.SIDE_LOOK):
            raise Exception("Only side-looking orientation of instrument supported for the synthetic aperture radar imaging.")
            
        # check supplied PRF range
        _PRFmin = d.get("minimumPRF", None)
        _PRFmax = d.get("maximumPRF", None)
        if(_PRFmin > _PRFmax):
            raise Exception("PRF minimum must be numerically less than or equal to PRF maximum.")

        # initialize the polarization configuration
        pol = d.get("polarization", None)
        dualPolPulseConfig = None
        dualPolPulseSep = None
        if(pol):
            polType = PolTypeSAR.get(pol.get("@type",None))
            if(polType == PolTypeSAR.DUAL):
                if(pol.get("pulseConfig", None)):
                    if(pol["pulseConfig"].get("@type", None) == DualPolPulseConfig.SMAP):
                        dualPolPulseConfig = DualPolPulseConfig.SMAP
                        dualPolPulseSep = pol["pulseConfig"].get("pulseSeparation") if pol["pulseConfig"].get("pulseSeparation", None) is not None else 0.5*d.get("pulseWidth")
                    elif(pol["pulseConfig"].get("@type", None) == DualPolPulseConfig.AIRSAR):
                        dualPolPulseConfig = DualPolPulseConfig.AIRSAR
                        dualPolPulseSep = None
                    else:
                        raise RuntimeError("Unknown/ invalid pulse configuration in dual-polarization specification.")
            elif(polType == PolTypeSAR.SINGLE or polType == PolTypeSAR.COMPACT):
                pass
            else:
                raise RuntimeError("Unknown/ invalid polarization specification.")
        else: # assign default polarization
            polType = PolTypeSAR.SINGLE            

        _scan = ScanTech.get(d.get("scanTechnique", None))
        if _scan is None:
            _scan = ScanTech.STRIPMAP # default scan technique
        # Only stripmap and scansar techniques are supported
        if(_scan == ScanTech.STRIPMAP) or (_scan == ScanTech.SCANSAR):
            
            if(_scan == ScanTech.SCANSAR):
                numSubSwaths = d.get("numSubSwaths",None)                
                if numSubSwaths is None:
                    numSubSwaths = 1 # default
            elif(_scan == ScanTech.STRIPMAP):
                numSubSwaths = d.get("numSubSwaths",None)               
                if numSubSwaths is not None:
                    warnings.warn("numSubSwaths parameter is not considered for Stripmap operation. It shall be ignored.")
                numSubSwaths = 1 # stripmap is equivalent to ScanSAR operation with one subswath

            # calculate instrument FOV based on antenna dimensions
            D_az_m = d.get("antennaHeight", None)
            D_elv_m = d.get("antennaWidth", None)
            opWavelength =  Constants.speedOfLight/ d.get("operatingFrequency", None)
            # calculate antenna beamwidth and hence instrument FOV [eqn41, 1].
            along_track_fov_deg = np.rad2deg(opWavelength/ D_az_m)
            cross_track_fov_deg = numSubSwaths*np.rad2deg(opWavelength/ D_elv_m) # number of subswaths x antenna cross-track beamwidth
            instru_fov_geom_dict = { "shape": "RECTANGULAR", "angleHeight":along_track_fov_deg, "angleWidth": cross_track_fov_deg } 
            
            scene_fov_geom_dict = d.get("sceneFieldOfViewGeometry", instru_fov_geom_dict)  # default sceneFOV geometry is the instrument FOV geometry
            # initialize the swath configuration
            fixedSwathSize = None
            swathConfig = d.get("swathConfig", None)
            
            if(swathConfig):
                swathType = SwathTypeSAR.get(swathConfig.get("@type", None))
                
                if(swathType == SwathTypeSAR.FIXED):
                    
                    fixedSwathSize = swathConfig.get("fixedSwathSize") if swathConfig.get("fixedSwathSize", None) is not None else 10 # default to 10km
                    if(_scan == ScanTech.SCANSAR):
                        fixedSwathSize = None
                        swathConfig = SwathTypeSAR.FULL  
                        warnings.warn("ScanSAR operation supports only FULL swath configuration. Specified FIXED swath configuration is ignored.")
                       
            else: # assign default
                swathType = SwathTypeSAR.FULL       
        else:
            raise RuntimeError("Unknown SAR scan technique specified.")
        
        # parse the pointing options as a list of Orientation objects.
        pnt_opt_dict = d.get("pointingOption", None)
        _pointing_option = None
        if pnt_opt_dict:
            # translate to a list of Orientation objects
            if isinstance(pnt_opt_dict, list):
                _pointing_option = [Orientation.from_dict(x) for x in pnt_opt_dict]
            else:
                _pointing_option = [Orientation.from_dict(pnt_opt_dict)]

        return SyntheticApertureRadarModel(
                        name = d.get("name", None),
                        mass = d.get("mass", None),
                        volume = d.get("volume", None),
                        power = d.get("power", None),
                        orientation = Orientation.from_dict(orien_dict),
                        fieldOfViewGeometry =  SphericalGeometry.from_json(instru_fov_geom_dict),
                        sceneFieldOfViewGeometry = SphericalGeometry.from_json(scene_fov_geom_dict),
                        maneuver =  Maneuver.from_json(d.get("maneuver", None)),
                        pointingOption = _pointing_option,
                        dataRate = d.get("dataRate", None),
                        bitsPerPixel = d.get("bitsPerPixel", None),
                        pulseWidth = d.get("pulseWidth", None),
                        antennaHeight = d.get("antennaHeight", None),
                        antennaWidth = d.get("antennaWidth", None),
                        antennaApertureEfficiency = d.get("antennaApertureEfficiency", None),
                        operatingFrequency = d.get("operatingFrequency", None),
                        peakTransmitPower = d.get("peakTransmitPower", None),
                        chirpBandwidth = d.get("chirpBandwidth", None),
                        minimumPRF = d.get("minimumPRF", None),
                        maximumPRF = d.get("maximumPRF", None),
                        radarLoss = d.get("radarLoss", None),
                        atmosLoss=d.get("atmosLoss", 2), # 2dB default value
                        sceneNoiseTemp = d.get("sceneNoiseTemp", None),
                        systemNoiseFigure = d.get("systemNoiseFigure", None),
                        polType=polType,
                        dualPolPulseConfig = dualPolPulseConfig,
                        dualPolPulseSep=dualPolPulseSep,
                        swathType=swathType,
                        scanTechnique = _scan,
                        fixedSwathSize=fixedSwathSize,
                        numSubSwaths=numSubSwaths,
                        _id = d.get("@id", uuid.uuid4())
                        )

    def to_dict(self):
        """ Translate the SyntheticApertureRadarModel object to a Python dictionary such that it can be uniquely reconstructed back from the dictionary.

        :returns: SyntheticApertureRadarModel specifications as python dictionary.
        :rtype: dict

        """
        fieldOfViewGeometry_dict = self.fieldOfView.sph_geom.to_dict() if self.fieldOfView is not None and isinstance(self.fieldOfView, ViewGeometry) else None
        sceneFieldOfViewGeometry_dict = self.sceneFieldOfView.sph_geom.to_dict() if self.sceneFieldOfView is not None and isinstance(self.sceneFieldOfView, ViewGeometry) else None
        orientation_dict = self.orientation.to_dict() if self.orientation is not None and isinstance(self.orientation, Orientation) else None
        maneuver_dict = self.maneuver.to_dict() if self.maneuver is not None and isinstance(self.maneuver, Maneuver) else None
        pointing_opt_dict = [Orientation.to_dict(x) for x in self.pointingOption] if self.pointingOption is not None else None
        return dict({
                "@type": "Synthetic Aperture Radar",
                "name":self.name,
                "mass":self.mass,
                "volume":self.volume,
                "power":self.power,
                "orientation":orientation_dict,
                "fieldOfViewGeometry":fieldOfViewGeometry_dict,
                "sceneFieldOfViewGeometry": sceneFieldOfViewGeometry_dict,                
                "maneuver":maneuver_dict,
                "pointingOption": pointing_opt_dict,
                "dataRate":self.dataRate,
                "bitsPerPixel": self.bitsPerPixel,
                "pulseWidth": self.pulseWidth,
                "antennaHeight": self.antennaHeight,
                "antennaWidth": self.antennaWidth,
                "antennaApertureEfficiency": self.antennaApertureEfficiency,
                "operatingFrequency": self.operatingFrequency,
                "peakTransmitPower": self.peakTransmitPower,
                "chirpBandwidth": self.chirpBandwidth,
                "minimumPRF": self.minimumPRF,
                "maximumPRF": self.maximumPRF,
                "radarLoss": self.radarLoss,
                "atmosLoss": self.atmosLoss,
                "sceneNoiseTemp": self.sceneNoiseTemp,
                "systemNoiseFigure": self.systemNoiseFigure,
                "polType": self.polType.value,
                "pulseConfig": self.dualPolPulseConfig.value if self.dualPolPulseConfig is not None else None,
                "pulseSeparation": self.dualPolPulseSep,
                "swathType": self.swathType.value,
                "scanTechnique":self.scanTechnique.value if self.dualPolPulseConfig is not None else None,
                "fixedSwathSize": self.fixedSwathSize,
                "numSubSwaths": self.numSubSwaths,
                "@id": self._id
                })

    def __repr__(self):
        return "SyntheticApertureRadarModel.from_dict({})".format(self.to_dict())

    @staticmethod
    def get_azimuthal_resolution(sc_speed, sc_gnd_speed, D_az):
        """ Calculate azimuthal resolution taking into consideration difference in spacecraft and footprint velocities.
            See eqn (5.3.6.3) in [2].

        :param sc_speed: [distance/time] Spacecraft speed.
        :paramtype sc_speed: float

        :param sc_gnd_speed: [distance/time] Spacecraft ground speed. Units must be consistent with the the spacecraft speed input.
        :paramtype sc_gnd_speed: float

        :param D_az: Length of antenna in meters along the azimuthal direction.
        :paramtype D_az: float

        :returns: Azimuthal resolution in meters.
        :rtype: float

        """
        return (D_az/2.0)*(sc_gnd_speed/ sc_speed)


    def calc_data_metrics(self, sc_orbit_state=None, target_coords=None, 
                          alt_km=None, sc_speed_kmps=None, sc_gnd_speed_kmps=None, inc_angle_deg=None, 
                          instru_look_angle_from_target_inc_angle=False):
        """ A wrapper function for calling the appropriate (depending on the input arguments) implementation of the data metrics calculator.
            Refer to the functions ``calc_data_metrics_impl1`` and ``calc_data_metrics_impl2`` for description of the function parameters.
            This function is invoked by the function ``Instrument.calc_data_metrics(.)`` class in the ``base`` module.

        :returns: Calculated observation data metrics.
        :rtype: dict

        """        
        if(sc_orbit_state is not None and target_coords is not None):
        
            obsv_metrics = SyntheticApertureRadarModel.calc_data_metrics_impl2(self, sc_orbit_state, target_coords, instru_look_angle_from_target_inc_angle)

        elif(alt_km is not None and sc_speed_kmps is not None and sc_gnd_speed_kmps is not None and inc_angle_deg is not None):
        
            obsv_metrics = SyntheticApertureRadarModel.calc_data_metrics_impl1(self, alt_km, sc_speed_kmps, sc_gnd_speed_kmps, inc_angle_deg, instru_look_angle_from_target_inc_angle)
        
        else:
            raise RuntimeError("Required set of arguments not present in the SyntheticApertureRadarModel.calc_data_metrics(.) function.")

        return obsv_metrics


    def calc_data_metrics_impl2(self, sc_orbit_state, target_coords, instru_look_angle_from_target_inc_angle=False):
        """ Calculate observation data metrics.

        :param sc_orbit_state: Spacecraft state at the time of observation.

        Dictionary keys are: 
        
        * :code:`time [JDUT1]` (:class:`float`), Time in Julian Day UT1. Corresponds to the time of observation. 
        * :code:`x [km]` (:class:`float`), :code:`y [km]` (:class:`float`), :code:`z [km]` (:class:`float`), Cartesian spatial coordinates of satellite in EARTH_CENTERED_INERTIAL frame at the time of observation.
        * :code:`vx [km/s]` (:class:`float`), :code:`vy [km/s]` (:class:`float`), :code:`vz [km/s]` (:class:`float`), Velocity of spacecraft in EARTH_CENTERED_INERTIAL frame at the time of observation.
        
        :paramtype sc_orbit_state: dict
        
        :param target_coords: Location of the observation. Also sometimes the Point-Of-Interest (POI).

        Dictionary keys are: 
        
        * :code:`lat [deg]` (:class:`float`), :code:`lon [deg]` (:class:`float`) indicating the corresponding ground-point accessed (latitude, longitude) in degrees.

        :paramtype target_coords: dict

        :param instru_look_angle_from_target_inc_angle: Flag (True or False) to indicate if the look angle to the middle of the swath is to be considered:  (1) using the nominal look-angle 
                                                        (specified in the ``orientation`` attribute of the instrument), 
                                                        OR
                                                        (2) the incidence angle at the target. 
                                                        Default is False.
        
        :paramtype instru_look_angle_from_target_inc_angle: bool

        :returns: Calculated observation data metrics.

        Dictionary keys are: 
    
        * :code:`NESZ [dB]` (:class:`float`)  The backscatter coefficient :math:`\\sigma_0` of a target for which the signal power level in final image is equal to the noise power level (units: decibels). **Numerically lesser is better.**
        * :code:`ground pixel along-track resolution [m]` (:class:`float`) Along-track resolution (meters) of an ground-pixel centered about observation point.
        * :code:`ground pixel cross-track resolution [m]` (:class:`float`) Cross-track resolution (meters) of an ground-pixel centered about observation point.
        * :code:`swath-width [m]` (:class:`float`) Swath-width (meters) of the strip of which the imaged pixel is part off.
        * :code:`incidence angle [deg]` (:class:`float`) Observation incidence angle (degrees) at the ground-pixel.
        * :code:`PRF [Hz]` (:class:`float`)  Highest Pulse Repetition Frequency (Hz) (within the specified PRF range) at which the observation is possible.

        :rtype: dict
                        
        """       
        # Observation time in Julian Day UT1
        tObs_JDUT1 = sc_orbit_state["time [JDUT1]"]

        # Calculate Target position in CARTESIAN_EARTH_CENTERED_INERTIAL frame
        target_pos_km = GeoUtilityFunctions.geo2eci([target_coords["lat [deg]"], target_coords["lon [deg]"], 0.0], tObs_JDUT1)
        # Spacecraft position in CARTESIAN_EARTH_CENTERED_INERTIAL frame
        sc_pos_km = np.array([sc_orbit_state["x [km]"], sc_orbit_state["y [km]"], sc_orbit_state["z [km]"]])  
        sc_vel_kmps = np.array([sc_orbit_state["vx [km/s]"], sc_orbit_state["vy [km/s]"], sc_orbit_state["vz [km/s]"]]) 
        sc_speed_kmps = np.linalg.norm(sc_vel_kmps)
        
        # Calculate range vector between spacecraft and POI (Target)
        range_vector_km = target_pos_km - sc_pos_km

        alt_km = np.linalg.norm(sc_pos_km) - Constants.radiusOfEarthInKM

        look_angle = np.arccos(np.dot(MathUtilityFunctions.normalize(range_vector_km), -1*MathUtilityFunctions.normalize(sc_pos_km)))
        incidence_angle_rad = np.arcsin(np.sin(look_angle)*(Constants.radiusOfEarthInKM + alt_km)/Constants.radiusOfEarthInKM)       


        sc_gnd_speed_kmps = 1e-3 * GeoUtilityFunctions.compute_satellite_footprint_speed(sc_pos_km*1e3, sc_vel_kmps*1e3) # This is approximation, since the image footprint velocity is not necessarily equal to the
                                        # satellite footprint speed. However it is reasonable approximation in case of low-altitudes and small look angles. TODO: Improve the model.
    
        #print("sc_speed_kmps", sc_speed_kmps)
        #print("sc_gnd_speed_kmps", sc_gnd_speed_kmps)
        obsv_metrics = SyntheticApertureRadarModel.calc_data_metrics_impl1(self, alt_km, sc_speed_kmps, sc_gnd_speed_kmps, np.rad2deg(incidence_angle_rad), instru_look_angle_from_target_inc_angle)

        return obsv_metrics

    def calc_data_metrics_impl1(self, alt_km, sc_speed_kmps, sc_gnd_speed_kmps, inc_angle_deg, instru_look_angle_from_target_inc_angle=False):
        """ Calculate the observation metrics.

        :param alt_km: Spacecraft altitude in kilometers.
        :paramtype alt_km: float

        :param sc_speed_kmps: Spacecraft speed in kilometers per second.
        :paramtype sc_speed_kmps: float

        :param sc_gnd_speed_kmps: Spacecraft ground-speed in kilometers per second.
        :paramtype sc_gnd_speed_kmps: float

        :param inc_angle_deg: Incidence angle in degrees at the target position (w.r.t the spacecraft position).
        :paramtype inc_angle_deg: float

        :param instru_look_angle_from_target_inc_angle: Flag (True or False) to indicate if the look angle to the middle of the swath is to be considered:  (1) using the nominal look-angle 
                                                        (specified in the ``orientation`` attribute of the instrument), 
                                                        OR
                                                        (2) the incidence angle at the target. 
                                                        Default is False.
        :paramtype instru_look_angle_from_target_inc_angle: bool

        :returns: Calculated observation data metrics. Refer to the return value of the function ``calc_data_metrics_impl2(.)``.

        :rtype: dict

        """
        inc_angle = np.deg2rad(inc_angle_deg)
        look_angle = np.arcsin(np.sin(inc_angle)*Constants.radiusOfEarthInKM/(Constants.radiusOfEarthInKM + alt_km)) 
        #print("inc_angle", inc_angle*180/np.pi)
        #print("look_angle", look_angle*180/np.pi)
        range_km = Constants.radiusOfEarthInKM * (np.sin(inc_angle - look_angle)/ np.sin(look_angle))

        if(instru_look_angle_from_target_inc_angle):
            instru_look_angle_rad = look_angle # instrument look angle from the target incidence angle
        else:
            instru_look_angle_rad = np.abs(np.deg2rad(self.orientation.euler_angle2)) # instrument look angle from the instrument orientation

        # Copying values into variables of more code-friendly variables
        c = Constants.speedOfLight
        k = Constants.Boltzmann
        tau_p = self.pulseWidth        
        B_T = self.chirpBandwidth
        P_T = self.peakTransmitPower
        D_az = self.antennaHeight
        D_elv = self.antennaWidth
        fc = self.operatingFrequency
        eta_ap = self.antennaApertureEfficiency               
        PRFmin_Hz = self.minimumPRF    
        PRFmax_Hz = self.maximumPRF   
        L_radar = 10.0**(self.radarLoss/10.0) # convert to linear units
        F_N = 10.0**(self.systemNoiseFigure/10.0) # convert to linear units
        L_atmos = 10.0**(self.atmosLoss/10.0) # convert to linear units
        L_r = SyntheticApertureRadarModel.L_r
        L_a = SyntheticApertureRadarModel.L_a
        a_wr = SyntheticApertureRadarModel.a_wr
        a_wa = SyntheticApertureRadarModel.a_wa
        T = self.sceneNoiseTemp       

        #print(tau_p)
        #print(fc)
        #print(self.polType)
        # Note that the nominal look angle is considered to evaluate the operable PRF.
        (f_P_master, W_gr_obs) = SyntheticApertureRadarModel.prf_constraint_eval(PRFmin_Hz, PRFmax_Hz, sc_speed_kmps, sc_gnd_speed_kmps, alt_km, 
                                                                             instru_look_angle_rad, tau_p, D_az, D_elv, fc,
                                                                             self.polType, self.dualPolPulseConfig, self.dualPolPulseSep, 
                                                                             self.swathType, self.fixedSwathSize, self.numSubSwaths)
        rho_a = None
        rho_y = None
        sigma_N_dB = None
        theta_i = None

        #print("f_P_master", f_P_master)

        if (f_P_master is not None): # Observation is (perhaps (since determined at nominal instrument look-angle, some parts of the swath may not be observable)) possible at PRF = f_P        
            
            R = range_km*1e3

            lamb = c/fc

            # Note that this may not be the same as incidence angle to middle of swath.
            theta_i = inc_angle
                            
            psi_g = np.pi/2.0 - theta_i # grazing angle   
            
            # [1] equation 17, find P_avg (average transmit power)
            T_eff = tau_p # approximate effective pulse duration to be the actual pulse duration, as in case of matched filter processing
            # In case of AIRSAR dual-pol configuration, the channel PRF is half of the master PRF.
            if (self.polType == PolTypeSAR.DUAL and self.dualPolPulseConfig == DualPolPulseConfig.AIRSAR):
                f_P_ch = 0.5*f_P_master
            else:
                f_P_ch = f_P_master
            d = T_eff * f_P_ch # [1] equation 17
            P_avg = d*P_T
            
            # [1] equation 8, find G_A
            A_A = D_elv*D_az
            G_A = 4.0*np.pi*eta_ap*A_A/lamb**2  
            
            # [1] equation 37 we can get the sigma_N. Note that the spacecraft speed is used in the equation and not the ground-speed, see [2] for explanation.              
            sigma_N = (265.0*np.pi**3*k*T / c)*(R**3 * (sc_speed_kmps*1e3) * np.cos(psi_g))*(B_T*F_N*L_radar*L_atmos/ (P_avg*G_A**2*lamb**3))*(L_r*L_a/(a_wr*a_wa))
            sigma_N_dB = 10.0*np.log10(sigma_N)

            # [1] equations 36, 23 we can get rho_y
            rho_y = a_wr*c/(2*B_T*np.cos(psi_g))

            # [1] equation 69 we get minimum possible azimuth resolution (for strip mapping)
            rho_a = SyntheticApertureRadarModel.get_azimuthal_resolution(sc_speed_kmps, sc_gnd_speed_kmps, D_az)
            # modify in case of scansar (multiple sub-swaths => trading off azimuthal resolution)
            rho_a = rho_a*self.numSubSwaths
             
        obsv_metrics = {}
        obsv_metrics["ground pixel along-track resolution [m]"] = round(rho_a, 2) if rho_a is not None else np.nan
        obsv_metrics["ground pixel cross-track resolution [m]"] = round(rho_y, 2) if rho_y is not None else np.nan
        obsv_metrics["NESZ [dB]"] = round(sigma_N_dB, 2) if sigma_N_dB is not None else np.nan
        obsv_metrics["incidence angle [deg]"] = round(np.rad2deg(theta_i), 2) if theta_i is not None else np.nan
        obsv_metrics["swath-width [km]"] = round(W_gr_obs/1e3, 1) if W_gr_obs is not None else np.nan        
        obsv_metrics["PRF [Hz]"] = f_P_master

        return obsv_metrics

    @staticmethod
    def prf_constraint_eval(f_P_min, f_P_max, sc_speed_kmps, sc_gnd_speed_kmps, alt_km, 
                            instru_look_angle_rad, tau_p_ch, D_az, D_elv, fc, 
                            pol_type=PolTypeSAR.SINGLE, dual_pol_conf=None, dual_pol_ps=None, 
                            swath_type=SwathTypeSAR.FULL, fixed_swath_size_km=10,
                            num_sub_swath=1):
        """ Function to find the highest possible pulse repetition frequency within the user supplied range of PRF, which 
            shall allow observation of target. Not all PRFs are valid and a valid PRF has to be chosen so that it meets all
            the below conditions:

            1. The length of the echo from illuminated/ desired swath is less than inter-pulse period.
            2. The PRFs should be high enough to allow for unambiguous detection of Doppler shifts.
            3. A transmit pulse does not overlap with the desired echo.
            4. The Nadir echoes from other transmit pulses do not overlap with the desired echo.

        [2] is the primary reference for this formulation, although some errors have been found (and corrected for the current
        implementation). [3] contains the corrections. The referenced formulation is further modified to incorporate the PRF constraints
        involving observations of multiple polarizations and fixed-swath (desired echo vs complete echo).

        Of all the available valid PRFs, the highest PRF is chosen since it improves the NESZ observation data-metric.
        The near-range and far-range calculations are based on the input instrument look-angle. 

        :param f_P_min: Minimum PRF in [Hz].
        :paramtype f_P_min: float

        :param f_P_max: Maximum PRF in [Hz].
        :paramtype f_P_max: float

        :param sc_speed_kmps: Satellite speed in [km/s].
        :paramtype sc_speed_kmps: float

        :param sc_gnd_speed_kmps: Satellite ground speed in [km/s].
        :paramtype sc_gnd_speed_kmps: float

        :param alt_km: Altitude of satellite in [km].
        :paramtype alt_km: float

        :param instru_look_angle_rad: Instrument look angle (to the middle of the swath) in [radians].
        :paramtype instru_look_angle_rad: float

        :param tau_p_ch: Pulse width in [s] per polarization (channel).
        :paramtype tau_p_ch: float

        :param D_az: Antenna dimension along cross-range direction in [m].
        :paramtype D_az: float

        :param D_elv: Antenna dimension along range direction in [m].
        :paramtype D_elv: float

        :param fc: Carrier center frequency in [Hz].
        :paramtype fc: float

        :param pol_type: SAR polarization type (default PolTypeSAR.SINGLE).
        :paramtype pol_type: :class:`PolTypeSAR`

        :cvar dual_pol_conf: In case of SMAP dual-pol configuration, this parameter indicates the pulse configuration.
        :vartype dual_pol_conf: :class:`DualPolPulseConfig` or None

        :param dual_pol_ps: In case of SMAP dual-pol configuration, this parameter indicates the pulse seperation in seconds.
        :paramtype dual_pol_ps: float or None

        :param swath_type: Desired swath type. For ScanSAR only "FULL" is accepted. (default = SwathTypeSAR.FULL)
        :paramtype swath_type: :class:`SwathTypeSAR` 

        :param fixed_swath_size_km: In case of fixed swath configuration this parameter indicates the size of the fixed swath. Not applicable for ScanSAR.
        :paramtype fixed_swath_size_km: float

        :param num_sub_swath: Number of subswaths (default = 1). In case of ScanTech.SCANSAR the number of sub-swaths maybe more than 1.
        :paramtype num_sub_swath: int

        :returns: Tuple with the highest possible master PRF which can be used for observation and observed swath-width.
        :rtype: tuple (int, float)

        """
        h = alt_km * 1e3
        Re = Constants.radiusOfEarthInKM * 1e3         
        c = Constants.speedOfLight
        Rs = Re + h # [2]  part of equation 5.1.2.3

        lamb = c/fc        

        theta_elv = lamb/ D_elv # null-to-null (full-beamwidth) illuminating the (sub)swath
        
        if(swath_type == SwathTypeSAR.FULL and num_sub_swath>1): # condition is true in case of ScanSAR operation with multiple sub-swaths
            """
                Consider the PRF evaluation for the farthest (off-nadir) sub-swath. In practice the different sub-swaths
                may have different operating PRFs. The PRF is expected to be most constrained for the farthest sub-swath.
            """
            y = instru_look_angle_rad + 0.5*num_sub_swath * theta_elv 
            gamma_m = y - 0.5 * theta_elv # modify the look angle to that of the farthest subswath
        
        elif(swath_type == SwathTypeSAR.FIXED and num_sub_swath>1):
            raise NotImplementedError

        elif(num_sub_swath == 1): # condition is reached for Stripmap operation or ScanSAR operation with 1 subswath
            gamma_m = instru_look_angle_rad
            
        else:
            raise RuntimeError("Unknown condition reached.")

        # NOTE: While calculating swath width (FULL or FIXED) the instrument look angle is used!!!! Not the Target look angle (which could be equal to or not equal to the instrument look angle).
        gamma_n_illum = gamma_m - 0.5*theta_elv
        gamma_f_illum = gamma_m + 0.5*theta_elv
        theta_in_illum = np.arcsin(np.sin(gamma_n_illum)*Rs/Re)
        theta_horizon = np.arcsin(Re/Rs)
        try:
            theta_if_illum = np.arcsin(np.sin(gamma_f_illum)*Rs/Re)
        except:
            # beyond horizon, hence set to horizon angle
            theta_if_illum = theta_horizon

        alpha_n_illum = theta_in_illum - gamma_n_illum
        alpha_f_illum = theta_if_illum - gamma_f_illum
        alpha_s_illum = alpha_f_illum - alpha_n_illum
        W_gr_illum = Re*alpha_s_illum  # illuminated (sub)swath

        Rn_illum = np.sqrt(Re**2 + Rs**2 - 2*Re*Rs*np.cos(alpha_n_illum)) # [2] equation 5.1.3.9 slant-range to near edge of swath
        Rf_illum = np.sqrt(Re**2 + Rs**2 - 2*Re*Rs*np.cos(alpha_f_illum)) # [2] equation 5.1.3.10 slant-range to far edge of swath
        
        tau_near_illum = 2*Rn_illum/c # [2] equation 5.1.3.11
        tau_far_illum = 2*Rf_illum/c # [2] equation 5.1.3.12

        if(swath_type == SwathTypeSAR.FULL):
            W_gr_obs = W_gr_illum  # desired (observed) swath

        elif(swath_type == SwathTypeSAR.FIXED): # fixed swath
            
            if(fixed_swath_size_km*1e3 < W_gr_illum):
                W_gr_obs = fixed_swath_size_km*1e3
            else:
                W_gr_obs = W_gr_illum

            alpha_s_obs = W_gr_obs/Re
            theta_im = np.arcsin(np.sin(gamma_m)*Rs/Re)
            alpha_m = theta_im - gamma_m  # [2] equation 5.1.3.5
            alpha_n_obs = alpha_m - alpha_s_obs/2.0 # [2] equation 5.1.3.7
            alpha_f_obs = alpha_m + alpha_s_obs/2.0 # [2] equation 5.1.3.8

            Rn_obs = np.sqrt(Re**2 + Rs**2 - 2*Re*Rs*np.cos(alpha_n_obs)) # [2] equation 5.1.3.9 slant-range to near edge of swath
            Rf_obs = np.sqrt(Re**2 + Rs**2 - 2*Re*Rs*np.cos(alpha_f_obs)) # [2] equation 5.1.3.10 slant-range to far edge of swath
            
            tau_near_obs = 2*Rn_obs/c # [2] equation 5.1.3.11
            tau_far_obs = 2*Rf_obs/c # [2] equation 5.1.3.12

        else:
            raise RuntimeError("Unknown swath configuration type.")

        PRFmin_f = 1 # Factor by which the minimum PRF constraint must be increased. Becomes significant in case of dual-pol.
        if(pol_type == PolTypeSAR.SINGLE or pol_type == PolTypeSAR.COMPACT):
            tau_p = tau_p_ch
        elif(pol_type == PolTypeSAR.DUAL):
            if(dual_pol_conf == DualPolPulseConfig.AIRSAR):
                tau_p = tau_p_ch
                PRFmin_f = 2
            elif(dual_pol_conf == DualPolPulseConfig.SMAP):
                #print(dual_pol_ps)
                tau_p = 2*tau_p_ch + dual_pol_ps            
            else:
                raise RuntimeError("Unknown dual-pol pulse configuration type.")

        else:
            raise RuntimeError("Unknown polarization type.")       
        
        if(swath_type == SwathTypeSAR.FULL):
            PRFmax = 1.0/(2.0*tau_p + tau_far_illum - tau_near_illum) # max allowable PRF [2] equation 5.1.3.13. This condition ensures that only one swath-echo is in a pulse period.
        elif(swath_type == SwathTypeSAR.FIXED):
            # PRF max constraint is loosened by considering a 'desired echo length' which is less than the illuminated echo length.
            PRFmax = 1.0/(2.0*tau_p + tau_far_obs - tau_near_obs) # max allowable PRF [2] equation 5.1.3.13        
        else:
            raise RuntimeError("Unknown swath configuration type.")
        
        # Note that the PRFmin is independent of the number of sub-swaths (PRFmin(stripmap) = PRFmin(scansar)) 
        PRFmin = PRFmin_f*sc_speed_kmps*1e3/SyntheticApertureRadarModel.get_azimuthal_resolution(sc_speed_kmps, sc_gnd_speed_kmps, D_az) # minimum allowable PRF to satisfy Nyquist sampling criteria [2] equation 5.1.2.1 modified to [2] equation (5.4.4.2)
        
        #print("PRFmax: ", PRFmax)
        #print("PRFmin: ", PRFmin)

        f_P = None   # Note that this is the Master PRF
        # Find the highest possible prf within the input prf range which allows for unambiguous echo detection. 
        for _f_P in range(int(f_P_max), int(f_P_min), -1): # step down in steps of 1 Hz
            PRFOK = True

            # perform first and second check of PRF validity
            if((_f_P< PRFmin) or (_f_P>PRFmax)):
                PRFOK = False 
                continue # goto next iteration of _f_P
            
            if(swath_type == SwathTypeSAR.FULL):
                # perform third check of PRF validity, check that target echo is not eclipsed by transmit pulse
                # inequality [2] 5.1.4.1 
                N = int(_f_P*2.0*Rn_illum/c) + 1
                if(not(((N-1)/(tau_near_illum-tau_p) < _f_P ) and (_f_P < N/(tau_far_illum + tau_p)))):
                    PRFOK = False 
                    continue # goto next iteration of _f_P
                            
                # perform fourth check of PRF validity, check that target echo is not eclipsed by nadir echo 
                # from any of the succeeding pulses (from the transmit pulse under consideration to the pulse 
                # just before the echo)
                # Inequality [2] 5.1.5.2 seems wrong. Refer to [3] Appendix Section A for the corrected version (R2 in eqn(38) is a typo, and must be replaced by Rn). 
                tau_nadir = 2.0*h/c
                if(tau_near_illum -tau_nadir - tau_p <= 0): # evaluate the condition independent of f_P
                    PRFOK = False # there is nadir echo overlap with desired echo for the mth pulse
                    break # break used since condition is independent of f_P
                M = int(_f_P*2.0*Rf_illum/c) + 1
                for m in range(1,M):
                    if(((m/(tau_near_illum - tau_p - tau_nadir)) > _f_P) and (_f_P > (m/(tau_far_illum + tau_p - tau_nadir)))):
                        PRFOK = False # there is nadir echo overlap with desired echo for the mth pulse
                        continue # goto next iteration of _f_P 

            elif(swath_type == SwathTypeSAR.FIXED):
                # If fixed swath check only in the desired echo window region

                # perform third check of PRF validity, check that target echo is not eclipsed by transmit pulse
                # inequality [2] 5.1.4.1 
                N = int(_f_P*2.0*Rn_obs/c) + 1
                if(not(((N-1)/(tau_near_obs-tau_p) < _f_P ) and (_f_P < N/(tau_far_obs + tau_p)))):
                    PRFOK = False 
                    continue # goto next iteration of _f_P     

                # perform fourth check of PRF validity, check that target echo is not eclipsed by nadir echo 
                # from any of the succeeding pulses (from the transmit pulse under consideration to the pulse 
                # just before the echo)
                # inequality [2] 5.1.5.2 seems wrong.  Refer to [3] Appendix Section A for the corrected version (R2 in eqn(38) is a type, and must be replaced by Rn). 
                
                tau_nadir = 2.0*h/c
                if(tau_near_obs -tau_nadir - tau_p <= 0): # evaluate the condition independent of f_P
                    PRFOK = False # there is nadir echo overlap with desired echo for the mth pulse
                    break # break used since condition is independent of f_P

                M = int(_f_P*2.0*Rf_obs/c) + 1
                for m in range(1,M):
                    if(((m/(tau_near_obs - tau_p - tau_nadir)) > _f_P) and (_f_P > (m/(tau_far_obs + tau_p - tau_nadir)))):
                        PRFOK = False # there is nadir echo overlap with desired echo 
                        continue # goto next iteration of _f_P 
                
                # To prevent range ambiguity the previous echos from the total illuminated swath should not 
                # overlap with the desired echo (current). This condition can be formulated similar to the Nadir
                # interference condition, albeit with the nadir-echo replaced in terms of main-lobe echo from pulses
                # occuring after the reference pulse until the desired echo window has ended.
                tau_ML_start = tau_near_illum # ML => Main Lobe
                tau_ML_len = tau_far_illum - tau_near_illum # length of main-lobe echo
                W = int(_f_P*2.0*Rf_obs/c) + 1
                for w in range(1,W):                    
                    """
                    # below code causes errors
                    print(1, (w/(tau_near_obs - tau_ML_len - tau_ML_start)) < _f_P)
                    print(2, _f_P < (w/(tau_far_obs + tau_p - tau_ML_start)))
                    if(((w/(tau_near_obs - tau_ML_len - tau_ML_start)) > _f_P) and (_f_P > (w/(tau_far_obs + tau_p - tau_ML_start)))):
                        PRFOK = False # there is main echo overlap with desired echo 
                        print("main echo overlap")
                        continue # goto next iteration of _f_P 
                    """
                    cond1 = (tau_near_obs > tau_ML_len + tau_ML_start + w/_f_P)                    
                    cond2 =  tau_ML_start > tau_far_obs + tau_p - w/_f_P
                    if( (not cond1) and (not cond2)):
                        PRFOK = False # there is main echo overlap with desired echo 
                        #print("main echo overlap")
                        continue # goto next iteration of _f_P 
                

            # If control has reached here and if PRFOK is still True, then this is a valid PRF
            if(PRFOK == True):
                f_P = _f_P
                break

        # compute swath-size
        if(swath_type == SwathTypeSAR.FULL and num_sub_swath>1): # condition is true in case of ScanSAR mode with multiple sub-swaths
            # imaged swath consists of multiple sub-swaths, each sub-swath size is different.
            swath_bw = num_sub_swath * theta_elv # in case of scan-sar, the num_sub_swath > 1
            gamma_n_illum = instru_look_angle_rad - 0.5*swath_bw
            gamma_f_illum = instru_look_angle_rad + 0.5*swath_bw
            theta_in_illum = np.arcsin(np.sin(gamma_n_illum)*Rs/Re)
            theta_horizon = np.arcsin(Re/Rs)
            if(abs(np.sin(gamma_f_illum)*Rs/Re)<1):
                theta_if_illum = np.arcsin(np.sin(gamma_f_illum)*Rs/Re)
            else:
                # beyond horizon, hence set to horizon angle
                theta_if_illum = theta_horizon
            alpha_n_illum = theta_in_illum - gamma_n_illum
            alpha_f_illum = theta_if_illum - gamma_f_illum
            alpha_s_illum = alpha_f_illum - alpha_n_illum
            swath_size = Re*alpha_s_illum  

        elif(swath_type == SwathTypeSAR.FIXED and num_sub_swath>1):
            raise NotImplementedError

        elif(num_sub_swath == 1): # condition is reached for Stripmap operation or ScanSAR operation with 1 subswath
            swath_size = W_gr_obs
            
        else:
            raise RuntimeError("Unknown condition reached.")
                        
        return (f_P, swath_size)

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