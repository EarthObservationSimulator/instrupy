""" 
.. module:: radiometer_model

:synopsis: *Module to handle radiometer.*

        The radiometer model is based on the reference listed below. The following system types are supported: total-power, 
        unbalanced-Dicke, balanced-Dicke and noise-adding. The following scan types are supported: fixed (no scan), cross-track 
        and conical.

        References: [1] Chapter 6,7 in "Microwave Radar and Radiometric Remote Sensing," David Gardner Long , Fawwaz T. Ulaby 2014 

The module consists of separate classes for each of the radiometric systems (listed in ``SystemType``), all with identical interface functions.
Similarly there are separate classes to handle different scan techniques (listed in ``ScanTech``), all with identical interface functions.  
By having identical interface functions, the functions can be invoked without consideration of the underlying radiometric system class.

Definition of the predetection stage:

From Pg 273, Fig.7-13 in [1], the predetection stage includes all subsystems between the antenna and the input terminals of the square-law detector.
The specifications of the radiometric system can be made by either defining the specification of the entire predetection stage or of their individual components.

.. todo:: Field-of-view for conical-scan radiometers.

"""
import copy
import uuid
import numpy as np
from collections import namedtuple
from instrupy.util import Entity, EnumEntity, Orientation,ReferenceFrame, SphericalGeometry, ViewGeometry, Maneuver, Antenna, GeoUtilityFunctions, MathUtilityFunctions, Constants, FileUtilityFunctions

PredetectionSectionParams = namedtuple("PredetectionSectionParams", ["G", "G_p", "G_m", "T_REC_q", "B"])
""" Function returns a namedtuple class to store predetection-section parameters: 
    
    * G: predetection-section gain
    * G_p: predetection gain +
    * G_m: predetection gain -
    * T_REC_q: Predetection stage (*including* the transmission line from antenna to the RF amplifier) input noise temperature. (Receiver noise temperature referred to the antenna terminals.)
    * B: Predetection bandwidth

"""

SystemParams = namedtuple("SystemParams", ["G_s_bar", "G_s_delta", "T_A", "T_SYS"])
""" Function returns a namedtuple class to store the entire-system parameters: 
    
    * G_s_bar: Average system gain.
    * G_s_delta: System gain variation.
    * T_A: Antenna (radiometric) temperature referred at the output terminal of the antenna.
    * T_SYS: System noise temperature.

"""
class SystemType(EnumEntity):
    """Enumeration of recognized radiometer systems.
    
    :cvar TOTAL_POWER: Total-power radiometer system (no reduction of effect of receiver gain-variation on radiometric-sensitivity).
    :vartype TOTAL_POWER: str

    :cvar UNBALANCED_DICKE: Unbalanced Dicke radiometer system. 
    :vartype UNBALANCED_DICKE: str

    :cvar BALANCED_DICKE: Balanced Dicke radiometer system. 
    :vartype BALANCED_DICKE: str

    :cvar NOISE_ADDING: Noise adding radiometer system with reduction of effect of receiver gain-variation on radiometric-sensitivity. Does "not" use Dicke switch. 
    :vartype NOISE_ADDING: str
    
    """
    TOTAL_POWER = "TOTAL_POWER",
    UNBALANCED_DICKE = "UNBALANCED_DICKE",
    BALANCED_DICKE = "BALANCED_DICKE",
    NOISE_ADDING = "NOISE_ADDING"

class TotalPowerRadiometerSystem(Entity):
    """ Class to handle total power radiometer system. Refer Section 7.4, 7.5 in [1].

    :ivar tlLoss: Transmission line loss in decibels.
    :vartype tlLoss: float

    :ivar tlPhyTemp: Transmission line physical temperature in Kelvin.
    :vartype tlPhyTemp: float

    :ivar rfAmpGain: RF amplifier gain in decibels.
    :vartype rfAmpGain: float

    :ivar rfAmpInpNoiseTemp: RF amplifier input noise temperature in Kelvin.
    :vartype rfAmpInpNoiseTemp: float

    :ivar rfAmpGainVariation: RF amplifier gain variation. Linear units.
    :vartype rfAmpGainVariation: float

    :ivar mixerGain: Mixer gain in decibels.
    :vartype mixerGain: float
    
    :ivar mixerInpNoiseAmp: Mixer input noise temperature in Kelvin.
    :vartype mixerInpNoiseAmp: float

    :ivar mixerGainVariation: Mixer gain variation. Linear units.
    :vartype mixerGainVariation: float

    :ivar ifAmpGain: Intermediate frequency amplifier gain in decibels.
    :vartype ifAmpGain: float

    :ivar ifAmpInpNoiseTemp: Intermediate frequency amplifier input noise temperature in Kelvin.
    :vartype ifAmpInpNoiseTemp: float

    :ivar ifAmpGainVariation: IF amplifier gain variation. Linear units.
    :vartype ifAmpGainVariation: float

    :ivar integratorVoltageGain: Integrator voltage gain (unitless).
    :vartype integratorVoltageGain: float

    :ivar predetectionGain: Pre-detection stage gain in decibels.
    :vartype predetectionGain: float

    :ivar predetectionInpNoiseTemp: Pre-detection input noise temperature in Kelvin.
    :vartype predetectionInpNoiseTemp: float

    :ivar predetectionGainVariation: Pre-detection stage gain variation. Linear units.
    :vartype predetectionGainVariation: float

    :ivar integrationTime: Integration time in seconds.
    :vartype integrationTime: float

    :ivar bandwidth: Pre-detection bandwidth in (Hertz).
    :vartype bandwidth: float

    :ivar _id: Unique identifier.
    :vartype _id: str or int

    """
    def __init__(self, tlLoss=None, tlPhyTemp=None, rfAmpGain=None, rfAmpInpNoiseTemp=None, rfAmpGainVariation=None, 
                 mixerGain=None, mixerInpNoiseTemp=None, mixerGainVariation=None, ifAmpGain=None, ifAmpInputNoiseTemp=None, ifAmpGainVariation=None, 
                 integratorVoltageGain=None, predetectionGain=None, predetectionInpNoiseTemp=None, predetectionGainVariation=None,
                 integrationTime=None, bandwidth=None, _id=None):
        
        self.tlLoss = float(tlLoss) if tlLoss is not None else None
        self.tlPhyTemp = float(tlPhyTemp) if tlPhyTemp is not None else None
        self.rfAmpGain = float(rfAmpGain) if rfAmpGain is not None else None
        self.rfAmpInpNoiseTemp = float(rfAmpInpNoiseTemp) if rfAmpInpNoiseTemp is not None else None
        self.rfAmpGainVariation = float(rfAmpGainVariation) if rfAmpGainVariation is not None else None
        self.mixerGain = float(mixerGain) if mixerGain is not None else None
        self.mixerInpNoiseTemp = float(mixerInpNoiseTemp) if mixerInpNoiseTemp is not None else None
        self.mixerGainVariation = float(mixerGainVariation) if mixerGainVariation is not None else None
        self.ifAmpGain = float(ifAmpGain) if ifAmpGain is not None else None
        self.ifAmpInputNoiseTemp = float(ifAmpInputNoiseTemp) if ifAmpInputNoiseTemp is not None else None
        self.ifAmpGainVariation = float(ifAmpGainVariation) if ifAmpGainVariation is not None else None
        self.integratorVoltageGain = float(integratorVoltageGain) if integratorVoltageGain is not None else None
        self.predetectionGain = float(predetectionGain) if predetectionGain is not None else None
        self.predetectionInpNoiseTemp = float(predetectionInpNoiseTemp) if predetectionInpNoiseTemp is not None else None
        self.predetectionGainVariation = float(predetectionGainVariation) if predetectionGainVariation is not None else None
        self.integrationTime = float(integrationTime) if integrationTime is not None else None
        self.bandwidth = float(bandwidth) if bandwidth is not None else None

        super(TotalPowerRadiometerSystem, self).__init__(_id, "TotalPowerRadiometerSystem")
    
    @staticmethod
    def from_dict(d):
        """Parses an ``TotalPowerRadiometerSystem`` object from a normalized JSON dictionary.
        
        :param d: Dictionary with the total-power radiometer system specifications.
        :paramtype d: dict

        :return: ``TotalPowerRadiometerSystem`` object.
        :rtype: :class:`instrupy.radiometer_model.TotalPowerRadiometerSystem`

        """             
        return TotalPowerRadiometerSystem(
                tlLoss = d.get("tlLoss", None),
                tlPhyTemp = d.get("tlPhyTemp", None),
                rfAmpGain = d.get("rfAmpGain", None),
                rfAmpInpNoiseTemp = d.get("rfAmpInpNoiseTemp", None),
                rfAmpGainVariation = d.get("rfAmpGainVariation", None),
                mixerGain = d.get("mixerGain", None),
                mixerInpNoiseTemp = d.get("mixerInpNoiseTemp", None),
                mixerGainVariation = d.get("mixerGainVariation", None),
                ifAmpGain = d.get("ifAmpGain", None),
                ifAmpInputNoiseTemp = d.get("ifAmpInputNoiseTemp", None),
                ifAmpGainVariation = d.get("ifAmpGainVariation", None),
                integratorVoltageGain = d.get("integratorVoltageGain", None),
                predetectionGain = d.get("predetectionGain", None),
                predetectionInpNoiseTemp = d.get("predetectionInpNoiseTemp", None),
                predetectionGainVariation = d.get("predetectionGainVariation", None),
                integrationTime = d.get("integrationTime", None),
                bandwidth = d.get("bandwidth", None),
                _id = d.get("@id", None)
                )
    
    def to_dict(self):
        """ Translate the ``TotalPowerRadiometerSystem`` object to a Python dictionary such that it can be uniquely reconstructed back from the dictionary.
        
        :return: ``TotalPowerRadiometerSystem`` object as python dictionary
        :rtype: dict

        """
        return dict({"tlLoss": self.tlLoss,
                     "tlPhyTemp": self.tlPhyTemp,
                     "rfAmpGain": self.rfAmpGain,
                     "rfAmpInpNoiseTemp": self.rfAmpInpNoiseTemp, 
                     "rfAmpGainVariation": self.rfAmpGainVariation,
                     "mixerGain,": self.mixerGain,
                     "mixerInpNoiseTemp": self.mixerInpNoiseTemp,
                     "mixerGainVariation": self.mixerGainVariation,
                     "ifAmpGain": self.ifAmpGain,
                     "ifAmpInputNoiseTemp": self.ifAmpInputNoiseTemp,
                     "ifAmpGainVariation": self.ifAmpGainVariation,
                     "integratorVoltageGain": self.integratorVoltageGain,
                     "predetectionGain": self.predetectionGain,
                     "predetectionInpNoiseTemp": self.predetectionInpNoiseTemp,
                     "predetectionGainVariation": self.predetectionGainVariation,
                     "integrationTime": self.integrationTime,
                     "bandwidth": self.bandwidth,
                     "@id": self._id,
                     "@type": SystemType.TOTAL_POWER.value
                    })
    
    def __repr__(self):
        return "TotalPowerRadiometerSystem.from_dict({})".format(self.to_dict())

    def __eq__(self, other):
        # Equality test is simple one which compares the data attributes.
        # note that _id data attribute may be different
        if(isinstance(self, other.__class__)):
            return (self.tlLoss==other.tlLoss) and (self.tlPhyTemp==other.tlPhyTemp) and (self.rfAmpGain==other.rfAmpGain) and \
                    (self.rfAmpInpNoiseTemp==other.rfAmpInpNoiseTemp) and (self.rfAmpGainVariation==other.rfAmpGainVariation) and \
                    (self.mixerGain==other.mixerGain) and (self.mixerInpNoiseTemp==other.mixerInpNoiseTemp) and (self.mixerGainVariation==other.mixerGainVariation) and (self.ifAmpGain==other.ifAmpGain) and (self.ifAmpInputNoiseTemp==other.ifAmpInputNoiseTemp) and \
                    (self.ifAmpGainVariation==other.ifAmpGainVariation) and (self.integratorVoltageGain==other.integratorVoltageGain) and (self.predetectionGain==other.predetectionGain) and \
                    (self.predetectionInpNoiseTemp==other.predetectionInpNoiseTemp) and (self.predetectionGainVariation==other.predetectionGainVariation) and \
                    (self.integrationTime==other.integrationTime)  and (self.bandwidth==other.bandwidth)                  
        else:
            return NotImplemented
        
    @staticmethod
    def compute_integration_time(td, integration_time_spec=None):
        """ Compute integration time.

        :param td: Available dwell time (i.e. total time available for integration) in seconds per ground-pixel.
        :paramtype td: float
        
        :param integration_time_spec: Integration time specification in seconds.
        :paramtype integration_time_spec: float        

        :return: Integration time in seconds.
        :rtype: float

        """
        # initialize the integration time (t_int) to be used for the calculations
        if integration_time_spec is None:
            t_int = td
        else:
            if td < integration_time_spec: 
                t_int = td # dwell time less than specified integration time => the specified integration time cannot be achieved.
            else:
                t_int = integration_time_spec
        
        return t_int

    @staticmethod
    def compute_predetection_sec_params(predetectionBandwidth, tlLoss=None, tlPhyTemp=None, predetectionGain=None, predetectionGainVariation=None, predetectionInpNoiseTemp=None,
                                            rfAmpGain=None, mixerGain=None, ifAmpGain=None, rfAmpGainVariation=None,  mixerGainVariation=None, ifAmpGainVariation=None,
                                            rfAmpInpNoiseTemp=None, mixerInpNoiseTemp=None, ifAmpInputNoiseTemp=None
                                            ):
        """ Compute predetection section parameters.

        :param predetectionBandwidth: Pre-detection bandwidth in (Hertz).
        :paramtype predetectionBandwidth: float

        :param tlLoss: Transmission line loss in decibels.
        :paramtype tlLoss: float        

        :param tlPhyTemp: Transmission line physical temperature in Kelvin.
        :paramtype tlPhyTemp: float

        :param rfAmpGain: RF amplifier gain in decibels.
        :paramtype rfAmpGain: float

        :param rfAmpInpNoiseTemp: RF amplifier input noise temperature in Kelvin.
        :paramtype rfAmpInpNoiseTemp: float

        :param rfAmpGainVariation: RF amplifier gain variation. Linear units.
        :paramtype rfAmpGainVariation: float

        :param mixerGain: Mixer gain in decibels.
        :paramtype mixerGain: float
        
        :param mixerInpNoiseAmp: Mixer input noise temperature in Kelvin.
        :paramtype mixerInpNoiseAmp: float

        :param mixerGainVariation: Mixer gain variation. Linear units.
        :paramtype mixerGainVariation: float

        :param ifAmpGain: Intermediate frequency amplifier gain in decibels.
        :paramtype ifAmpGain: float

        :param ifAmpInpNoiseTemp: Intermediate frequency amplifier input noise temperature in Kelvin.
        :paramtype ifAmpInpNoiseTemp: float

        :param ifAmpGainVariation: IF amplifier gain variation. Linear units.
        :paramtype ifAmpGainVariation: float

        :return: Predetection section parameters. Namedtuple <gain, gain minus, gain plus, receiver (including transmission line) noise temperature, bandwidth>
        :rtype: :class:`instrupy.radiometer_model.PredetectionSectionParams` 

        """        
        # Check if the predetection section specifications are directly available, if so then use them to set the local predetection section variables.
        # Else use the individual component (of the predetection stage) specifications to calculate the predetection section parameters.

        # predetection gain
        if predetectionGain is not None:
            predetection_gain = 10**(predetectionGain/10) # convert to linear units
            predetection_gain_minus = predetection_gain - 0.5*predetectionGainVariation
            predetection_gain_plus = predetection_gain + 0.5*predetectionGainVariation
        else:
            try: #try to compute the pre-detection section gain from the component specifications
                # Fig.7-9 in [1] describes the gain of the transmission line as 1/L, where L is the transmission line loss.
                L = 10**(tlLoss/ 10)   # Transmission line loss in linear units 
                G_TL = 1/L # transmission line "gain"
                G_RF = 10**(rfAmpGain/10)
                G_MIX = 10**(mixerGain/10)
                G_IF = 10**(ifAmpGain/10)
                predetection_gain = G_TL * G_RF * G_MIX * G_IF
                predetection_gain_minus = G_TL * (G_RF - 0.5*rfAmpGainVariation) * (G_MIX - 0.5*mixerGainVariation) * (G_IF - 0.5*ifAmpGainVariation) 
                predetection_gain_plus = G_TL * (G_RF + 0.5*rfAmpGainVariation) * (G_MIX + 0.5*mixerGainVariation) * (G_IF + 0.5*ifAmpGainVariation) 
            except:
                raise RuntimeError("Missing specification of one or more component specifications in the radiometer predetection section.")        

        # predetection noise temperature
        if predetectionInpNoiseTemp is not None:
            T_REC_q = predetectionInpNoiseTemp
        else:
            try: #try to compute the predetection noise temperature from the component specifications     
                L = 10**(tlLoss/ 10)   # Transmission line loss in linear units 
                # See Secion 7-3.1 in [1] for example calculation of noise temperature from cascaded stages.                          
                G_RF = 10**(rfAmpGain/10)
                G_IF = 10**(ifAmpGain/10) 
                T_REC = rfAmpInpNoiseTemp + mixerInpNoiseTemp/ G_RF + ifAmpInputNoiseTemp/ (G_RF*G_IF)  # Eqn 7.29 in [1]
                T_REC_q = (L-1)*tlPhyTemp + L*T_REC # receiver (including the transmission line) input noise-temperature eqn 7.28 in [1]. 'q' stands for quote.)
            except:
                raise RuntimeError("Missing specification of one or more component specifications in the radiometer predetection section.")
        
        return PredetectionSectionParams(predetection_gain, predetection_gain_plus, predetection_gain_minus,  T_REC_q, predetectionBandwidth)
    
    @staticmethod
    def compute_system_params(antenna, pd_sec_params, integratorVoltageGain, T_A_q):
        """ Calculate the radiometric system parameters from the antenna, predetection section, integrator specifications and the target brightness temperature. 

        :param antenna: Antenna specifications.
        :paramtype antenna: :class:`instrupy.util.Antenna`

        :param pd_sec_params: Predetection section parameters. Namedtuple <gain, gain minus, gain plus, receiver (including transmission line) noise temperature, bandwidth>
        :paramtype pd_sec_params: :class:`instrupy.radiometer_model.PredetectionSectionParams` 

        :param integratorVoltageGain: Integrator voltage gain (unitless).
        :paramtype integratorVoltageGain: float

        :param T_A_q: Brightness temperature from target (mainlobe and sidelobes of antenna eqn 6.37 in [1]).
        :paramtype T_A_q: float

        :return: Radiometric system parameters as a namedtuple <average system gain, system gain variation, antenna radiometric temperature, system radiometric temperature>.
        :rtype: :class:`instrupy.radiometer_model.SystemParams` 

        """
        # calculate system gain factor (eqn 7.43 in [1])
        G_s = 2*integratorVoltageGain*pd_sec_params.G*Constants.Boltzmann*pd_sec_params.B        

        # calculate the system gain variation
        G_s_minus = 2*integratorVoltageGain*pd_sec_params.G_m*Constants.Boltzmann*pd_sec_params.B
        G_s_plus = 2*integratorVoltageGain*pd_sec_params.G_p*Constants.Boltzmann*pd_sec_params.B
        G_s_delta = G_s_plus - G_s_minus        

        G_s_bar = G_s # average system power gain, TODO: check

        # calculate system temperature                 
        psi = antenna.radiationEfficiency # antenna radiation efficiency = 1/ antenna loss
        T_A = psi*T_A_q + (1-psi)*antenna.phyTemp
        T_SYS = T_A + pd_sec_params.T_REC_q # eqn 7.31 in [1]

        return SystemParams(G_s_bar, G_s_delta, T_A, T_SYS)
    
    def compute_radiometric_resolution(self, td, antenna, T_A_q):
        """ Compute the radiometric resolution of a total power radiometer.

        :param td: available dwell time in seconds per ground-pixel
        :paramtype td: float

        :param antenna: Antenna specifications.
        :paramtype antenna: :class:`instrupy.util.Antenna`

        :param T_A_q: Brightness temperature from target (mainlobe and sidelobes of antenna eqn 6.37 in [1]).
        :paramtype T_A_q: float

        :return: Radiometric resolution in Kelvin.
        :rtype: float

        """        
        t_int = TotalPowerRadiometerSystem.compute_integration_time(td, self.integrationTime)       

        pd_sec_params = TotalPowerRadiometerSystem.compute_predetection_sec_params(self.bandwidth, self.tlLoss, self.tlPhyTemp, self.predetectionGain, self.predetectionGainVariation, self.predetectionInpNoiseTemp,
                                            self.rfAmpGain, self.mixerGain, self.ifAmpGain, self.rfAmpGainVariation, self.mixerGainVariation, self.ifAmpGainVariation,
                                            self.rfAmpInpNoiseTemp, self.mixerInpNoiseTemp, self.ifAmpInputNoiseTemp
                                            )   
        
        sys_params = TotalPowerRadiometerSystem.compute_system_params(antenna, pd_sec_params, self.integratorVoltageGain, T_A_q)

        # copy values to variables with shorter names
        T_SYS = sys_params.T_SYS
        G_s_bar = sys_params.G_s_bar
        G_s_delta = sys_params.G_s_delta
        B = pd_sec_params.B
        print("T_SYS",T_SYS)
        print("B",B)
        print("t_int",t_int)
        print("g_s_delta",G_s_delta)
        print("g_s_bar",G_s_bar)
        # calculate the radiometric resolution, eqn 7.58 in [1]
        resolution = T_SYS*(1/(B*t_int) + ((G_s_delta/G_s_bar)**2))**0.5
        print("resolution",resolution)
        return resolution        

class UnbalancedDikeRadiometerSystem(Entity):
    """ Class to handle unbalanced Dicke radiometer system. Refer Section 7.6 in [1].

    :ivar tlLoss: Transmission line loss in decibels.
    :vartype tlLoss: float

    :ivar tlPhyTemp: Transmission line physical temperature in Kelvin.
    :vartype tlPhyTemp: float

    :ivar rfAmpGain: RF amplifier gain in decibels.
    :vartype rfAmpGain: float

    :ivar rfAmpInpNoiseTemp: RF amplifier input noise temperature in Kelvin.
    :vartype rfAmpInpNoiseTemp: float

    :ivar rfAmpGainVariation: RF amplifier gain variation. Linear units.
    :vartype rfAmpGainVariation: float

    :ivar mixerGain: Mixer gain in decibels.
    :vartype mixerGain: float
    
    :ivar mixerInpNoiseAmp: Mixer input noise temperature in Kelvin.
    :vartype mixerInpNoiseAmp: float

    :ivar mixerGainVariation: Mixer gain variation. Linear units.
    :vartype mixerGainVariation: float

    :ivar ifAmpGain: Intermediate frequency amplifier gain in decibels.
    :vartype ifAmpGain: float

    :ivar ifAmpInpNoiseTemp: Intermediate frequency amplifier input noise temperature in Kelvin.
    :vartype ifAmpInpNoiseTemp: float

    :ivar ifAmpGainVariation: IF amplifier gain variation. Linear units.
    :vartype ifAmpGainVariation: float

    :ivar dickeSwitchOutputNoiseTemperature: Dicke switch noise temperature in Kelvin referenced to the output port.
    :vartype dickeSwitchOutputNoiseTemperature: float

    :ivar referenceTemperature: Reference source noise temperature in Kelvin.
    :vartype referenceTemperature: float

    :ivar integratorVoltageGain: Integrator voltage gain (unitless).
    :vartype integratorVoltageGain: float

    :ivar predetectionGain: Pre-detection stage gain in decibels.
    :vartype predetectionGain: float

    :ivar predetectionInpNoiseTemp: Pre-detection input noise temperature in Kelvin.
    :vartype predetectionInpNoiseTemp: float

    :ivar predetectionGainVariation: Pre-detection stage gain variation. Linear units.
    :vartype predetectionGainVariation: float

    :ivar integrationTime: Integration time in seconds.
    :vartype integrationTime: float

    :ivar bandwidth: Pre-detection bandwidth in (Hertz).
    :vartype bandwidth: float

    :ivar _id: Unique identifier.
    :vartype _id: str or int

    """
    def __init__(self, tlLoss=None, tlPhyTemp=None, rfAmpGain=None, rfAmpInpNoiseTemp=None, rfAmpGainVariation=None, 
                 mixerGain=None, mixerInpNoiseTemp=None, mixerGainVariation=None, ifAmpGain=None, ifAmpInputNoiseTemp=None, ifAmpGainVariation=None, 
                 dickeSwitchOutputNoiseTemperature=None, referenceTemperature=None,
                 integratorVoltageGain=None, predetectionGain=None, predetectionInpNoiseTemp=None, predetectionGainVariation=None,
                 integrationTime=None, bandwidth=None, _id=None):
        
        self.tlLoss = float(tlLoss) if tlLoss is not None else None
        self.tlPhyTemp = float(tlPhyTemp) if tlPhyTemp is not None else None
        self.rfAmpGain = float(rfAmpGain) if rfAmpGain is not None else None
        self.rfAmpInpNoiseTemp = float(rfAmpInpNoiseTemp) if rfAmpInpNoiseTemp is not None else None
        self.rfAmpGainVariation = float(rfAmpGainVariation) if rfAmpGainVariation is not None else None
        self.mixerGain = float(mixerGain) if mixerGain is not None else None
        self.mixerInpNoiseTemp = float(mixerInpNoiseTemp) if mixerInpNoiseTemp is not None else None
        self.mixerGainVariation = float(mixerGainVariation) if mixerGainVariation is not None else None
        self.ifAmpGain = float(ifAmpGain) if ifAmpGain is not None else None
        self.ifAmpInputNoiseTemp = float(ifAmpInputNoiseTemp) if ifAmpInputNoiseTemp is not None else None
        self.ifAmpGainVariation = float(ifAmpGainVariation) if ifAmpGainVariation is not None else None
        self.dickeSwitchOutputNoiseTemperature = float(dickeSwitchOutputNoiseTemperature) if dickeSwitchOutputNoiseTemperature is not None else None
        self.referenceTemperature = float(referenceTemperature) if referenceTemperature is not None else None
        self.integratorVoltageGain = float(integratorVoltageGain) if integratorVoltageGain is not None else None
        self.predetectionGain = float(predetectionGain) if predetectionGain is not None else None
        self.predetectionInpNoiseTemp = float(predetectionInpNoiseTemp) if predetectionInpNoiseTemp is not None else None
        self.predetectionGainVariation = float(predetectionGainVariation) if predetectionGainVariation is not None else None
        self.integrationTime = float(integrationTime) if integrationTime is not None else None
        self.bandwidth = float(bandwidth) if bandwidth is not None else None

        super(UnbalancedDikeRadiometerSystem, self).__init__(_id, "UnbalancedDikeRadiometerSystem")
    
    @staticmethod
    def from_dict(d):
        """Parses an ``UnbalancedDikeRadiometerSystem`` object from a normalized JSON dictionary.
        
        :param d: Dictionary with the unbalanced Dicke radiometer system specifications.
        :paramtype d: dict

        :return: ``UnbalancedDikeRadiometerSystem`` object.
        :rtype: :class:`instrupy.radiometer_model.UnbalancedDikeRadiometerSystem`

        """             
        return UnbalancedDikeRadiometerSystem(
                tlLoss = d.get("tlLoss", None),
                tlPhyTemp = d.get("tlPhyTemp", None),
                rfAmpGain = d.get("rfAmpGain", None),
                rfAmpInpNoiseTemp = d.get("rfAmpInpNoiseTemp", None),
                rfAmpGainVariation = d.get("rfAmpGainVariation", None),
                mixerGain = d.get("mixerGain", None),
                mixerInpNoiseTemp = d.get("mixerInpNoiseTemp", None),
                mixerGainVariation = d.get("mixerGainVariation", None),
                ifAmpGain = d.get("ifAmpGain", None),
                ifAmpInputNoiseTemp = d.get("ifAmpInputNoiseTemp", None),
                ifAmpGainVariation = d.get("ifAmpGainVariation", None),
                dickeSwitchOutputNoiseTemperature = d.get("dickeSwitchOutputNoiseTemperature", None),
                referenceTemperature = d.get("referenceTemperature", None),
                integratorVoltageGain = d.get("integratorVoltageGain", None),
                predetectionGain = d.get("predetectionGain", None),
                predetectionInpNoiseTemp = d.get("predetectionInpNoiseTemp", None),
                predetectionGainVariation = d.get("predetectionGainVariation", None),
                integrationTime = d.get("integrationTime", None),
                bandwidth = d.get("bandwidth", None),
                _id = d.get("@id", None)
                )
    
    def to_dict(self):
        """ Translate the ``UnbalancedDikeRadiometerSystem`` object to a Python dictionary such that it can be uniquely reconstructed back from the dictionary.
        
        :return: ``UnbalancedDikeRadiometerSystem`` object as python dictionary
        :rtype: dict

        """
        return dict({"tlLoss": self.tlLoss,
                     "tlPhyTemp": self.tlPhyTemp,
                     "rfAmpGain": self.rfAmpGain,
                     "rfAmpInpNoiseTemp": self.rfAmpInpNoiseTemp, 
                     "rfAmpGainVariation": self.rfAmpGainVariation,
                     "mixerGain,": self.mixerGain,
                     "mixerInpNoiseTemp": self.mixerInpNoiseTemp,
                     "mixerGainVariation": self.mixerGainVariation,
                     "ifAmpGain": self.ifAmpGain,
                     "ifAmpInputNoiseTemp": self.ifAmpInputNoiseTemp,
                     "ifAmpGainVariation": self.ifAmpGainVariation,
                     "dickeSwitchOutputNoiseTemperature":self.dickeSwitchOutputNoiseTemperature,
                     "referenceTemperature": self.referenceTemperature,
                     "integratorVoltageGain": self.integratorVoltageGain,
                     "predetectionGain": self.predetectionGain,
                     "predetectionInpNoiseTemp": self.predetectionInpNoiseTemp,
                     "predetectionGainVariation": self.predetectionGainVariation,
                     "integrationTime": self.integrationTime,
                     "bandwidth": self.bandwidth,
                     "@id": self._id,
                     "@type": SystemType.UNBALANCED_DICKE.value
                    })
    
    def __repr__(self):
        return "UnbalancedDikeRadiometerSystem.from_dict({})".format(self.to_dict())

    def __eq__(self, other):
        # Equality test is simple one which compares the data attributes.
        # note that _id data attribute may be different
        if(isinstance(self, other.__class__)):
            return (self.tlLoss==other.tlLoss) and (self.tlPhyTemp==other.tlPhyTemp) and (self.rfAmpGain==other.rfAmpGain) and \
                    (self.rfAmpInpNoiseTemp==other.rfAmpInpNoiseTemp) and (self.rfAmpGainVariation==other.rfAmpGainVariation) and \
                    (self.mixerGain==other.mixerGain) and (self.mixerInpNoiseTemp==other.mixerInpNoiseTemp) and (self.mixerGainVariation==other.mixerGainVariation) and (self.ifAmpGain==other.ifAmpGain) and (self.ifAmpInputNoiseTemp==other.ifAmpInputNoiseTemp) and \
                    (self.ifAmpGainVariation==other.ifAmpGainVariation) and (self.dickeSwitchOutputNoiseTemperature==other.dickeSwitchOutputNoiseTemperature) and (self.referenceTemperature==other.referenceTemperature) and \
                    (self.integratorVoltageGain==other.integratorVoltageGain) and (self.predetectionGain==other.predetectionGain) and \
                    (self.predetectionInpNoiseTemp==other.predetectionInpNoiseTemp) and (self.predetectionGainVariation==other.predetectionGainVariation) and \
                    (self.integrationTime==other.integrationTime)  and (self.bandwidth==other.bandwidth)         
        else:
            return NotImplemented

    def compute_radiometric_resolution(self, td, antenna, T_A_q):
        """ Compute the radiometric resolution of an unbalanced Dicke radiometer.
        
        Note that the predetection stage in Dicke radiometers include the noise from the Dicke switch referenced to the switch output port.

        :param td: available dwell time in seconds per ground-pixel
        :paramtype td: float

        :param antenna: Antenna specifications.
        :paramtype antenna: :class:`instrupy.util.Antenna`

        :param T_A_q: Brightness temperature from target (mainlobe and sidelobes of antenna eqn 6.37 in [1]).
        :paramtype T_A_q: float

        :return: Radiometric resolution in Kelvin.
        :rtype: float

        """
        t_int = TotalPowerRadiometerSystem.compute_integration_time(td, self.integrationTime)

        pd_sec_params = TotalPowerRadiometerSystem.compute_predetection_sec_params(self.bandwidth, self.tlLoss, self.tlPhyTemp, self.predetectionGain, self.predetectionGainVariation, self.predetectionInpNoiseTemp,
                                            self.rfAmpGain, self.mixerGain, self.ifAmpGain, self.rfAmpGainVariation, self.mixerGainVariation, self.ifAmpGainVariation,
                                            self.rfAmpInpNoiseTemp, self.mixerInpNoiseTemp, self.ifAmpInputNoiseTemp
                                            )
        if self.dickeSwitchOutputNoiseTemperature: # self.dickeSwitchOutputNoiseTemperature could be none in case the entire predetection system parameters are specified instead of the component-wise specification.
            pd_sec_params = pd_sec_params._replace(T_REC_q=pd_sec_params.T_REC_q + self.dickeSwitchOutputNoiseTemperature) # add the Dicke switch output noise temperature
        sys_params = TotalPowerRadiometerSystem.compute_system_params(antenna, pd_sec_params, self.integratorVoltageGain, T_A_q)
        # copy values to variables with shorter names
        T_REF = self.referenceTemperature
        T_SYS = sys_params.T_SYS
        G_s_bar = sys_params.G_s_bar
        G_s_delta = sys_params.G_s_delta
        T_A = sys_params.T_A
        T_REC_q = pd_sec_params.T_REC_q
        B = pd_sec_params.B
        
        # calculate the radiometric resolution, eqn 7.68 in [1]
        resolution = (2*T_SYS**2 + 2*(T_REF + T_REC_q)**2)/(B*t_int) + ((G_s_delta/G_s_bar)**2)*(T_A - T_REF)**2
        resolution = resolution**0.5
        return resolution

class BalancedDikeRadiometerSystem(Entity):
    """ Class to handle balanced Dicke radiometer System. Refer Section 7.6 and 7.7 in [1].

    :ivar tlLoss: Transmission line loss in decibels.
    :vartype tlLoss: float

    :ivar tlPhyTemp: Transmission line physical temperature in Kelvin.
    :vartype tlPhyTemp: float

    :ivar rfAmpGain: RF amplifier gain in decibels.
    :vartype rfAmpGain: float

    :ivar rfAmpInpNoiseTemp: RF amplifier input noise temperature in Kelvin.
    :vartype rfAmpInpNoiseTemp: float

    :ivar rfAmpGainVariation: RF amplifier gain variation. Linear units.
    :vartype rfAmpGainVariation: float

    :ivar mixerGain: Mixer gain in decibels.
    :vartype mixerGain: float
    
    :ivar mixerInpNoiseAmp: Mixer input noise temperature in Kelvin.
    :vartype mixerInpNoiseAmp: float

    :ivar mixerGainVariation: Mixer gain variation. Linear units.
    :vartype mixerGainVariation: float

    :ivar ifAmpGain: Intermediate frequency amplifier gain in decibels.
    :vartype ifAmpGain: float

    :ivar ifAmpInpNoiseTemp: Intermediate frequency amplifier input noise temperature in Kelvin.
    :vartype ifAmpInpNoiseTemp: float

    :ivar ifAmpGainVariation: IF amplifier gain variation. Linear units.
    :vartype ifAmpGainVariation: float

    :ivar dickeSwitchOutputNoiseTemperature: Dicke switch noise temperature in Kelvin referenced to the output port.
    :vartype dickeSwitchOutputNoiseTemperature: float

    :ivar integratorVoltageGain: Integrator voltage gain (unitless).
    :vartype integratorVoltageGain: float

    :ivar predetectionGain: Pre-detection stage gain in decibels.
    :vartype predetectionGain: float

    :ivar predetectionInpNoiseTemp: Pre-detection input noise temperature in Kelvin.
    :vartype predetectionInpNoiseTemp: float

    :ivar predetectionGainVariation: Pre-detection stage gain variation. Linear units.
    :vartype predetectionGainVariation: float

    :ivar integrationTime: Integration time in seconds.
    :vartype integrationTime: float

    :ivar bandwidth: Pre-detection bandwidth in (Hertz).
    :vartype bandwidth: float

    :ivar _id: Unique identifier.
    :vartype _id: str or int

    """
    def __init__(self, tlLoss=None, tlPhyTemp=None, rfAmpGain=None, rfAmpInpNoiseTemp=None, rfAmpGainVariation=None, 
                 mixerGain=None, mixerInpNoiseTemp=None, mixerGainVariation=None, ifAmpGain=None, ifAmpInputNoiseTemp=None, ifAmpGainVariation=None, 
                 dickeSwitchOutputNoiseTemperature=None,
                 integratorVoltageGain=None, predetectionGain=None, predetectionInpNoiseTemp=None, predetectionGainVariation=None,
                 integrationTime=None, bandwidth=None, _id=None):
        
        self.tlLoss = float(tlLoss) if tlLoss is not None else None
        self.tlPhyTemp = float(tlPhyTemp) if tlPhyTemp is not None else None
        self.rfAmpGain = float(rfAmpGain) if rfAmpGain is not None else None
        self.rfAmpInpNoiseTemp = float(rfAmpInpNoiseTemp) if rfAmpInpNoiseTemp is not None else None
        self.rfAmpGainVariation = float(rfAmpGainVariation) if rfAmpGainVariation is not None else None
        self.mixerGain = float(mixerGain) if mixerGain is not None else None
        self.mixerInpNoiseTemp = float(mixerInpNoiseTemp) if mixerInpNoiseTemp is not None else None
        self.mixerGainVariation = float(mixerGainVariation) if mixerGainVariation is not None else None
        self.ifAmpGain = float(ifAmpGain) if ifAmpGain is not None else None
        self.ifAmpInputNoiseTemp = float(ifAmpInputNoiseTemp) if ifAmpInputNoiseTemp is not None else None
        self.ifAmpGainVariation = float(ifAmpGainVariation) if ifAmpGainVariation is not None else None
        self.dickeSwitchOutputNoiseTemperature = float(dickeSwitchOutputNoiseTemperature) if dickeSwitchOutputNoiseTemperature is not None else None
        self.integratorVoltageGain = float(integratorVoltageGain) if integratorVoltageGain is not None else None
        self.predetectionGain = float(predetectionGain) if predetectionGain is not None else None
        self.predetectionInpNoiseTemp = float(predetectionInpNoiseTemp) if predetectionInpNoiseTemp is not None else None
        self.predetectionGainVariation = float(predetectionGainVariation) if predetectionGainVariation is not None else None
        self.integrationTime = float(integrationTime) if integrationTime is not None else None
        self.bandwidth = float(bandwidth) if bandwidth is not None else None

        super(BalancedDikeRadiometerSystem, self).__init__(_id, "BalancedDikeRadiometerSystem")
    
    @staticmethod
    def from_dict(d):
        """Parses an ``BalancedDikeRadiometerSystem`` object from a normalized JSON dictionary.
        
        :param d: Dictionary with the balanced Dicke radiometer system specifications.
        :paramtype d: dict

        :return: ``BalancedDikeRadiometerSystem`` object.
        :rtype: :class:`instrupy.radiometer_model.BalancedDikeRadiometerSystem`

        """             
        return BalancedDikeRadiometerSystem(
                tlLoss = d.get("tlLoss", None),
                tlPhyTemp = d.get("tlPhyTemp", None),
                rfAmpGain = d.get("rfAmpGain", None),
                rfAmpInpNoiseTemp = d.get("rfAmpInpNoiseTemp", None),
                rfAmpGainVariation = d.get("rfAmpGainVariation", None),
                mixerGain = d.get("mixerGain", None),
                mixerInpNoiseTemp = d.get("mixerInpNoiseTemp", None),
                mixerGainVariation = d.get("mixerGainVariation", None),
                ifAmpGain = d.get("ifAmpGain", None),
                ifAmpInputNoiseTemp = d.get("ifAmpInputNoiseTemp", None),
                ifAmpGainVariation = d.get("ifAmpGainVariation", None),
                dickeSwitchOutputNoiseTemperature = d.get("dickeSwitchOutputNoiseTemperature", None),
                integratorVoltageGain = d.get("integratorVoltageGain", None),
                predetectionGain = d.get("predetectionGain", None),
                predetectionInpNoiseTemp = d.get("predetectionInpNoiseTemp", None),
                predetectionGainVariation = d.get("predetectionGainVariation", None),
                integrationTime = d.get("integrationTime", None),
                bandwidth = d.get("bandwidth", None),
                _id = d.get("@id", None)
                )
    
    def to_dict(self):
        """ Translate the ``BalancedDikeRadiometerSystem`` object to a Python dictionary such that it can be uniquely reconstructed back from the dictionary.
        
        :return: ``BalancedDikeRadiometerSystem`` object as python dictionary
        :rtype: dict

        """
        return dict({"tlLoss": self.tlLoss,
                     "tlPhyTemp": self.tlPhyTemp,
                     "rfAmpGain": self.rfAmpGain,
                     "rfAmpInpNoiseTemp": self.rfAmpInpNoiseTemp, 
                     "rfAmpGainVariation": self.rfAmpGainVariation,
                     "mixerGain,": self.mixerGain,
                     "mixerInpNoiseTemp": self.mixerInpNoiseTemp,
                     "mixerGainVariation": self.mixerGainVariation,
                     "ifAmpGain": self.ifAmpGain,
                     "ifAmpInputNoiseTemp": self.ifAmpInputNoiseTemp,
                     "ifAmpGainVariation": self.ifAmpGainVariation,
                     "dickeSwitchOutputNoiseTemperature":self.dickeSwitchOutputNoiseTemperature,
                     "integratorVoltageGain": self.integratorVoltageGain,
                     "predetectionGain": self.predetectionGain,
                     "predetectionInpNoiseTemp": self.predetectionInpNoiseTemp,
                     "predetectionGainVariation": self.predetectionGainVariation,
                     "integrationTime": self.integrationTime,
                     "bandwidth": self.bandwidth,
                     "@id": self._id,
                     "@type": SystemType.BALANCED_DICKE.value
                    })
    
    def __repr__(self):
        return "BalancedDikeRadiometerSystem.from_dict({})".format(self.to_dict())

    def __eq__(self, other):
        # Equality test is simple one which compares the data attributes.
        # note that _id data attribute may be different
        if(isinstance(self, other.__class__)):
            return (self.tlLoss==other.tlLoss) and (self.tlPhyTemp==other.tlPhyTemp) and (self.rfAmpGain==other.rfAmpGain) and \
                    (self.rfAmpInpNoiseTemp==other.rfAmpInpNoiseTemp) and (self.rfAmpGainVariation==other.rfAmpGainVariation) and \
                    (self.mixerGain==other.mixerGain) and (self.mixerInpNoiseTemp==other.mixerInpNoiseTemp) and (self.mixerGainVariation==other.mixerGainVariation) and (self.ifAmpGain==other.ifAmpGain) and (self.ifAmpInputNoiseTemp==other.ifAmpInputNoiseTemp) and \
                    (self.ifAmpGainVariation==other.ifAmpGainVariation) and (self.dickeSwitchOutputNoiseTemperature==other.dickeSwitchOutputNoiseTemperature) and \
                    (self.integratorVoltageGain==other.integratorVoltageGain) and (self.predetectionGain==other.predetectionGain) and \
                    (self.predetectionInpNoiseTemp==other.predetectionInpNoiseTemp) and (self.predetectionGainVariation==other.predetectionGainVariation)  and \
                    (self.integrationTime==other.integrationTime)  and (self.bandwidth==other.bandwidth)         
        else:
            return NotImplemented
    
    def compute_radiometric_resolution(self, td, antenna, T_A_q):
        """ Compute the radiometric resolution of an unbalanced Dicke radiometer.

        :param td: available dwell time in seconds per ground-pixel
        :paramtype td: float

        :param antenna: Antenna specifications.
        :paramtype antenna: :class:`instrupy.util.Antenna`

        :param T_A_q: Brightness temperature from target (mainlobe and sidelobes of antenna eqn 6.37 in [1]).
        :paramtype T_A_q: float

        :return: Radiometric resolution in Kelvin.
        :rtype: float

        """
        t_int = TotalPowerRadiometerSystem.compute_integration_time(td, self.integrationTime)
        
        pd_sec_params = TotalPowerRadiometerSystem.compute_predetection_sec_params(self.bandwidth, self.tlLoss, self.tlPhyTemp, self.predetectionGain, self.predetectionGainVariation, self.predetectionInpNoiseTemp,
                                    self.rfAmpGain, self.mixerGain, self.ifAmpGain, self.rfAmpGainVariation, self.mixerGainVariation, self.ifAmpGainVariation,
                                    self.rfAmpInpNoiseTemp, self.mixerInpNoiseTemp, self.ifAmpInputNoiseTemp
                                    )
        if self.dickeSwitchOutputNoiseTemperature: # self.dickeSwitchOutputNoiseTemperature could be none in case the entire predetection system parameters are specified instead of the component-wise specification.
            pd_sec_params = pd_sec_params._replace(T_REC_q=pd_sec_params.T_REC_q + self.dickeSwitchOutputNoiseTemperature) # add the Dicke switch output noise temperature

        sys_params = TotalPowerRadiometerSystem.compute_system_params(antenna, pd_sec_params, self.integratorVoltageGain, T_A_q)
        
        # copy values to variables with shorter names
        T_SYS = sys_params.T_SYS
        B = pd_sec_params.B

        # calculate the radiometric resolution, eqn 7.69 in [1]
        resolution = 2*T_SYS/(B*t_int)**0.5

        return resolution

class NoiseAddingRadiometerSystem(Entity):
    """ Class to handle noise-adding radiometer system. Refer Section 7.9 in [1].

    :ivar tlLoss: Transmission line loss in decibels.
    :vartype tlLoss: float

    :ivar tlPhyTemp: Transmission line physical temperature in Kelvin.
    :vartype tlPhyTemp: float

    :ivar rfAmpGain: RF amplifier gain in decibels.
    :vartype rfAmpGain: float

    :ivar rfAmpInpNoiseTemp: RF amplifier input noise temperature in Kelvin.
    :vartype rfAmpInpNoiseTemp: float

    :ivar rfAmpGainVariation: RF amplifier gain variation. Linear units.
    :vartype rfAmpGainVariation: float

    :ivar mixerGain: Mixer gain in decibels.
    :vartype mixerGain: float
    
    :ivar mixerInpNoiseAmp: Mixer input noise temperature in Kelvin.
    :vartype mixerInpNoiseAmp: float

    :ivar mixerGainVariation: Mixer gain variation. Linear units.
    :vartype mixerGainVariation: float

    :ivar ifAmpGain: Intermediate frequency amplifier gain in decibels.
    :vartype ifAmpGain: float

    :ivar ifAmpInpNoiseTemp: Intermediate frequency amplifier input noise temperature in Kelvin.
    :vartype ifAmpInpNoiseTemp: float

    :ivar ifAmpGainVariation: IF amplifier gain variation. Linear units.
    :vartype ifAmpGainVariation: float

    :ivar excessNoiseTemperature: Excess noise temperature (added noise to the receiver input during the diode ON half-cycle) in Kelvin referenced to the output port.
    :vartype excessNoiseTemperature: float

    :ivar integratorVoltageGain: Integrator voltage gain (unitless).
    :vartype integratorVoltageGain: float

    :ivar predetectionGain: Pre-detection stage gain in decibels.
    :vartype predetectionGain: float

    :ivar predetectionInpNoiseTemp: Pre-detection input noise temperature in Kelvin.
    :vartype predetectionInpNoiseTemp: float

    :ivar predetectionGainVariation: Pre-detection stage gain variation. Linear units.
    :vartype predetectionGainVariation: float

    :ivar integrationTime: Integration time in seconds.
    :vartype integrationTime: float

    :ivar bandwidth: Pre-detection bandwidth in (Hertz).
    :vartype bandwidth: float

    :ivar _id: Unique identifier.
    :vartype _id: str or int

    """
    def __init__(self, tlLoss=None, tlPhyTemp=None, rfAmpGain=None, rfAmpInpNoiseTemp=None, rfAmpGainVariation=None, 
                 mixerGain=None, mixerInpNoiseTemp=None, mixerGainVariation=None, ifAmpGain=None, ifAmpInputNoiseTemp=None, ifAmpGainVariation=None, 
                 excessNoiseTemperature=None,
                 integratorVoltageGain=None, predetectionGain=None, predetectionInpNoiseTemp=None, predetectionGainVariation=None,
                 integrationTime=None, bandwidth=None, _id=None):
        
        self.tlLoss = float(tlLoss) if tlLoss is not None else None
        self.tlPhyTemp = float(tlPhyTemp) if tlPhyTemp is not None else None
        self.rfAmpGain = float(rfAmpGain) if rfAmpGain is not None else None
        self.rfAmpInpNoiseTemp = float(rfAmpInpNoiseTemp) if rfAmpInpNoiseTemp is not None else None
        self.rfAmpGainVariation = float(rfAmpGainVariation) if rfAmpGainVariation is not None else None
        self.mixerGain = float(mixerGain) if mixerGain is not None else None
        self.mixerInpNoiseTemp = float(mixerInpNoiseTemp) if mixerInpNoiseTemp is not None else None
        self.mixerGainVariation = float(mixerGainVariation) if mixerGainVariation is not None else None
        self.ifAmpGain = float(ifAmpGain) if ifAmpGain is not None else None
        self.ifAmpInputNoiseTemp = float(ifAmpInputNoiseTemp) if ifAmpInputNoiseTemp is not None else None
        self.ifAmpGainVariation = float(ifAmpGainVariation) if ifAmpGainVariation is not None else None
        self.excessNoiseTemperature = float(excessNoiseTemperature) if excessNoiseTemperature is not None else None
        self.integratorVoltageGain = float(integratorVoltageGain) if integratorVoltageGain is not None else None
        self.predetectionGain = float(predetectionGain) if predetectionGain is not None else None
        self.predetectionInpNoiseTemp = float(predetectionInpNoiseTemp) if predetectionInpNoiseTemp is not None else None
        self.predetectionGainVariation = float(predetectionGainVariation) if predetectionGainVariation is not None else None
        self.integrationTime = float(integrationTime) if integrationTime is not None else None
        self.bandwidth = float(bandwidth) if bandwidth is not None else None

        super(NoiseAddingRadiometerSystem, self).__init__(_id, "NoiseAddingRadiometerSystem")
    
    @staticmethod
    def from_dict(d):
        """Parses an ``NoiseAddingRadiometerSystem`` object from a normalized JSON dictionary.
        
        :param d: Dictionary with the noise-adding radiometer system specifications.
        :paramtype d: dict

        :return: ``NoiseAddingRadiometerSystem`` object.
        :rtype: :class:`instrupy.util.NoiseAddingRadiometerSystem`

        """             
        return NoiseAddingRadiometerSystem(
                tlLoss = d.get("tlLoss", None),
                tlPhyTemp = d.get("tlPhyTemp", None),
                rfAmpGain = d.get("rfAmpGain", None),
                rfAmpInpNoiseTemp = d.get("rfAmpInpNoiseTemp", None),
                rfAmpGainVariation = d.get("rfAmpGainVariation", None),
                mixerGain = d.get("mixerGain", None),
                mixerInpNoiseTemp = d.get("mixerInpNoiseTemp", None),
                mixerGainVariation = d.get("mixerGainVariation", None),
                ifAmpGain = d.get("ifAmpGain", None),
                ifAmpInputNoiseTemp = d.get("ifAmpInputNoiseTemp", None),
                ifAmpGainVariation = d.get("ifAmpGainVariation", None),
                excessNoiseTemperature = d.get("excessNoiseTemperature", None),
                integratorVoltageGain = d.get("integratorVoltageGain", None),
                predetectionGain = d.get("predetectionGain", None),
                predetectionInpNoiseTemp = d.get("predetectionInpNoiseTemp", None),
                predetectionGainVariation = d.get("predetectionGainVariation", None),
                integrationTime  = d.get("integrationTime", None),
                bandwidth = d.get("bandwidth", None),
                _id = d.get("@id", None)
                )
    
    def to_dict(self):
        """ Translate the ``NoiseAddingRadiometerSystem`` object to a Python dictionary such that it can be uniquely reconstructed back from the dictionary.
        
        :return: ``NoiseAddingRadiometerSystem`` object as python dictionary
        :rtype: dict

        """
        return dict({"tlLoss": self.tlLoss,
                     "tlPhyTemp": self.tlPhyTemp,
                     "rfAmpGain": self.rfAmpGain,
                     "rfAmpInpNoiseTemp": self.rfAmpInpNoiseTemp, 
                     "rfAmpGainVariation": self.rfAmpGainVariation,
                     "mixerGain,": self.mixerGain,
                     "mixerInpNoiseTemp": self.mixerInpNoiseTemp,
                     "mixerGainVariation": self.mixerGainVariation,
                     "ifAmpGain": self.ifAmpGain,
                     "ifAmpInputNoiseTemp": self.ifAmpInputNoiseTemp,
                     "ifAmpGainVariation": self.ifAmpGainVariation,
                     "excessNoiseTemperature":self.excessNoiseTemperature,
                     "integratorVoltageGain": self.integratorVoltageGain,
                     "predetectionGain": self.predetectionGain,
                     "predetectionInpNoiseTemp": self.predetectionInpNoiseTemp,
                     "predetectionGainVariation": self.predetectionGainVariation,
                     "integrationTime": self.integrationTime,
                     "bandwidth": self.bandwidth,
                     "@id": self._id,
                     "@type": SystemType.NOISE_ADDING.value
                    })
    
    def __repr__(self):
        return "NoiseAddingRadiometerSystem.from_dict({})".format(self.to_dict())

    def __eq__(self, other):
        # Equality test is simple one which compares the data attributes.
        # note that _id data attribute may be different
        if(isinstance(self, other.__class__)):
            return (self.tlLoss==other.tlLoss) and (self.tlPhyTemp==other.tlPhyTemp) and (self.rfAmpGain==other.rfAmpGain) and \
                    (self.rfAmpInpNoiseTemp==other.rfAmpInpNoiseTemp) and (self.rfAmpGainVariation==other.rfAmpGainVariation) and \
                    (self.mixerGain==other.mixerGain) and (self.mixerInpNoiseTemp==other.mixerInpNoiseTemp) and (self.mixerGainVariation==other.mixerGainVariation) and (self.ifAmpGain==other.ifAmpGain) and (self.ifAmpInputNoiseTemp==other.ifAmpInputNoiseTemp) and \
                    (self.ifAmpGainVariation==other.ifAmpGainVariation) and (self.excessNoiseTemperature==other.excessNoiseTemperature) and \
                    (self.integratorVoltageGain==other.integratorVoltageGain) and (self.predetectionGain==other.predetectionGain) and \
                    (self.predetectionInpNoiseTemp==other.predetectionInpNoiseTemp) and (self.predetectionGainVariation==other.predetectionGainVariation) and \
                    (self.integrationTime==other.integrationTime)  and (self.bandwidth==other.bandwidth)         
        else:
            return NotImplemented
    
    def compute_radiometric_resolution(self, td, antenna, T_A_q):
        """ Compute the radiometric resolution of a total power radiometer.

        :param td: available dwell time in seconds per ground-pixel
        :paramtype td: float

        :param antenna: Antenna specifications.
        :paramtype antenna: :class:`instrupy.util.Antenna`

        :param T_A_q: Brightness temperature from target (mainlobe and sidelobes of antenna eqn 6.37 in [1]).
        :paramtype T_A_q: float

        :return: Radiometric resolution in Kelvin.
        :rtype: float

        """
        t_int = TotalPowerRadiometerSystem.compute_integration_time(td, self.integrationTime)
        
        pd_sec_params = TotalPowerRadiometerSystem.compute_predetection_sec_params(self.bandwidth, self.tlLoss, self.tlPhyTemp, self.predetectionGain, self.predetectionGainVariation, self.predetectionInpNoiseTemp,
                                    self.rfAmpGain, self.mixerGain, self.ifAmpGain, self.rfAmpGainVariation, self.mixerGainVariation, self.ifAmpGainVariation,
                                    self.rfAmpInpNoiseTemp, self.mixerInpNoiseTemp, self.ifAmpInputNoiseTemp
                                    )

        sys_params = TotalPowerRadiometerSystem.compute_system_params(antenna, pd_sec_params, self.integratorVoltageGain, T_A_q)
        
        # copy values to variables with shorter names
        T_SYS = sys_params.T_SYS
        B = pd_sec_params.B

        # calculate the radiometric resolution, eqn 7.87 in [1]
        resolution = 2*T_SYS/(B*t_int)**0.5 * (1 + 2*T_SYS/self.excessNoiseTemperature)

        return resolution   

class ScanTech(EnumEntity):
    """Enumeration of recognized radiometer scanning techniques. See Section 7.12 in [1].
    
    :cvar FIXED: No scan.
    :vartype FIXED: str

    :cvar CROSS_TRACK: Scan along the cross-track direction.
    :vartype CROSS_TRACK: str

    :cvar CONICAL: Scan with a fixed of-nadir angle from one clock angle to another.
    :vartype CONICAL: str
    
    """
    FIXED = "FIXED",
    CROSS_TRACK = "CROSS_TRACK",
    CONICAL = "CONICAL"

class FixedScan(Entity):
    """ Class to handle fixed scan.

    :ivar _id: Unique identifier.
    :vartype _id: str or int

    """
    def __init__(self, _id=None):

        super(FixedScan, self).__init__(_id, "FixedScan")
    
    @staticmethod
    def from_dict(d):
        """Parses an ``FixedScan`` object from a normalized JSON dictionary.
        
        :param d: Dictionary with the fixed-scan specifications.
        :paramtype d: dict

        :return: ``FixedScan`` object.
        :rtype: :class:`instrupy.radiometer_model.FixedScan`

        """             
        return FixedScan(_id = d.get("@id", None))
    
    def to_dict(self):
        """ Translate the ``FixedScan`` object to a Python dictionary such that it can be uniquely reconstructed back from the dictionary.
        
        :return: ``FixedScan`` object as python dictionary
        :rtype: dict

        """
        return dict({"@id": self._id,  "@type": ScanTech.FIXED.value})
    
    def __repr__(self):
        return "FixedScan.from_dict({})".format(self.to_dict())

    def __eq__(self, other):
        # Equality test is simple one which compares the data attributes.
        # note that _id data attribute may be different
        if(isinstance(self, other.__class__)):
            return True
        else:
            return NotImplemented
    
    def compute_instru_field_of_view(self, antenna_fov_sph_geom, instru_orientation):
        """ Compute the instrument field-of-view from the antenna FOV spherical geometry and the instrument orientation 
            specifications. The instrument FOV depends on the scan-type and parameters.
        
            :param antenna_fov_sph_geom: Antenna FOV spherical geometry specification.
            :paramtype antenna_fov_sph_geom: :class:`instrupy.util.SphericalGeometry`

            :param instru_orientation: Orientation of the instrument.
            :paramtype instru_orientation: :class:`instrupy.util.Orientation`

            :return: Instrument field-of-view.
            :rtype: :class:`instrupy.util.ViewGeometry`
        
        """
        # for FIXED scan-type, the instrument FOV spherical geometry = antenna FOV spherical geometry
        return ViewGeometry(orien=instru_orientation, sph_geom=antenna_fov_sph_geom)

    def compute_dwell_time_per_ground_pixel(self, res_AT_m, sat_speed_kmps, **kwargs):
        """ Get the available dwell time per ground-pixel. The integration time 
            is set to be around the dwell time.

        :param res_AT_m: Along track pixel resolution in meters.
        :paramtype res_AT_m: float

        :param sat_speed_kmps: Satellite speed in kilometers per second.
        :paramtype sat_speed_kmps: float

        :return: Ground-pixel dwell time.
        :rtype: float

        """
        return res_AT_m/(sat_speed_kmps*1e3)

    def compute_swath_width(self, alt_km, instru_look_angle_deg, antenna_fov_sph_geom):
        """ Obtain the swath-width.
            In case of fixed-scan mode, there is only 1 imaged ground-pixel per swath. 
            Swath-width is computed to be equal to the antenna-footprint cross-track size. 
            See Fig.5.1.3.1 in Spaceborne SAR Study: LDRD 92 Final Report SANDIA Report March 1993.

        :param alt_km: Altitude of observer in kilometers.
        :paramtype alt_km: float

        :param instru_look_angle_deg: Instrument look angle in degrees. This corresponds to the off-nadir angle at which the ground-pixel is imaged.
        :paramtype instru_look_angle_deg: float

        :param antenna_fov_sph_geom: Antenna FOV spherical geometry specification.
        :paramtype antenna_fov_sph_geom: :class:`instrupy.util.SphericalGeometry`

        :return: Swath-width in kilometers.
        :rtype: float

        """
        h = alt_km * 1e3
        Re = Constants.radiusOfEarthInKM * 1e3         
        Rs = Re + h 
        gamma_m = np.deg2rad(instru_look_angle_deg)

        if antenna_fov_sph_geom.shape == SphericalGeometry.Shape.RECTANGULAR:
            iFOV_CT = np.deg2rad(antenna_fov_sph_geom.angle_width)
        elif antenna_fov_sph_geom.shape == SphericalGeometry.Shape.CIRCULAR:
            iFOV_CT = np.deg2rad(antenna_fov_sph_geom.diameter)
        else:
            raise NotImplementedError

        gamma_n_illum = gamma_m - 0.5*iFOV_CT
        gamma_f_illum = gamma_m + 0.5*iFOV_CT
        theta_horizon = np.arcsin(Re/Rs)
        
        theta_in_illum = np.arcsin(np.sin(abs(gamma_n_illum))*Rs/Re)
        if np.isnan(theta_in_illum):            
            theta_in_illum = theta_horizon # beyond horizon, hence set to horizon angle

        theta_if_illum = np.arcsin(np.sin(abs(gamma_f_illum))*Rs/Re)
        if np.isnan(theta_if_illum):
            theta_if_illum = theta_horizon # beyond horizon, hence set to horizon angle

        alpha_n_illum = theta_in_illum - abs(gamma_n_illum)
        alpha_f_illum = theta_if_illum - abs(gamma_f_illum)

        if gamma_n_illum <= 0: # NOT completely sidelooking, footprint falls on the nadir-direction as well.
            alpha_s_illum = alpha_f_illum + alpha_n_illum
        else:
            alpha_s_illum = abs(alpha_f_illum - alpha_n_illum)

        W_gr = Re*alpha_s_illum  # swath

        return W_gr*1e-3

class CrossTrackScan(Entity):
    """ Class to handle cross-track scan.

    :ivar scanWidth: Angular scan-width in degrees.
    :vartype scanWidth: float

    :ivar interScanOverheadTime: Time in seconds taken from ending current scan to starting next scan. Significant in case of mechanical scanning.
    :vartype interScanOverheadTime: float

    :ivar _id: Unique identifier.
    :vartype _id: str or int

    """
    def __init__(self, scanWidth=None, interScanOverheadTime=None, _id=None):

        self.scanWidth = float(scanWidth) if scanWidth else None
        self.interScanOverheadTime = float(interScanOverheadTime) if interScanOverheadTime else None
        super(CrossTrackScan, self).__init__(_id, "CrossTrackScan")
    
    @staticmethod
    def from_dict(d):
        """Parses an ``CrossTrackScan`` object from a normalized JSON dictionary.
        
        :param d: Dictionary with the cross-track scan specifications.

        The following default values are assigned to the object instance parameters in case of 
        :class:`None` values or missing key/value pairs in the input dictionary.

        .. csv-table:: Default values
            :header: Parameter, Default Value
            :widths: 10,40

            interScanOverheadTime, 0

        :paramtype d: dict

        :return: ``CrossTrackScan`` object.
        :rtype: :class:`instrupy.radiometer_model.CrossTrackScan`

        """             
        return CrossTrackScan( scanWidth = d.get("scanWidth", None),
                               interScanOverheadTime = d.get("interScanOverheadTime", 0),
                               _id = d.get("@id", None),)
    
    def to_dict(self):
        """ Translate the ``CrossTrackScan`` object to a Python dictionary such that it can be uniquely reconstructed back from the dictionary.
        
        :return: ``CrossTrackScan`` object as python dictionary.
        :rtype: dict

        """
        return dict({"scanWidth": self.scanWidth,
                     "interScanOverheadTime": self.interScanOverheadTime,
                     "@id": self._id,
                     "@type": ScanTech.CROSS_TRACK.value
                    })
    
    def __repr__(self):
        return "CrossTrackScan.from_dict({})".format(self.to_dict())

    def __eq__(self, other):
        # Equality test is simple one which compares the data attributes.
        # note that _id data attribute may be different
        if(isinstance(self, other.__class__)):
            return (self.scanWidth==other.scanWidth) and (self.interScanOverheadTime==other.interScanOverheadTime)
        else:
            return NotImplemented
    
    def compute_instru_field_of_view(self, antenna_fov_sph_geom, instru_orientation):
        """ Compute the instrument field-of-view from the antenna FOV spherical geometry and the instrument orientation 
            specifications. The instrument FOV depends on the scan-type and parameters.
        
            :param antenna_fov_sph_geom: Antenna FOV spherical geometry specification.
            :paramtype antenna_fov_sph_geom: :class:`instrupy.util.SphericalGeometry`

            :param instru_orientation: Orientation of the instrument.
            :paramtype instru_orientation: :class:`instrupy.util.Orientation`

            :return: Instrument field-of-view.
            :rtype: :class:`instrupy.util.ViewGeometry`
        
        """
        if antenna_fov_sph_geom.shape == SphericalGeometry.Shape.RECTANGULAR:
            iFOV_AT_deg = antenna_fov_sph_geom.angle_height
            iFOV_CT_deg = antenna_fov_sph_geom.angle_width
        elif antenna_fov_sph_geom.shape == SphericalGeometry.Shape.CIRCULAR:
            iFOV_AT_deg = antenna_fov_sph_geom.diameter
            iFOV_CT_deg = antenna_fov_sph_geom.diameter
        else:
            raise NotImplementedError

        # for CROSS_TRACK scan-type, the instrument FOV spherical geometry is RECTANGULAR shape
        instru_sph_geom = SphericalGeometry.from_rectangular_specs(angle_height=iFOV_AT_deg , angle_width= self.scanWidth + iFOV_CT_deg) 

        return ViewGeometry(orien=instru_orientation, sph_geom=instru_sph_geom)

    def compute_dwell_time_per_ground_pixel(self, res_AT_m, sat_speed_kmps, **kwargs):
        """ Get the available dwell time per ground-pixel. THe integration time 
            is set to be around the dwell time.

        :param res_AT_m: Along track pixel resolution in meters.
        :paramtype res_AT_m: float

        :param sat_speed_kmps: Satellite speed in kilometers per second.
        :paramtype sat_speed_kmps: float

        :param iFOV_CT_deg: IFOV (FOV corresponding to the ground-pixel, in degrees) in the cross-track direction.
        :paramtype iFOV_CT_deg: float

        :return: Ground-pixel dwell time.
        :rtype: float

        """
        iFOV_CT_deg = kwargs.get('iFOV_CT_deg')
        num_ground_pixels_per_strip = self.scanWidth / iFOV_CT_deg
        time_avail_for_scan = res_AT_m/(sat_speed_kmps*1e3) - self.interScanOverheadTime
        if time_avail_for_scan > 0:
            return time_avail_for_scan/num_ground_pixels_per_strip
        else:
            return 0

    def compute_swath_width(self, alt_km, instru_look_angle_deg, antenna_fov_sph_geom):
        """ Obtain the swath-width.
            In case of cross-track-scan mode, there are multiple imaged ground-pixels per swath along the cross-track.
            The instru_look_angle is assumed to correspond to a (pure) roll. 
            See Fig.5.1.3.1 in Spaceborne SAR Study: LDRD 92 Final Report SANDIA Report March 1993

        :param alt_km: Altitude of observer in kilometers.
        :paramtype alt_km: float

        :param instru_look_angle_deg: Instrument look angle in degrees. 
        :paramtype instru_look_angle_deg: float
        
        :param antenna_fov_sph_geom: Antenna FOV spherical geometry specification.
        :paramtype antenna_fov_sph_geom: :class:`instrupy.util.SphericalGeometry`

        :return: Swath-width in kilometers.
        :rtype: float

        """
        h = alt_km * 1e3
        Re = Constants.radiusOfEarthInKM * 1e3         
        Rs = Re + h 
        gamma_m = np.deg2rad(instru_look_angle_deg)

        if antenna_fov_sph_geom.shape == SphericalGeometry.Shape.RECTANGULAR:
            iFOV_CT_deg = antenna_fov_sph_geom.angle_width
        elif antenna_fov_sph_geom.shape == SphericalGeometry.Shape.CIRCULAR:
            iFOV_CT_deg = antenna_fov_sph_geom.diameter
        else:
            raise NotImplementedError

        strip_CT_rad = np.deg2rad(self.scanWidth + iFOV_CT_deg)

        gamma_n_illum = gamma_m - 0.5*strip_CT_rad
        gamma_f_illum = gamma_m + 0.5*strip_CT_rad
        theta_horizon = np.arcsin(Re/Rs)
        
        theta_in_illum = np.arcsin(np.sin(abs(gamma_n_illum))*Rs/Re)
        if np.isnan(theta_in_illum):            
            theta_in_illum = theta_horizon # beyond horizon, hence set to horizon angle

        theta_if_illum = np.arcsin(np.sin(abs(gamma_f_illum))*Rs/Re)
        if np.isnan(theta_if_illum):
            theta_if_illum = theta_horizon # beyond horizon, hence set to horizon angle

        alpha_n_illum = theta_in_illum - abs(gamma_n_illum)
        alpha_f_illum = theta_if_illum - abs(gamma_f_illum)

        if gamma_n_illum <= 0: # NOT completely sidelooking, footprint falls on the nadir-direction as well.
            alpha_s_illum = alpha_f_illum + alpha_n_illum
        else:
            alpha_s_illum = abs(alpha_f_illum - alpha_n_illum)

        W_gr = Re*abs(alpha_s_illum)  # swath

        return W_gr*1e-3

class ConicalScan(Entity):
    """ Class to handle conical scan. 
        For illustration of off-nadir angle and clock angles see Fig.7 in T. Kawanishi et al., "The Advanced Microwave Scanning Radiometer for the Earth Observing System (AMSR-E), NASDA's contribution to the EOS for global energy and water cycle studies," in IEEE Transactions on Geoscience and Remote Sensing, vol. 41, no. 2, pp. 184-194, Feb. 2003, doi: 10.1109/TGRS.2002.808331.

    :ivar offNadirAngle: Off-nadir angle (i.e. the half-cone angle of the conical scan) in degrees.
    :vartype offNadirAngle: float

    :ivar clockAngleRange: Scan clock angle range in degrees.
    :vartype clockAngleRange: float

    :ivar interScanOverheadTime: Time in seconds taken from ending current scan to starting next scan. Significant in case of mechanical scanning.
    
    The following default values are assigned to the object instance parameters in case of 
    :class:`None` values or missing key/value pairs in the input dictionary.

    .. csv-table:: Default values
        :header: Parameter, Default Value
        :widths: 10,40

        interScanOverheadTime, 0
    
    :vartype interScanOverheadTime: float

    :ivar _id: Unique identifier.
    :vartype _id: str or int

    """
    def __init__(self, offNadirAngle=None, clockAngleRange=None, interScanOverheadTime=None, _id=None):

        self.offNadirAngle = float(offNadirAngle) if offNadirAngle else None
        self.clockAngleRange = float(clockAngleRange) if clockAngleRange else None
        self.interScanOverheadTime = float(interScanOverheadTime) if interScanOverheadTime else None
        super(ConicalScan, self).__init__(_id, "ConicalScan")
    
    @staticmethod
    def from_dict(d):
        """Parses an ``ConicalScan`` object from a normalized JSON dictionary.
        
        :param d: Dictionary with the cross-track scan specifications.
        :paramtype d: dict

        :return: ``ConicalScan`` object.
        :rtype: :class:`instrupy.radiometer_model.ConicalScan`

        """             
        return ConicalScan( offNadirAngle = d.get("offNadirAngle", None),
                            clockAngleRange = d.get("clockAngleRange", None),
                            interScanOverheadTime = d.get("interScanOverheadTime", None),
                            _id = d.get("@id", None))
    
    def to_dict(self):
        """ Translate the ``ConicalScan`` object to a Python dictionary such that it can be uniquely reconstructed back from the dictionary.
        
        :return: ``ConicalScan`` object as python dictionary
        :rtype: dict

        """
        return dict({"offNadirAngle": self.offNadirAngle,
                     "clockAngleRange": self.clockAngleRange,
                     "interScanOverheadTime": self.interScanOverheadTime,
                     "@id": self._id,
                     "@type": ScanTech.CONICAL.value
                    })
    
    def __repr__(self):
        return "ConicalScan.from_dict({})".format(self.to_dict())

    def __eq__(self, other):
        # Equality test is simple one which compares the data attributes.
        # note that _id data attribute may be different
        if(isinstance(self, other.__class__)):
            return (self.offNadirAngle==other.offNadirAngle) and (self.clockAngleRange==other.clockAngleRange) and (self.interScanOverheadTime==other.interScanOverheadTime)
        else:
            return NotImplemented

    def compute_instru_field_of_view(self, antenna_fov_sph_geom, instru_orientation):
        """ Compute the instrument field-of-view from the antenna FOV spherical geometry and the instrument orientation 
            specifications. The instrument FOV depends on the scan-type and parameters.

            .. todo:: Need to implement this function. A new spherical geometry type needs to be defined (arc).
        
            :param antenna_fov_sph_geom: Antenna FOV spherical geometry specification.
            :paramtype antenna_fov_sph_geom: :class:`instrupy.util.SphericalGeometry`

            :param instru_orientation: Orientation of the instrument.
            :paramtype instru_orientation: :class:`instrupy.util.Orientation`

            :return: Instrument field-of-view.
            :rtype: :class:`instrupy.util.ViewGeometry`
        
        """
        raise NotImplementedError

    def compute_dwell_time_per_ground_pixel(self, res_AT_m, sat_speed_kmps, **kwargs):
        """ Get the available dwell time per ground-pixel. The integration time 
            is set to be around the dwell time.

            :param res_AT_m: Along track pixel resolution in meters.
            :paramtype res_AT_m: float

            :param sat_speed_kmps: Satellite speed in kilometers per second.
            :paramtype sat_speed_kmps: float

            :param iFOV_CT_deg: IFOV (FOV corresponding to the ground-pixel, in degrees) in the cross-track direction, 
            :paramtype iFOV_CT_deg: float

            :return: Ground-pixel dwell time.
            :rtype: float

        """
        iFOV_CT_deg = kwargs.get('iFOV_CT_deg')
        num_ground_pixels_per_strip = self.clockAngleRange / iFOV_CT_deg
        time_avail_for_scan = res_AT_m/(sat_speed_kmps*1e3) - self.interScanOverheadTime
        if time_avail_for_scan > 0:
            return time_avail_for_scan/num_ground_pixels_per_strip
        else:
            return 0

    def compute_swath_width(self, alt_km, instru_look_angle_deg, antenna_fov_sph_geom=None):
        """ Obtain the swath-width.
            In case of conical-scan mode, there are multiple imaged ground-pixels per swath along the scanned-strip.
            The "swath" is considered to be the length of the strip-arc. This is different from the length of the scene
            along the cross-track direction.

            :param alt_km: Altitude of observer in kilometers.
            :paramtype alt_km: float

            :param instru_look_angle_deg: Instrument look angle in degrees. This correspond to the off-nadir angle at which the ground-pixel is imaged.
            :paramtype instru_look_angle_deg: float
    
            :param antenna_fov_sph_geom: Antenna FOV spherical geometry specification. Not used in this function.
            :paramtype antenna_fov_sph_geom: :class:`instrupy.util.SphericalGeometry`

            :return: Swath-width in kilometers.
            :rtype: float

        """
        if instru_look_angle_deg != 0:
            raise RuntimeError("Only a 0 deg instrument look angle is supported for the case of conical scan.")

        # calculate the radius of the small-circle on the Earth surface on which the imaged arc lies.
        h = alt_km * 1e3
        Re = Constants.radiusOfEarthInKM * 1e3         
        Rs = Re + h 

        gamma = np.deg2rad(self.offNadirAngle) 
        theta_i = np.arcsin(np.sin(gamma)*Rs/Re)
        alpha = theta_i - gamma # here alpha is the earth centered angle

        # small circle radius
        r = Re * np.sin(alpha)

        arc_len = np.deg2rad(self.clockAngleRange)*r        

        return arc_len*1e-3

class RadiometerModel(Entity):
    """A radiometer class estimating observation data-metrics.      
      
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
    
    :ivar pointing_option: List of ``Orientation`` objects which specify the orientations into which the instrument-axis can be maneuvered. 
                            The orientation must be specified in the NADIR_POINTING frame.
    :vartype pointing_option: list, :class:`orbitpy.util.Orientation` 

    :ivar dataRate: Rate of data recorded (Mega-bits-per-sec) during nominal operations.
    :vartype dataRate: float  

    :ivar bitsPerPixel: Number of bits encoded per pixel of image.
    :vartype bitsPerPixel: int    

    :ivar antenna: Antenna specifications.
    :vartype antenna: :class:`instrupy.util.Antenna`

    :ivar operatingFrequency: Operating radar center frequency in (Hertz).
    :vartype operatingFrequency: float
    
    :ivar system: Radiometer system object.
    :varytype system: :class:`instrupy.radiometer_model.TotalPowerRadiometerSystem` or :class:`instrupy.radiometer_model.UnbalancedDikeRadiometerSystem` or :class:`instrupy.radiometer_model.BalancedDikeRadiometerSystem` or :class:`instrupy.radiometer_model.NoiseAddingRadiometerSystem`
    
    :ivar scan: Scan object.
    :vartype scan: :class:`instrupy.radiometer_model.FixedScan` or :class:`instrupy.radiometer_model.CrossTrackScan` or :class:`instrupy.radiometer_model.ConicalScan`.

    :ivar targetBrightnessTemp: Target brightness temperature in Kelvin.
    :vartype targetBrightnessTemp: float

    :ivar _id: Unique identifier.
    :vartype _id: str or int
          
    """
    def __init__(self, name=None, mass=None, volume=None, power=None,  orientation=None, 
            fieldOfViewGeometry=None, sceneFieldOfViewGeometry=None, maneuver=None, pointingOption=None, 
            dataRate=None, bitsPerPixel=None, antenna=None, operatingFrequency=None, 
            system=None, scan=None, targetBrightnessTemp=None, _id=None):
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
        self.antenna = copy.deepcopy(antenna) if antenna is not None and isinstance(antenna, Antenna) else None
        self.operatingFrequency = float(operatingFrequency) if operatingFrequency is not None else None
        self.system = copy.deepcopy(system) if system is not None and (isinstance(system, TotalPowerRadiometerSystem) or isinstance(system, UnbalancedDikeRadiometerSystem) or isinstance(system, BalancedDikeRadiometerSystem) or isinstance(system, NoiseAddingRadiometerSystem))  else None
        self.scan = copy.deepcopy(scan) if scan is not None and (isinstance(scan, FixedScan) or isinstance(scan, CrossTrackScan) or isinstance(scan, ConicalScan)) else None
        self.targetBrightnessTemp = float(targetBrightnessTemp) if targetBrightnessTemp is not None else None     

        super(RadiometerModel,self).__init__(_id, "Radiometer")
        
    @staticmethod
    def from_dict(d):
        """ Parses an ``RadiometerModel`` object from a normalized JSON dictionary.

        The following default values are assigned to the object instance parameters in case of 
        :class:`None` values or missing key/value pairs in the input dictionary.

        .. csv-table:: Default values
            :header: Parameter, Default Value
            :widths: 10,40

            scanTech, ScanTech.Fixed
            orientation, Orientation.Convention.REF_FRAME_ALIGNED and ReferenceFrame.SC_BODY_FIXED
            sceneFieldOfViewGeometry, (Instrument) fieldOfViewGeometry
            targetBrightnessTemp, 290 Kelvin
            _id, random string
        
        :param d: Normalized JSON dictionary with the corresponding model specifications. 
        :paramtype d: dict

        :returns: ``RadiometerModel object`` initialized with the input specifications.
        :rtype: :class:`instrupy.RadiometerModel`

        """
        # Only side-looking orientation of instrument supported for synthetic aperture radar 
        orien_dict = d.get("orientation", {"referenceFrame": "SC_BODY_FIXED", "convention": "REF_FRAME_ALIGNED"})
        orientation = Orientation.from_dict(orien_dict)

        # parse maneuver
        maneuver = Maneuver.from_json(d.get("maneuver", None))

        # parse antenna object and get the field-of-view geometry
        antenna_dict = d.get("antenna", None)
        if antenna_dict:
            antenna = Antenna.from_dict(antenna_dict)
            fieldOfViewGeometry = antenna.get_spherical_geometry(d.get("operatingFrequency", None))            
        else:
            antenna = None
            fieldOfViewGeometry = None

        sceneFieldOfViewGeometry = SphericalGeometry.from_dict( d.get("sceneFieldOfViewGeometry")) if  d.get("sceneFieldOfViewGeometry") else fieldOfViewGeometry    
        
        # parse the pointing options as a list of Orientation objects.
        pnt_opt_dict = d.get("pointingOption", None)
        _pointing_option = None
        if pnt_opt_dict:
            # translate to a list of Orientation objects
            if isinstance(pnt_opt_dict, list):
                _pointing_option = [Orientation.from_dict(x) for x in pnt_opt_dict]
            else:
                _pointing_option = [Orientation.from_dict(pnt_opt_dict)]

        # parse the system and the systemType objects
        system_lookup = { 'TOTAL_POWER' : TotalPowerRadiometerSystem,  'UNBALANCED_DICKE': UnbalancedDikeRadiometerSystem, 
                          'BALANCED_DICKE': BalancedDikeRadiometerSystem, 'NOISE_ADDING': NoiseAddingRadiometerSystem}
        system_dict = d.get("system", None)
        if system_dict:
            systemType = SystemType.get(system_dict["@type"])
            system = system_lookup[systemType].from_dict(system_dict)
        else:
            systemType = None
            system = None

        # parse the scan and the scanTechnique objects
        scan_lookup = { 'FIXED':FixedScan, 'CROSS_TRACK':CrossTrackScan, 'CONICAL':ConicalScan}
        scan_dict = d.get("scan", None)
        if scan_dict:
            scanTechnique = ScanTech.get(scan_dict["@type"])
            scan = scan_lookup[scanTechnique].from_dict(scan_dict)
        else:
            scanTechnique = None
            scan = None
        # TODO: check for compatibility of the scan technique with the instrument orientation, maneuver specifications

        return RadiometerModel(
                        name = d.get("name", None),
                        mass = d.get("mass", None),
                        volume = d.get("volume", None),
                        power = d.get("power", None),
                        orientation = orientation,
                        fieldOfViewGeometry =  fieldOfViewGeometry,
                        sceneFieldOfViewGeometry = sceneFieldOfViewGeometry,
                        maneuver =  maneuver,
                        pointingOption = _pointing_option,
                        dataRate = d.get("dataRate", None),
                        bitsPerPixel = d.get("bitsPerPixel", None),
                        antenna = antenna,
                        operatingFrequency = d.get("operatingFrequency", None),
                        system = system,
                        scan = scan,
                        targetBrightnessTemp = d.get("targetBrightnessTemp", 290),
                        _id = str(d.get("@id", uuid.uuid4()))
                        )

    def to_dict(self):
        """ Translate the ``RadiometerModel`` object to a Python dictionary such that it can be uniquely reconstructed back from the dictionary.

        :returns: ``RadiometerModel`` specifications as python dictionary.
        :rtype: dict

        """
        fieldOfViewGeometry_dict = self.fieldOfView.sph_geom.to_dict() if self.fieldOfView is not None and isinstance(self.fieldOfView, ViewGeometry) else None
        sceneFieldOfViewGeometry_dict = self.sceneFieldOfView.sph_geom.to_dict() if self.sceneFieldOfView is not None and isinstance(self.sceneFieldOfView, ViewGeometry) else None
        orientation_dict = self.orientation.to_dict() if self.orientation is not None and isinstance(self.orientation, Orientation) else None
        maneuver_dict = self.maneuver.to_dict() if self.maneuver is not None and isinstance(self.maneuver, Maneuver) else None
        pointing_opt_dict = [Orientation.to_dict(x) for x in self.pointingOption] if self.pointingOption is not None else None
        antenna_dict = self.antenna.to_dict() if self.antenna is not None and isinstance(self.antenna, Antenna) else None
        system_dict = self.system.to_dict() if self.system is not None and  (isinstance(self.system, TotalPowerRadiometerSystem) or isinstance(self.system, UnbalancedDikeRadiometerSystem) or isinstance(self.system, BalancedDikeRadiometerSystem) or isinstance(self.system, NoiseAddingRadiometerSystem)) else None
        scan_dict = self.scan.to_dict() if self.scan is not None and (isinstance(self.scan, FixedScan) or isinstance(self.scan, CrossTrackScan) or isinstance(self.scan, ConicalScan)) else None
        return dict({
                "@type": "Radiometer",
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
                "antenna": antenna_dict,
                "operatingFrequency": self.operatingFrequency,
                "system": system_dict,
                "scan": scan_dict,
                "targetBrightnessTemp": self.targetBrightnessTemp,
                "@id": self._id
                })

    def __repr__(self):
        return "RadiometerModel.from_dict({})".format(self.to_dict())

    def calc_data_metrics(self, alt_km, gnd_spd):
        """ Calculate typical observation data metrics. This function is invoked by the function ``Instrument.calc_data_metrics(.)`` class in the ``base`` module.

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
                
                    * :code:`radiometric res [K]` (:class:`float`) Radiometric resolution/ sensitivity in Kelvin.
                    * :code:`ground pixel along-track resolution [m]` (:class:`float`) Along-track resolution (meters) of an ground-pixel centered about observation point.
                    * :code:`ground pixel cross-track resolution [m]` (:class:`float`) Cross-track resolution (meters) of an ground-pixel centered about observation point.
                    * :code:`swath-width [m]` (:class:`float`) Swath-width (meters) of the strip of which the imaged pixel is part off.
                    * :code:`beam efficiency` (:class:`float`) Beam efficiency.
                    * :code:`incidence angle [deg]` (:class:`float`) Observation incidence angle (degrees) at the ground-pixel.

                .. todo:: The along-track and cross-track pixel resolutions are accurate only for pixels imaged at strictly sidelooking geometry (roll-only, no pitch). Needs revision.

        :rtype: dict 

        """
        range_km = alt_km

        # Calculate look angle to the target location
        look_angle = 0.0
        
        # Look angle to corresponding incidence angle conversion for spherical Earth (incidence angle at the target location)
        incidence_angle = 0.0

        ############## Calculate the pixel resolution. ##############
        # The size of the antenna footprint at the target-location corresponds to the pixel dimensions. It is assumed that
        # in case of rectangular antennas the along-track resolution corresponds to the antenna "height" dimension while the cross-track resolution corresponds to
        # the antenna "width" dimension.
        if self.fieldOfView.sph_geom.shape == SphericalGeometry.Shape.RECTANGULAR:
            iFOV_AT_deg = self.fieldOfView.sph_geom.angle_height
            iFOV_CT_deg = self.fieldOfView.sph_geom.angle_width
        elif self.fieldOfView.sph_geom.shape == SphericalGeometry.Shape.CIRCULAR:
            iFOV_AT_deg = self.fieldOfView.sph_geom.diameter
            iFOV_CT_deg = iFOV_AT_deg
        else:
            raise NotImplementedError

        print('iFOV_AT_deg', iFOV_AT_deg)
        print('iFOV_CT_deg', iFOV_CT_deg)
        res_AT_m = np.deg2rad(iFOV_AT_deg)*range_km*1.0e3
        res_CT_m = np.deg2rad(iFOV_CT_deg)*range_km*1.0e3/np.cos(incidence_angle)
        
        ############## calculate the radiometric resolution ##############
        # calculate the dwell time per ground-pixel
        sat_speed_mps = gnd_spd
        td = self.scan.compute_dwell_time_per_ground_pixel(res_AT_m=res_AT_m, sat_speed_kmps=sat_speed_mps*1e-3, iFOV_CT_deg=iFOV_CT_deg)

        print('td', td)
        rad_res = self.system.compute_radiometric_resolution(td, self.antenna, self.targetBrightnessTemp)

        ############## calculate beam-efficiency ##############
        be = 0.0

        ############## Calculate the swath-width. ##############
        # The swath-width is calculated based on the instrument look angle, which may or may-not be equal to the target look angle. 
        # If the instrument look-angle = target look angle, then it implies that the target is at the middle of the swath.
        instru_look_angle = 0.0
                
        obsv_metrics = {}
        obsv_metrics["fov"] = iFOV_CT_deg if iFOV_CT_deg is not None else np.nan
        obsv_metrics["ground pixel along-track resolution [m]"] = round(res_AT_m, 2) if res_AT_m is not None else np.nan
        obsv_metrics["ground pixel cross-track resolution [m]"] = round(res_CT_m, 2) if res_CT_m is not None else np.nan
        obsv_metrics["sensitivity [K]"] = round(rad_res, 3) if rad_res is not None else np.nan
        obsv_metrics["incidence angle [deg]"] = round(np.rad2deg(incidence_angle), 2) if incidence_angle is not None else np.nan
        obsv_metrics["beam efficiency"] = round(be, 2) if be is not np.nan else np.nan
        
        return obsv_metrics

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
