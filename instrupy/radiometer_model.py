""" 
.. module:: radiometer_model

:synopsis: *Module to handle radiometer.*

        References: [1] Chapter 6,7 in "Microwave Radar and Radiometric Remote Sensing," David Gardner Long , Fawwaz T. Ulaby 2014 

"""
import json
import copy
from typing import Coroutine
import uuid
import numpy as np
import warnings
from instrupy.util import Entity, EnumEntity, Orientation, SphericalGeometry, ViewGeometry, Maneuver, Antenna, GeoUtilityFunctions, MathUtilityFunctions, Constants, FileUtilityFunctions

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
    """ Container to handle total power radiometer system. Refer Section 7.4, 7.5 in [1].

    :ivar tlLoss: Transmission line loss in decibels.
    :vartype tlLoss: float

    :ivar tlPhyTemp: Transmission line physical temperature in Kelvin.
    :vartype tlPhyTemp: float

    :ivar rfAmpGain: RF amplifier gain in decibels.
    :vartype rfAmpGain: float

    :ivar rfAmpInpNoiseTemp: RF amplifier input noise temperature in Kelvin.
    :vartype rfAmpInpNoiseTemp: float

    :ivar rfAmpGainVariation: RF amplifier gain variation.
    :vartype rfAmpGainVariation: float

    :ivar mixerInpNoiseAmp: Mixer input noise temperature in Kelvin.
    :vartype mixerInpNoiseAmp: float

    :ivar ifAmpGain: Intermediate frequency amplifier gain in decibels.
    :vartype ifAmpGain: float

    :ivar ifAmpInpNoiseTemp: Intermediate frequency amplifier input noise temperature in Kelvin.
    :vartype ifAmpInpNoiseTemp: float

    :ivar ifAmpGainVariation: IF amplifier gain variation.
    :vartype ifAmpGainVariation: float

    :ivar integratorVoltageGain: Integrator voltage gain (unitless).
    :vartype integratorVoltageGain: float

    :ivar predetectionGain: Pre-detection stage gain in decibels.
    :vartype predetectionGain: float

    :ivar predetectionInpNoiseTemp: Pre-detection input noise temperature in Kelvin.
    :vartype predetectionInpNoiseTemp: float

    :ivar predetectionGainVariation: Pre-detection stage gain variation.
    :vartype predetectionGainVariation: float

    :ivar _id: Unique identifier.
    :vartype _id: str or int

    """
    def __init__(self, tlLoss=None, tlPhyTemp=None, rfAmpGain=None, rfAmpInputNoiseTemp=None, rfAmpGainVariation=None, 
                 mixerInputNoiseAmp=None, ifAmpGain=None, ifAmpInputNoiseTemp=None, ifAmpGainVariation=None, 
                 integratorVoltageGain=None, predetectionGain=None, predetectionInpNoiseTemp=None, predetectionGainVariation=None,
                 _id=None):
        
        self.tlLoss = float(tlLoss) if tlLoss is not None else None
        self.tlPhyTemp = float(tlPhyTemp) if tlPhyTemp is not None else None
        self.rfAmpGain = float(rfAmpGain) if rfAmpGain is not None else None
        self.rfAmpInputNoiseTemp = float(rfAmpInputNoiseTemp) if rfAmpInputNoiseTemp is not None else None
        self.rfAmpGainVariation = float(rfAmpGainVariation) if rfAmpGainVariation is not None else None
        self.mixerInputNoiseAmp = float(mixerInputNoiseAmp) if mixerInputNoiseAmp is not None else None
        self.ifAmpGain = float(ifAmpGain) if ifAmpGain is not None else None
        self.ifAmpInputNoiseTemp = float(ifAmpInputNoiseTemp) if ifAmpInputNoiseTemp is not None else None
        self.ifAmpGainVariation = float(ifAmpGainVariation) if ifAmpGainVariation is not None else None
        self.integratorVoltageGain = float(integratorVoltageGain) if integratorVoltageGain is not None else None
        self.predetectionGain = float(predetectionGain) if predetectionGain is not None else None
        self.predetectionInpNoiseTemp = float(predetectionInpNoiseTemp) if predetectionInpNoiseTemp is not None else None
        self.predetectionGainVariation = float(predetectionGainVariation) if predetectionGainVariation is not None else None

        super(TotalPowerRadiometerSystem, self).__init__(_id, "TotalPowerRadiometerSystem")
    
    @staticmethod
    def from_dict(d):
        """Parses an TotalPowerRadiometerSystem object from a normalized JSON dictionary.
        
        :param d: Dictionary with the total-power radiometer system specifications.
        :paramtype d: dict

        :return: TotalPowerRadiometerSystem object.
        :rtype: :class:`instrupy.radiometer_model.TotalPowerRadiometerSystem`

        """             
        return TotalPowerRadiometerSystem(
                tlLoss = d.get("tlLoss", None),
                tlPhyTemp = d.get("tlPhyTemp", None),
                rfAmpGain = d.get("rfAmpGain", None),
                rfAmpInputNoiseTemp = d.get("rfAmpInputNoiseTemp", None),
                rfAmpGainVariation = d.get("rfAmpGainVariation", None),
                mixerInputNoiseAmp = d.get("mixerInputNoiseAmp", None),
                ifAmpGain = d.get("ifAmpGain", None),
                ifAmpInputNoiseTemp = d.get("ifAmpInputNoiseTemp", None),
                ifAmpGainVariation = d.get("ifAmpGainVariation", None),
                integratorVoltageGain = d.get("integratorVoltageGain", None),
                predetectionGain = d.get("predetectionGain", None),
                predetectionInpNoiseTemp = d.get("predetectionInpNoiseTemp", None),
                predetectionGainVariation = d.get("predetectionGainVariation", None),
                _id = d.get("@id", None)
                )
    
    def to_dict(self):
        """ Translate the TotalPowerRadiometerSystem object to a Python dictionary such that it can be uniquely reconstructed back from the dictionary.
        
        :return: TotalPowerRadiometerSystem object as python dictionary
        :rtype: dict

        """
        return dict({"tlLoss": self.tlLoss,
                     "tlPhyTemp": self.tlPhyTemp,
                     "rfAmpGain": self.rfAmpGain,
                     "rfAmpInputNoiseTemp": self.rfAmpInputNoiseTemp, 
                     "rfAmpGainVariation": self.rfAmpGainVariation,
                     "mixerInputNoiseAmp": self.mixerInputNoiseAmp,
                     "ifAmpGain": self.ifAmpGain,
                     "ifAmpInputNoiseTemp": self.ifAmpInputNoiseTemp,
                     "ifAmpGainVariation": self.ifAmpGainVariation,
                     "integratorVoltageGain": self.integratorVoltageGain,
                     "predetectionGain": self.predetectionGain,
                     "predetectionInpNoiseTemp": self.predetectionInpNoiseTemp,
                     "predetectionGainVariation": self.predetectionGainVariation,
                     "@id": self._id
                    })
    
    def __repr__(self):
        return "TotalPowerRadiometerSystem.from_dict({})".format(self.to_dict())

    def __eq__(self, other):
        # Equality test is simple one which compares the data attributes.
        # note that _id data attribute may be different
        if(isinstance(self, other.__class__)):
            return (self.tlLoss==other.tlLoss) and (self.tlPhyTemp==other.tlPhyTemp) and (self.rfAmpGain==other.rfAmpGain) and \
                    (self.rfAmpInputNoiseTemp==other.rfAmpInputNoiseTemp) and (self.rfAmpGainVariation==other.rfAmpGainVariation) and \
                    (self.mixerInputNoiseAmp==other.mixerInputNoiseAmp) and (self.ifAmpGain==other.ifAmpGain) and (self.ifAmpInputNoiseTemp==other.ifAmpInputNoiseTemp) and \
                    (self.ifAmpGainVariation==other.ifAmpGainVariation) and (self.integratorVoltageGain==other.integratorVoltageGain) and (self.predetectionGain==other.predetectionGain) and \
                    (self.predetectionInpNoiseTemp==other.predetectionInpNoiseTemp) and (self.predetectionGainVariation==other.predetectionGainVariation)                   
        else:
            return NotImplemented

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

    :ivar rfAmpGainVariation: RF amplifier gain variation.
    :vartype rfAmpGainVariation: float

    :ivar mixerInpNoiseAmp: Mixer input noise temperature in Kelvin.
    :vartype mixerInpNoiseAmp: float

    :ivar ifAmpGain: Intermediate frequency amplifier gain in decibels.
    :vartype ifAmpGain: float

    :ivar ifAmpInpNoiseTemp: Intermediate frequency amplifier input noise temperature in Kelvin.
    :vartype ifAmpInpNoiseTemp: float

    :ivar ifAmpGainVariation: IF amplifier gain variation.
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

    :ivar predetectionGainVariation: Pre-detection stage gain variation.
    :vartype predetectionGainVariation: float

    :ivar _id: Unique identifier.
    :vartype _id: str or int

    """
    def __init__(self, tlLoss=None, tlPhyTemp=None, rfAmpGain=None, rfAmpInputNoiseTemp=None, rfAmpGainVariation=None, 
                 mixerInputNoiseAmp=None, ifAmpGain=None, ifAmpInputNoiseTemp=None, ifAmpGainVariation=None, 
                 dickeSwitchOutputNoiseTemperature=None, referenceTemperature=None,
                 integratorVoltageGain=None, predetectionGain=None, predetectionInpNoiseTemp=None, predetectionGainVariation=None,
                 _id=None):
        
        self.tlLoss = float(tlLoss) if tlLoss is not None else None
        self.tlPhyTemp = float(tlPhyTemp) if tlPhyTemp is not None else None
        self.rfAmpGain = float(rfAmpGain) if rfAmpGain is not None else None
        self.rfAmpInputNoiseTemp = float(rfAmpInputNoiseTemp) if rfAmpInputNoiseTemp is not None else None
        self.rfAmpGainVariation = float(rfAmpGainVariation) if rfAmpGainVariation is not None else None
        self.mixerInputNoiseAmp = float(mixerInputNoiseAmp) if mixerInputNoiseAmp is not None else None
        self.ifAmpGain = float(ifAmpGain) if ifAmpGain is not None else None
        self.ifAmpInputNoiseTemp = float(ifAmpInputNoiseTemp) if ifAmpInputNoiseTemp is not None else None
        self.ifAmpGainVariation = float(ifAmpGainVariation) if ifAmpGainVariation is not None else None
        self.dickeSwitchOutputNoiseTemperature = float(dickeSwitchOutputNoiseTemperature) if dickeSwitchOutputNoiseTemperature is not None else None
        self.referenceTemperature = float(referenceTemperature) if referenceTemperature is not None else None
        self.integratorVoltageGain = float(integratorVoltageGain) if integratorVoltageGain is not None else None
        self.predetectionGain = float(predetectionGain) if predetectionGain is not None else None
        self.predetectionInpNoiseTemp = float(predetectionInpNoiseTemp) if predetectionInpNoiseTemp is not None else None
        self.predetectionGainVariation = float(predetectionGainVariation) if predetectionGainVariation is not None else None

        super(UnbalancedDikeRadiometerSystem, self).__init__(_id, "UnbalancedDikeRadiometerSystem")
    
    @staticmethod
    def from_dict(d):
        """Parses an UnbalancedDikeRadiometerSystem object from a normalized JSON dictionary.
        
        :param d: Dictionary with the unbalanced Dicke radiometer system specifications.
        :paramtype d: dict

        :return: UnbalancedDikeRadiometerSystem object.
        :rtype: :class:`instrupy.radiometer_model.UnbalancedDikeRadiometerSystem`

        """             
        return UnbalancedDikeRadiometerSystem(
                tlLoss = d.get("tlLoss", None),
                tlPhyTemp = d.get("tlPhyTemp", None),
                rfAmpGain = d.get("rfAmpGain", None),
                rfAmpInputNoiseTemp = d.get("rfAmpInputNoiseTemp", None),
                rfAmpGainVariation = d.get("rfAmpGainVariation", None),
                mixerInputNoiseAmp = d.get("mixerInputNoiseAmp", None),
                ifAmpGain = d.get("ifAmpGain", None),
                ifAmpInputNoiseTemp = d.get("ifAmpInputNoiseTemp", None),
                ifAmpGainVariation = d.get("ifAmpGainVariation", None),
                dickeSwitchOutputNoiseTemperature = d.get("dickeSwitchOutputNoiseTemperature", None),
                referenceTemperature = d.get("referenceTemperature", None),
                integratorVoltageGain = d.get("integratorVoltageGain", None),
                predetectionGain = d.get("predetectionGain", None),
                predetectionInpNoiseTemp = d.get("predetectionInpNoiseTemp", None),
                predetectionGainVariation = d.get("predetectionGainVariation", None),
                _id = d.get("@id", None)
                )
    
    def to_dict(self):
        """ Translate the UnbalancedDikeRadiometerSystem object to a Python dictionary such that it can be uniquely reconstructed back from the dictionary.
        
        :return: UnbalancedDikeRadiometerSystem object as python dictionary
        :rtype: dict

        """
        return dict({"tlLoss": self.tlLoss,
                     "tlPhyTemp": self.tlPhyTemp,
                     "rfAmpGain": self.rfAmpGain,
                     "rfAmpInputNoiseTemp": self.rfAmpInputNoiseTemp, 
                     "rfAmpGainVariation": self.rfAmpGainVariation,
                     "mixerInputNoiseAmp": self.mixerInputNoiseAmp,
                     "ifAmpGain": self.ifAmpGain,
                     "ifAmpInputNoiseTemp": self.ifAmpInputNoiseTemp,
                     "ifAmpGainVariation": self.ifAmpGainVariation,
                     "dickeSwitchOutputNoiseTemperature":self.dickeSwitchOutputNoiseTemperature,
                     "referenceTemperature": self.referenceTemperature,
                     "integratorVoltageGain": self.integratorVoltageGain,
                     "predetectionGain": self.predetectionGain,
                     "predetectionInpNoiseTemp": self.predetectionInpNoiseTemp,
                     "predetectionGainVariation": self.predetectionGainVariation,
                     "@id": self._id
                    })
    
    def __repr__(self):
        return "UnbalancedDikeRadiometerSystem.from_dict({})".format(self.to_dict())

    def __eq__(self, other):
        # Equality test is simple one which compares the data attributes.
        # note that _id data attribute may be different
        if(isinstance(self, other.__class__)):
            return (self.tlLoss==other.tlLoss) and (self.tlPhyTemp==other.tlPhyTemp) and (self.rfAmpGain==other.rfAmpGain) and \
                    (self.rfAmpInputNoiseTemp==other.rfAmpInputNoiseTemp) and (self.rfAmpGainVariation==other.rfAmpGainVariation) and \
                    (self.mixerInputNoiseAmp==other.mixerInputNoiseAmp) and (self.ifAmpGain==other.ifAmpGain) and (self.ifAmpInputNoiseTemp==other.ifAmpInputNoiseTemp) and \
                    (self.ifAmpGainVariation==other.ifAmpGainVariation) and (self.dickeSwitchOutputNoiseTemperature==other.dickeSwitchOutputNoiseTemperature) and (self.referenceTemperature==other.referenceTemperature) and \
                    (self.integratorVoltageGain==other.integratorVoltageGain) and (self.predetectionGain==other.predetectionGain) and \
                    (self.predetectionInpNoiseTemp==other.predetectionInpNoiseTemp) and (self.predetectionGainVariation==other.predetectionGainVariation)                   
        else:
            return NotImplemented

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

    :ivar rfAmpGainVariation: RF amplifier gain variation.
    :vartype rfAmpGainVariation: float

    :ivar mixerInpNoiseAmp: Mixer input noise temperature in Kelvin.
    :vartype mixerInpNoiseAmp: float

    :ivar ifAmpGain: Intermediate frequency amplifier gain in decibels.
    :vartype ifAmpGain: float

    :ivar ifAmpInpNoiseTemp: Intermediate frequency amplifier input noise temperature in Kelvin.
    :vartype ifAmpInpNoiseTemp: float

    :ivar ifAmpGainVariation: IF amplifier gain variation.
    :vartype ifAmpGainVariation: float

    :ivar dickeSwitchOutputNoiseTemperature: Dicke switch noise temperature in Kelvin referenced to the output port.
    :vartype dickeSwitchOutputNoiseTemperature: float

    :ivar integratorVoltageGain: Integrator voltage gain (unitless).
    :vartype integratorVoltageGain: float

    :ivar predetectionGain: Pre-detection stage gain in decibels.
    :vartype predetectionGain: float

    :ivar predetectionInpNoiseTemp: Pre-detection input noise temperature in Kelvin.
    :vartype predetectionInpNoiseTemp: float

    :ivar predetectionGainVariation: Pre-detection stage gain variation.
    :vartype predetectionGainVariation: float

    :ivar _id: Unique identifier.
    :vartype _id: str or int

    """
    def __init__(self, tlLoss=None, tlPhyTemp=None, rfAmpGain=None, rfAmpInputNoiseTemp=None, rfAmpGainVariation=None, 
                 mixerInputNoiseAmp=None, ifAmpGain=None, ifAmpInputNoiseTemp=None, ifAmpGainVariation=None, 
                 dickeSwitchOutputNoiseTemperature=None,
                 integratorVoltageGain=None, predetectionGain=None, predetectionInpNoiseTemp=None, predetectionGainVariation=None,
                 _id=None):
        
        self.tlLoss = float(tlLoss) if tlLoss is not None else None
        self.tlPhyTemp = float(tlPhyTemp) if tlPhyTemp is not None else None
        self.rfAmpGain = float(rfAmpGain) if rfAmpGain is not None else None
        self.rfAmpInputNoiseTemp = float(rfAmpInputNoiseTemp) if rfAmpInputNoiseTemp is not None else None
        self.rfAmpGainVariation = float(rfAmpGainVariation) if rfAmpGainVariation is not None else None
        self.mixerInputNoiseAmp = float(mixerInputNoiseAmp) if mixerInputNoiseAmp is not None else None
        self.ifAmpGain = float(ifAmpGain) if ifAmpGain is not None else None
        self.ifAmpInputNoiseTemp = float(ifAmpInputNoiseTemp) if ifAmpInputNoiseTemp is not None else None
        self.ifAmpGainVariation = float(ifAmpGainVariation) if ifAmpGainVariation is not None else None
        self.dickeSwitchOutputNoiseTemperature = float(dickeSwitchOutputNoiseTemperature) if dickeSwitchOutputNoiseTemperature is not None else None
        self.integratorVoltageGain = float(integratorVoltageGain) if integratorVoltageGain is not None else None
        self.predetectionGain = float(predetectionGain) if predetectionGain is not None else None
        self.predetectionInpNoiseTemp = float(predetectionInpNoiseTemp) if predetectionInpNoiseTemp is not None else None
        self.predetectionGainVariation = float(predetectionGainVariation) if predetectionGainVariation is not None else None

        super(BalancedDikeRadiometerSystem, self).__init__(_id, "BalancedDikeRadiometerSystem")
    
    @staticmethod
    def from_dict(d):
        """Parses an BalancedDikeRadiometerSystem object from a normalized JSON dictionary.
        
        :param d: Dictionary with the balanced Dicke radiometer system specifications.
        :paramtype d: dict

        :return: BalancedDikeRadiometerSystem object.
        :rtype: :class:`instrupy.radiometer_model.BalancedDikeRadiometerSystem`

        """             
        return BalancedDikeRadiometerSystem(
                tlLoss = d.get("tlLoss", None),
                tlPhyTemp = d.get("tlPhyTemp", None),
                rfAmpGain = d.get("rfAmpGain", None),
                rfAmpInputNoiseTemp = d.get("rfAmpInputNoiseTemp", None),
                rfAmpGainVariation = d.get("rfAmpGainVariation", None),
                mixerInputNoiseAmp = d.get("mixerInputNoiseAmp", None),
                ifAmpGain = d.get("ifAmpGain", None),
                ifAmpInputNoiseTemp = d.get("ifAmpInputNoiseTemp", None),
                ifAmpGainVariation = d.get("ifAmpGainVariation", None),
                dickeSwitchOutputNoiseTemperature = d.get("dickeSwitchOutputNoiseTemperature", None),
                integratorVoltageGain = d.get("integratorVoltageGain", None),
                predetectionGain = d.get("predetectionGain", None),
                predetectionInpNoiseTemp = d.get("predetectionInpNoiseTemp", None),
                predetectionGainVariation = d.get("predetectionGainVariation", None),
                _id = d.get("@id", None)
                )
    
    def to_dict(self):
        """ Translate the BalancedDikeRadiometerSystem object to a Python dictionary such that it can be uniquely reconstructed back from the dictionary.
        
        :return: BalancedDikeRadiometerSystem object as python dictionary
        :rtype: dict

        """
        return dict({"tlLoss": self.tlLoss,
                     "tlPhyTemp": self.tlPhyTemp,
                     "rfAmpGain": self.rfAmpGain,
                     "rfAmpInputNoiseTemp": self.rfAmpInputNoiseTemp, 
                     "rfAmpGainVariation": self.rfAmpGainVariation,
                     "mixerInputNoiseAmp": self.mixerInputNoiseAmp,
                     "ifAmpGain": self.ifAmpGain,
                     "ifAmpInputNoiseTemp": self.ifAmpInputNoiseTemp,
                     "ifAmpGainVariation": self.ifAmpGainVariation,
                     "dickeSwitchOutputNoiseTemperature":self.dickeSwitchOutputNoiseTemperature,
                     "integratorVoltageGain": self.integratorVoltageGain,
                     "predetectionGain": self.predetectionGain,
                     "predetectionInpNoiseTemp": self.predetectionInpNoiseTemp,
                     "predetectionGainVariation": self.predetectionGainVariation,
                     "@id": self._id
                    })
    
    def __repr__(self):
        return "BalancedDikeRadiometerSystem.from_dict({})".format(self.to_dict())

    def __eq__(self, other):
        # Equality test is simple one which compares the data attributes.
        # note that _id data attribute may be different
        if(isinstance(self, other.__class__)):
            return (self.tlLoss==other.tlLoss) and (self.tlPhyTemp==other.tlPhyTemp) and (self.rfAmpGain==other.rfAmpGain) and \
                    (self.rfAmpInputNoiseTemp==other.rfAmpInputNoiseTemp) and (self.rfAmpGainVariation==other.rfAmpGainVariation) and \
                    (self.mixerInputNoiseAmp==other.mixerInputNoiseAmp) and (self.ifAmpGain==other.ifAmpGain) and (self.ifAmpInputNoiseTemp==other.ifAmpInputNoiseTemp) and \
                    (self.ifAmpGainVariation==other.ifAmpGainVariation) and (self.dickeSwitchOutputNoiseTemperature==other.dickeSwitchOutputNoiseTemperature) and \
                    (self.integratorVoltageGain==other.integratorVoltageGain) and (self.predetectionGain==other.predetectionGain) and \
                    (self.predetectionInpNoiseTemp==other.predetectionInpNoiseTemp) and (self.predetectionGainVariation==other.predetectionGainVariation)                   
        else:
            return NotImplemented

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

    :ivar rfAmpGainVariation: RF amplifier gain variation.
    :vartype rfAmpGainVariation: float

    :ivar mixerInpNoiseAmp: Mixer input noise temperature in Kelvin.
    :vartype mixerInpNoiseAmp: float

    :ivar ifAmpGain: Intermediate frequency amplifier gain in decibels.
    :vartype ifAmpGain: float

    :ivar ifAmpInpNoiseTemp: Intermediate frequency amplifier input noise temperature in Kelvin.
    :vartype ifAmpInpNoiseTemp: float

    :ivar ifAmpGainVariation: IF amplifier gain variation.
    :vartype ifAmpGainVariation: float

    :ivar excessNoiseTemperature: Excess noise temperature (added noise to the receiver input during the diode ON half-cycle) in Kelvin referenced to the output port.
    :vartype excessNoiseTemperature: float

    :ivar integratorVoltageGain: Integrator voltage gain (unitless).
    :vartype integratorVoltageGain: float

    :ivar predetectionGain: Pre-detection stage gain in decibels.
    :vartype predetectionGain: float

    :ivar predetectionInpNoiseTemp: Pre-detection input noise temperature in Kelvin.
    :vartype predetectionInpNoiseTemp: float

    :ivar predetectionGainVariation: Pre-detection stage gain variation.
    :vartype predetectionGainVariation: float

    :ivar _id: Unique identifier.
    :vartype _id: str or int

    """
    def __init__(self, tlLoss=None, tlPhyTemp=None, rfAmpGain=None, rfAmpInputNoiseTemp=None, rfAmpGainVariation=None, 
                 mixerInputNoiseAmp=None, ifAmpGain=None, ifAmpInputNoiseTemp=None, ifAmpGainVariation=None, 
                 excessNoiseTemperature=None,
                 integratorVoltageGain=None, predetectionGain=None, predetectionInpNoiseTemp=None, predetectionGainVariation=None,
                 _id=None):
        
        self.tlLoss = float(tlLoss) if tlLoss is not None else None
        self.tlPhyTemp = float(tlPhyTemp) if tlPhyTemp is not None else None
        self.rfAmpGain = float(rfAmpGain) if rfAmpGain is not None else None
        self.rfAmpInputNoiseTemp = float(rfAmpInputNoiseTemp) if rfAmpInputNoiseTemp is not None else None
        self.rfAmpGainVariation = float(rfAmpGainVariation) if rfAmpGainVariation is not None else None
        self.mixerInputNoiseAmp = float(mixerInputNoiseAmp) if mixerInputNoiseAmp is not None else None
        self.ifAmpGain = float(ifAmpGain) if ifAmpGain is not None else None
        self.ifAmpInputNoiseTemp = float(ifAmpInputNoiseTemp) if ifAmpInputNoiseTemp is not None else None
        self.ifAmpGainVariation = float(ifAmpGainVariation) if ifAmpGainVariation is not None else None
        self.excessNoiseTemperature = float(excessNoiseTemperature) if excessNoiseTemperature is not None else None
        self.integratorVoltageGain = float(integratorVoltageGain) if integratorVoltageGain is not None else None
        self.predetectionGain = float(predetectionGain) if predetectionGain is not None else None
        self.predetectionInpNoiseTemp = float(predetectionInpNoiseTemp) if predetectionInpNoiseTemp is not None else None
        self.predetectionGainVariation = float(predetectionGainVariation) if predetectionGainVariation is not None else None

        super(NoiseAddingRadiometerSystem, self).__init__(_id, "NoiseAddingRadiometerSystem")
    
    @staticmethod
    def from_dict(d):
        """Parses an NoiseAddingRadiometerSystem object from a normalized JSON dictionary.
        
        :param d: Dictionary with the noise-adding radiometer system specifications.
        :paramtype d: dict

        :return: NoiseAddingRadiometerSystem object.
        :rtype: :class:`instrupy.util.NoiseAddingRadiometerSystem`

        """             
        return NoiseAddingRadiometerSystem(
                tlLoss = d.get("tlLoss", None),
                tlPhyTemp = d.get("tlPhyTemp", None),
                rfAmpGain = d.get("rfAmpGain", None),
                rfAmpInputNoiseTemp = d.get("rfAmpInputNoiseTemp", None),
                rfAmpGainVariation = d.get("rfAmpGainVariation", None),
                mixerInputNoiseAmp = d.get("mixerInputNoiseAmp", None),
                ifAmpGain = d.get("ifAmpGain", None),
                ifAmpInputNoiseTemp = d.get("ifAmpInputNoiseTemp", None),
                ifAmpGainVariation = d.get("ifAmpGainVariation", None),
                excessNoiseTemperature = d.get("excessNoiseTemperature", None),
                integratorVoltageGain = d.get("integratorVoltageGain", None),
                predetectionGain = d.get("predetectionGain", None),
                predetectionInpNoiseTemp = d.get("predetectionInpNoiseTemp", None),
                predetectionGainVariation = d.get("predetectionGainVariation", None),
                _id = d.get("@id", None)
                )
    
    def to_dict(self):
        """ Translate the NoiseAddingRadiometerSystem object to a Python dictionary such that it can be uniquely reconstructed back from the dictionary.
        
        :return: NoiseAddingRadiometerSystem object as python dictionary
        :rtype: dict

        """
        return dict({"tlLoss": self.tlLoss,
                     "tlPhyTemp": self.tlPhyTemp,
                     "rfAmpGain": self.rfAmpGain,
                     "rfAmpInputNoiseTemp": self.rfAmpInputNoiseTemp, 
                     "rfAmpGainVariation": self.rfAmpGainVariation,
                     "mixerInputNoiseAmp": self.mixerInputNoiseAmp,
                     "ifAmpGain": self.ifAmpGain,
                     "ifAmpInputNoiseTemp": self.ifAmpInputNoiseTemp,
                     "ifAmpGainVariation": self.ifAmpGainVariation,
                     "excessNoiseTemperature":self.excessNoiseTemperature,
                     "integratorVoltageGain": self.integratorVoltageGain,
                     "predetectionGain": self.predetectionGain,
                     "predetectionInpNoiseTemp": self.predetectionInpNoiseTemp,
                     "predetectionGainVariation": self.predetectionGainVariation,
                     "@id": self._id
                    })
    
    def __repr__(self):
        return "NoiseAddingRadiometerSystem.from_dict({})".format(self.to_dict())

    def __eq__(self, other):
        # Equality test is simple one which compares the data attributes.
        # note that _id data attribute may be different
        if(isinstance(self, other.__class__)):
            return (self.tlLoss==other.tlLoss) and (self.tlPhyTemp==other.tlPhyTemp) and (self.rfAmpGain==other.rfAmpGain) and \
                    (self.rfAmpInputNoiseTemp==other.rfAmpInputNoiseTemp) and (self.rfAmpGainVariation==other.rfAmpGainVariation) and \
                    (self.mixerInputNoiseAmp==other.mixerInputNoiseAmp) and (self.ifAmpGain==other.ifAmpGain) and (self.ifAmpInputNoiseTemp==other.ifAmpInputNoiseTemp) and \
                    (self.ifAmpGainVariation==other.ifAmpGainVariation) and (self.excessNoiseTemperature==other.excessNoiseTemperature) and \
                    (self.integratorVoltageGain==other.integratorVoltageGain) and (self.predetectionGain==other.predetectionGain) and \
                    (self.predetectionInpNoiseTemp==other.predetectionInpNoiseTemp) and (self.predetectionGainVariation==other.predetectionGainVariation)                   
        else:
            return NotImplemented

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
        """Parses an FixedScan object from a normalized JSON dictionary.
        
        :param d: Dictionary with the fixed-scan specifications.
        :paramtype d: dict

        :return: FixedScan object.
        :rtype: :class:`instrupy.radiometer_model.FixedScan`

        """             
        return FixedScan(_id = d.get("@id", None))
    
    def to_dict(self):
        """ Translate the FixedScan object to a Python dictionary such that it can be uniquely reconstructed back from the dictionary.
        
        :return: FixedScan object as python dictionary
        :rtype: dict

        """
        return dict({"@id": self._id})
    
    def __repr__(self):
        return "FixedScan.from_dict({})".format(self.to_dict())

    def __eq__(self, other):
        # Equality test is simple one which compares the data attributes.
        # note that _id data attribute may be different
        if(isinstance(self, other.__class__)):
            return True
        else:
            return NotImplemented

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
        """Parses an CrossTrackScan object from a normalized JSON dictionary.
        
        :param d: Dictionary with the cross-track scan specifications.
        :paramtype d: dict

        :return: CrossTrackScan object.
        :rtype: :class:`instrupy.radiometer_model.CrossTrackScan`

        """             
        return CrossTrackScan( scanWidth = d.get("scanWidth", None),
                               interScanOverheadTime = d.get("interScanOverheadTime", None),
                               _id = d.get("@id", None))
    
    def to_dict(self):
        """ Translate the CrossTrackScan object to a Python dictionary such that it can be uniquely reconstructed back from the dictionary.
        
        :return: CrossTrackScan object as python dictionary
        :rtype: dict

        """
        return dict({"scanWidth": self.scanWidth,
                     "interScanOverheadTime": self.interScanOverheadTime,
                     "@id": self._id
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

class ConicalScan(Entity):
    """ Class to handle conical scan.

    :ivar offNadirAngle: Off-nadir angle (i.e. the half-cone angle of the conical scan) in degrees.
    :vartype offNadirAngle: float

    :ivar clockAngleRange: Scan clock angle range in degrees.
    :vartype clockAngleRange: float

    :ivar interScanOverheadTime: Time in seconds taken from ending current scan to starting next scan. Significant in case of mechanical scanning.
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
        """Parses an ConicalScan object from a normalized JSON dictionary.
        
        :param d: Dictionary with the cross-track scan specifications.
        :paramtype d: dict

        :return: ConicalScan object.
        :rtype: :class:`instrupy.radiometer_model.ConicalScan`

        """             
        return ConicalScan( offNadirAngle = d.get("offNadirAngle", None),
                            clockAngleRange = d.get("clockAngleRange", None),
                            interScanOverheadTime = d.get("interScanOverheadTime", None),
                            _id = d.get("@id", None))
    
    def to_dict(self):
        """ Translate the ConicalScan object to a Python dictionary such that it can be uniquely reconstructed back from the dictionary.
        
        :return: ConicalScan object as python dictionary
        :rtype: dict

        """
        return dict({"offNadirAngle": self.offNadirAngle,
                     "clockAngleRange": self.clockAngleRange,
                     "interScanOverheadTime": self.interScanOverheadTime,
                     "@id": self._id
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
    
    :ivar pointing_option: List of ``Orientation`` objects which specify the orientation of the instrument pointing axis into which the instrument-axis can be maneuvered. 
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

    :ivar bandwidth: Pre-detection bandwidth in (Hertz).
    :vartype bandwidth: float

    :ivar systemType: Radiometer system type.
    :vartype: :class:`instrupy.radiometer_model.SystemType`
    
    :ivar system: Radiometer system object.
    :varytype system: :class:`instrupy.radiometer_model.TotalPowerRadiometerSystem` or :class:`instrupy.radiometer_model.UnbalancedDikeRadiometerSystem` or :class:`instrupy.radiometer_model.BalancedDikeRadiometerSystem` or :class:`instrupy.radiometer_model.NoiseAddingRadiometerSystem`
    
    :ivar scanTechnique: Scanning technique.
    :vartype scanTechnique: :class:`instrupy.radiometer_model.ScanTech`
    
    :ivar scan: Scan object.
    :vartype scan: :class:`instrupy.radiometer_model.FixedScan` or :class:`instrupy.radiometer_model.CrossTrackScan` or :class:`instrupy.radiometer_model.ConicalScan`.

    :ivar _id: Unique identifier.
    :vartype _id: str or int
          
    """
    def __init__(self, name=None, mass=None, volume=None, power=None,  orientation=None, 
            fieldOfViewGeometry=None, sceneFieldOfViewGeometry=None, maneuver=None, pointingOption=None, 
            dataRate=None, bitsPerPixel=None, antenna=None, operatingFrequency=None, 
            bandwidth=None, systemType=None, system=None, scanTechnique=None, 
            scan=None, _id=None):
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
        self.bandwidth = float(bandwidth) if bandwidth is not None else None
        self.systemType = SystemType.get(systemType) if systemType is not None else None
        self.system = copy.deepcopy(system) if system is not None and (isinstance(system, TotalPowerRadiometerSystem) or isinstance(system, UnbalancedDikeRadiometerSystem) or isinstance(system, BalancedDikeRadiometerSystem) or isinstance(system, NoiseAddingRadiometerSystem))  else None
        self.scanTechnique = ScanTech.get(scanTechnique) if scanTechnique is not None else None
        self.scan = copy.deepcopy(scan) if scan is not None and (isinstance(scan, FixedScan) or isinstance(scan, CrossTrackScan) or isinstance(scan, ConicalScan)) else None

        super(RadiometerModel,self).__init__(_id, "Radiometer")
        
    @staticmethod
    def from_dict(d):
        """ Parses an radiometer instrument from a normalized JSON dictionary.

        The following default values are assigned to the object instance parameters in case of 
        :class:`None` values or missing key/value pairs in the input dictionary.

        .. csv-table:: Default values
            :header: Parameter, Default Value
            :widths: 10,40

            scanTech, ScanTech.Fixed
            orientation, Orientation.Convention.REF_FRAME_ALIGNED
            sceneFieldOfViewGeometry, (Instrument) fieldOfViewGeometry
            _id, random string
        
        :param d: Normalized JSON dictionary with the corresponding model specifications. 
        :paramtype d: dict

        :returns: RadiometerModel object initialized with the input specifications.
        :rtype: :class:`instrupy.RadiometerModel`

        """
        # Only side-looking orientation of instrument supported for synthetic aperture radar 
        orien_dict = d.get("orientation", {"convention": "REF_FRAME_ALIGNED"})
        orientation = Orientation.from_dict(orien_dict)

        # parse maneuver
        maneuver_dict =  d.get("maneuver", None)
        maneuver = Maneuver.from_dict(maneuver_dict)

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

        # parse the system object
        system_lookup = { 'TOTAL_POWER' : TotalPowerRadiometerSystem,  'UNBALANCED_DICKE': UnbalancedDikeRadiometerSystem, 
                          'BALANCED_DICKE': BalancedDikeRadiometerSystem, 'NOISE_ADDING': NoiseAddingRadiometerSystem}
        system_dict = d.get("system", None)
        if system_dict:
            systemType = SystemType.get(system_dict["@type"])
            system = system_lookup[systemType].from_dict(system_dict)
        else:
            systemType = None
            system = None

        # parse the scan object
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
                        bandwidth = d.get("bandwidth", None),
                        systemType = systemType,
                        system = system,
                        scanTechnique = scanTechnique,
                        scan = scan,
                        _id = d.get("@id", uuid.uuid4())
                        )

    def to_dict(self):
        """ Translate the RadiometerModel object to a Python dictionary such that it can be uniquely reconstructed back from the dictionary.

        :returns: RadiometerModel specifications as python dictionary.
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
                "bandwidth": self.bandwidth,
                "system": system_dict,
                "scan": scan_dict,
                "@id": self._id
                })

    def __repr__(self):
        return "RadiometerModel.from_dict({})".format(self.to_dict())

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