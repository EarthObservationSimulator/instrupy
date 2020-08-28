import json
import numpy
import copy
import pandas, csv
from instrupy.util import Entity, Orientation, FieldOfView, MathUtilityFunctions, Constants, FileUtilityFunctions, EnumEntity

class PolTypeSAR(EnumEntity):
    """Enumeration of recognized SAR polarization types"""
    SINGLE = "SINGLE",
    COMPACT = "COMPACT",
    DUAL = "DUAL"

class DualPolPulseConfig(EnumEntity):
    """Enumeration of recognized dual-polarization pulse configurations"""
    AIRSAR = "AIRSAR",
    SMAP = "SMAP"

class SwathTypeSAR(EnumEntity):
    """Enumeration of recognized SAR swath imaging configurations"""
    FULL = "FULL",
    FIXED = "FIXED"
