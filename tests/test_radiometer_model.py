"""Unit tests for instrupy.radiometer_model.

References: [1] Chapter 6,7 in "Microwave Radar and Radiometric Remote Sensing," David Gardner Long , Fawwaz T. Ulaby 2014 

"""

import unittest
import json
import numpy as np
import sys, os

from instrupy.radiometer_model import PredetectionSectionParams, SystemParams
from instrupy.radiometer_model import RadiometerModel, SystemType, TotalPowerRadiometerSystem, UnbalancedDikeRadiometerSystem, \
                                      BalancedDikeRadiometerSystem, NoiseAddingRadiometerSystem, \
                                      ScanTech, FixedScan, CrossTrackScan, ConicalScan
from instrupy.util import Orientation, SphericalGeometry, ViewGeometry, FileUtilityFunctions, Maneuver


class TestTotalPowerRadiometerSystem(unittest.TestCase):   

    @classmethod
    def setUpClass(cls):
        cls.tpr_sys1_json = '{"tlLoss": 0.5,' \
                             '"tlPhyTemp": 290,' \
                             '"rfAmpGain": 30,' \
                             '"rfAmpInpNoiseTemp": 200,' \
                             '"rfAmpGainVariation": 10,' \
                             '"mixerGain": 23,' \
                             '"mixerInpNoiseTemp": 1200,' \
                             '"mixerGainVariation": 2,' \
                             '"ifAmpGain": 30,' \
                             '"ifAmpInputNoiseTemp": 100,' \
                             '"ifAmpGainVariation": 10,' \
                             '"integratorVoltageGain": 1,' \
                             '"integrationTime": 100e-3,' \
                             '"bandwidth": 10e6,' \
                             '"@id": 121}'
        
        cls.tpr_sys2_json = '{"predetectionGain": 83,' \
                             '"predetectionInpNoiseTemp": 200,' \
                             '"predetectionGainVariation": 2000000,' \
                             '"integrationTime": 100e-3,' \
                             '"bandwidth": 10e6}'

    def test_from_json(self):
        """ Test typical initialization of the total power radiometer system.
        """   
        # See [1] Section 7-3.1 for the source of some of the receiver specs specified below. Section 7.5 lists a normalaized gain variation specs of 10^-2.
        o = TotalPowerRadiometerSystem.from_json(self.tpr_sys1_json) 
        self.assertEqual(o._id, 121) 
        self.assertEqual(o._type, "TotalPowerRadiometerSystem")
        self.assertEqual(o.tlLoss, 0.5)
        self.assertEqual(o.tlPhyTemp, 290)
        self.assertEqual(o.rfAmpGain, 30)
        self.assertEqual(o.rfAmpInpNoiseTemp, 200)
        self.assertEqual(o.rfAmpGainVariation, 10)
        self.assertEqual(o.mixerGain, 23)
        self.assertEqual(o.mixerInpNoiseTemp, 1200)
        self.assertEqual(o.mixerGainVariation, 2)
        self.assertEqual(o.ifAmpGain, 30)
        self.assertEqual(o.ifAmpInputNoiseTemp, 100)
        self.assertEqual(o.ifAmpGainVariation, 10)
        self.assertEqual(o.integratorVoltageGain, 1)
        self.assertIsNone(o.predetectionGain)
        self.assertIsNone(o.predetectionInpNoiseTemp)
        self.assertIsNone(o.predetectionGainVariation)
        self.assertEqual(o.integrationTime, 100e-3)
        self.assertEqual(o.bandwidth, 10e6)

        o = TotalPowerRadiometerSystem.from_json(self.tpr_sys2_json) 
        self.assertIsNone(o._id) 
        self.assertEqual(o._type, "TotalPowerRadiometerSystem")
        self.assertIsNone(o.tlLoss)
        self.assertIsNone(o.tlPhyTemp)
        self.assertIsNone(o.rfAmpGain)
        self.assertIsNone(o.rfAmpInpNoiseTemp)
        self.assertIsNone(o.rfAmpGainVariation)
        self.assertIsNone(o.mixerGain)
        self.assertIsNone(o.mixerInpNoiseTemp)
        self.assertIsNone(o.mixerGainVariation)
        self.assertIsNone(o.ifAmpGain)
        self.assertIsNone(o.ifAmpInputNoiseTemp)
        self.assertIsNone(o.ifAmpGainVariation)
        self.assertIsNone(o.integratorVoltageGain)
        self.assertEqual(o.predetectionGain, 83)
        self.assertEqual(o.predetectionInpNoiseTemp, 200)
        self.assertEqual(o.predetectionGainVariation, 2000000)
        self.assertEqual(o.integrationTime, 100e-3)
        self.assertEqual(o.bandwidth, 10e6)

    def test_to_dict(self):
        o = TotalPowerRadiometerSystem.from_json(self.tpr_sys1_json) 
        self.assertEqual(o.to_dict(), {'tlLoss': 0.5, 'tlPhyTemp': 290.0, 'rfAmpGain': 30.0, 'rfAmpInpNoiseTemp': 200.0, 'rfAmpGainVariation': 10.0, 
                                        'mixerGain,': 23.0, 'mixerInpNoiseTemp': 1200.0, 'mixerGainVariation': 2.0, 'ifAmpGain': 30.0, 'ifAmpInputNoiseTemp': 100.0, 
                                        'ifAmpGainVariation': 10.0, 'integratorVoltageGain': 1.0, 'predetectionGain': None, 'predetectionInpNoiseTemp': None, 
                                        'predetectionGainVariation': None, 'integrationTime': 0.1, 'bandwidth': 10000000.0, '@id': 121, '@type': 'TOTAL_POWER'}
                        )

        o = TotalPowerRadiometerSystem.from_json(self.tpr_sys2_json) 
        self.assertEqual(o.to_dict(), {'tlLoss': None, 'tlPhyTemp': None, 'rfAmpGain': None, 'rfAmpInpNoiseTemp': None, 'rfAmpGainVariation': None, 
                                        'mixerGain,': None, 'mixerInpNoiseTemp': None, 'mixerGainVariation': None, 'ifAmpGain': None, 'ifAmpInputNoiseTemp': None, 
                                        'ifAmpGainVariation': None, 'integratorVoltageGain': None, 'predetectionGain': 83.0, 'predetectionInpNoiseTemp': 200.0, 
                                        'predetectionGainVariation': 2000000.0, 'integrationTime': 0.1, 'bandwidth': 10000000.0, '@id': None, '@type': 'TOTAL_POWER'}
                        )
        
        

    def test_compute_integration_time(self):

        self.assertEqual(TotalPowerRadiometerSystem.compute_integration_time(td=1.5, integration_time_spec=0.5), 0.5)
        self.assertEqual(TotalPowerRadiometerSystem.compute_integration_time(td=1.5, integration_time_spec=2), 1.5)
        self.assertEqual(TotalPowerRadiometerSystem.compute_integration_time(td=1.5, integration_time_spec=None), 1.5)

    def test_compute_predetection_sec_params(self):
        x = TotalPowerRadiometerSystem.compute_predetection_sec_params(predetectionBandwidth=10e6, tlLoss=0.5, tlPhyTemp=290,
                                            rfAmpGain=30, mixerGain=23, ifAmpGain=30, rfAmpGainVariation=10,  mixerGainVariation=2, ifAmpGainVariation=10,
                                            rfAmpInpNoiseTemp=200, mixerInpNoiseTemp=1200, ifAmpInputNoiseTemp=100)
        self.assertAlmostEqual(x.G, 177827941.00389218)
        self.assertAlmostEqual(x.G_p, 180510851.84124476)
        self.assertAlmostEqual(x.G_m, 175171746.5823525)
        self.assertAlmostEqual(x.T_REC_q, 261.1355769549698)
        self.assertAlmostEqual(x.B, 10000000.0)

        x = TotalPowerRadiometerSystem.compute_predetection_sec_params(predetectionBandwidth=15e6, predetectionGain=90, predetectionGainVariation=10000000, predetectionInpNoiseTemp=300)
        self.assertAlmostEqual(x.G, 1000000000)
        self.assertAlmostEqual(x.G_p, 1005000000)
        self.assertAlmostEqual(x.G_m, 995000000)
        self.assertAlmostEqual(x.T_REC_q, 300)
        self.assertAlmostEqual(x.B, 15000000.0)

        # no RF amplifier
        x = TotalPowerRadiometerSystem.compute_predetection_sec_params(predetectionBandwidth=10e6, tlLoss=0.5, tlPhyTemp=290,
                                            rfAmpGain=1, mixerGain=23, ifAmpGain=30, rfAmpGainVariation=0,  mixerGainVariation=2, ifAmpGainVariation=10,
                                            rfAmpInpNoiseTemp=0, mixerInpNoiseTemp=1200, ifAmpInputNoiseTemp=100)
        self.assertAlmostEqual(x.G, 223872.1138568339)
        self.assertAlmostEqual(x.G_p, 226119.10297269153)
        self.assertAlmostEqual(x.G_m, 221636.34492551928)
        self.assertAlmostEqual(x.T_REC_q, 1104.9756026018772)
        self.assertAlmostEqual(x.B, 10000000.0)                                

    def test_compute_system_params(self):
        pass

    def test_compute_radiometric_resolution(self):
        pass


class TestUnbalancedDikeRadiometerSystem(unittest.TestCase):   
    pass
class TestBalancedDikeRadiometerSystem(unittest.TestCase):   
    pass
class TestNoiseAddingRadiometerSystem(unittest.TestCase):   
    pass
class TestFixedScan(unittest.TestCase):   
    pass
class TestCrossTrackScan(unittest.TestCase):   
    pass
class TestConicalScan(unittest.TestCase):   
    pass
class TestRadiometerModel(unittest.TestCase):   
    pass