"""Unit tests for instrupy.radiometer_model.

References: [1] Chapter 6,7 in "Microwave Radar and Radiometric Remote Sensing," David Gardner Long , Fawwaz T. Ulaby 2014 

@TODO Include rectangular antenna tests

"""

import unittest
import json
import numpy as np
import sys, os

from deepdiff import DeepDiff

from instrupy.radiometer_model import PredetectionSectionParams, SystemParams
from instrupy.radiometer_model import RadiometerModel, SystemType, TotalPowerRadiometerSystem, UnbalancedDikeRadiometerSystem, \
                                      BalancedDikeRadiometerSystem, NoiseAddingRadiometerSystem, \
                                      ScanTech, FixedScan, CrossTrackScan, ConicalScan
from instrupy.util import Antenna, Orientation, SphericalGeometry, ViewGeometry, FileUtilityFunctions, Maneuver

class TestTotalPowerRadiometerSystem(unittest.TestCase):   

    @classmethod
    def setUpClass(cls):
        # See [1] Section 7-3.1 for the source of some of the receiver specs specified below. Section 7.5 lists a normalaized gain variation specs of 10^-2.
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
                             '"bandwidth": 10e6,' \
                             '"integratorVoltageGain": 1 }'

    def test_from_json(self):
        """ Test typical initialization of the total power radiometer system.
        """   
        o = TotalPowerRadiometerSystem.from_json(self.tpr_sys1_json) 
        self.assertIsInstance(o, TotalPowerRadiometerSystem)
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
        self.assertIsInstance(o, TotalPowerRadiometerSystem)
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
        self.assertEqual(o.integratorVoltageGain, 1)
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
                                        'ifAmpGainVariation': None, 'integratorVoltageGain': 1.0, 'predetectionGain': 83.0, 'predetectionInpNoiseTemp': 200.0, 
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
        self.assertIsInstance(x, PredetectionSectionParams)
        self.assertAlmostEqual(x.G, 177827941.00389218)
        self.assertAlmostEqual(x.G_p, 180510851.84124476)
        self.assertAlmostEqual(x.G_m, 175171746.5823525)
        self.assertAlmostEqual(x.T_REC_q, 261.1355769549698)
        self.assertAlmostEqual(x.B, 10000000.0)

        x = TotalPowerRadiometerSystem.compute_predetection_sec_params(predetectionBandwidth=15e6, predetectionGain=90, predetectionGainVariation=10000000, predetectionInpNoiseTemp=300)
        self.assertIsInstance(x, PredetectionSectionParams)
        self.assertAlmostEqual(x.G, 1000000000)
        self.assertAlmostEqual(x.G_p, 1005000000)
        self.assertAlmostEqual(x.G_m, 995000000)
        self.assertAlmostEqual(x.T_REC_q, 300)
        self.assertAlmostEqual(x.B, 15000000.0)

        # no RF amplifier
        x = TotalPowerRadiometerSystem.compute_predetection_sec_params(predetectionBandwidth=10e6, tlLoss=0.5, tlPhyTemp=290,
                                            rfAmpGain=1, mixerGain=23, ifAmpGain=30, rfAmpGainVariation=0,  mixerGainVariation=2, ifAmpGainVariation=10,
                                            rfAmpInpNoiseTemp=0, mixerInpNoiseTemp=1200, ifAmpInputNoiseTemp=100)
        self.assertIsInstance(x, PredetectionSectionParams)
        self.assertAlmostEqual(x.G, 223872.1138568339)
        self.assertAlmostEqual(x.G_p, 226119.10297269153)
        self.assertAlmostEqual(x.G_m, 221636.34492551928)
        self.assertAlmostEqual(x.T_REC_q, 1104.9756026018772)
        self.assertAlmostEqual(x.B, 10000000.0)                                

    def test_compute_system_params(self):
        antenna = Antenna.from_dict({"radiationEfficiency": 0.8, "phyTemp": 270})
        pd_sec_params = TotalPowerRadiometerSystem.compute_predetection_sec_params(predetectionBandwidth=10e6, tlLoss=0.5, tlPhyTemp=290,
                                            rfAmpGain=30, mixerGain=23, ifAmpGain=30, rfAmpGainVariation=10,  mixerGainVariation=2, ifAmpGainVariation=10,
                                            rfAmpInpNoiseTemp=200, mixerInpNoiseTemp=1200, ifAmpInputNoiseTemp=100)
        G = 180000000
        pd_sec_params = PredetectionSectionParams(G=G, G_p=G+0.01*G, G_m=G-0.01*G, T_REC_q=260, B=10e6)
        x = TotalPowerRadiometerSystem.compute_system_params(antenna, pd_sec_params, integratorVoltageGain=1000, T_A_q=290)
        self.assertIsInstance(x, SystemParams)
        self.assertAlmostEqual(x.G_s_delta/x.G_s_bar, 0.02)
        self.assertAlmostEqual(x.T_A, 286)
        self.assertAlmostEqual(x.T_SYS, 546)

    
    def test_compute_radiometric_resolution(self):
        # system 1
        antenna = Antenna.from_dict({"radiationEfficiency": 0.8, "phyTemp": 270})

        o = TotalPowerRadiometerSystem.from_json(self.tpr_sys1_json) # note that these is 100ms integration time specification
        self.assertAlmostEqual(o.compute_radiometric_resolution(td=200e-3, antenna=antenna, T_A_q=300), 16.676630237262927) 
        self.assertAlmostEqual(o.compute_radiometric_resolution(td=200e-3, antenna=antenna, T_A_q=600), 23.886384796495147)
        self.assertAlmostEqual(o.compute_radiometric_resolution(td=500e-3, antenna=antenna, T_A_q=300), 16.676630237262927)        
        self.assertAlmostEqual(o.compute_radiometric_resolution(td=50e-3, antenna=antenna, T_A_q=300), 16.685867420640534) 

        antenna = Antenna.from_dict({"radiationEfficiency": 0.8, "phyTemp": 350})        
        self.assertAlmostEqual(o.compute_radiometric_resolution(td=200e-3, antenna=antenna, T_A_q=300), 17.15728054121174) 

        antenna = Antenna.from_dict({"radiationEfficiency": 0.5, "phyTemp": 270})        
        self.assertAlmostEqual(o.compute_radiometric_resolution(td=200e-3, antenna=antenna, T_A_q=300), 16.406264441291718) # reduced radiantion-efficiency appears to make the radiometer more sensitive

        # system 2
        o = TotalPowerRadiometerSystem.from_json(self.tpr_sys2_json) # note that these is 100ms integration time specification
        antenna = Antenna.from_dict({"radiationEfficiency": 0.8, "phyTemp": 270})
        self.assertAlmostEqual(o.compute_radiometric_resolution(td=200e-3, antenna=antenna, T_A_q=300), 4.976310348842347) 
    
class TestUnbalancedDikeRadiometerSystem(unittest.TestCase):   
    @classmethod
    def setUpClass(cls):
        cls.udr_sys1_json = '{"tlLoss": 0.5,' \
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
                             '"dickeSwitchOutputNoiseTemperature": 90,' \
                             '"referenceTemperature": 300,' \
                             '"integratorVoltageGain": 1,' \
                             '"integrationTime": 1,' \
                             '"bandwidth": 100e6,' \
                             '"@id": "abc"}'
        
        # See Section 7-6, end of Pg. 282.
        cls.udr_sys2_json = '{"predetectionGain": 83,' \
                             '"predetectionInpNoiseTemp": 700,' \
                             '"predetectionGainVariation": 1995262.314968883,' \
                             '"integrationTime": 1,' \
                             '"bandwidth": 100e6,' \
                             '"referenceTemperature": 300,' \
                             '"integratorVoltageGain": 1 }'

    def test_from_json(self):
        """ Test typical initialization of the unbalanced Dicke radiometer system.
        """   
        o = UnbalancedDikeRadiometerSystem.from_json(self.udr_sys1_json) 
        self.assertIsInstance(o, UnbalancedDikeRadiometerSystem)
        self.assertEqual(o._id, "abc") 
        self.assertEqual(o._type, "UnbalancedDikeRadiometerSystem")
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
        self.assertEqual(o.dickeSwitchOutputNoiseTemperature, 90)
        self.assertEqual(o.referenceTemperature, 300)
        self.assertEqual(o.integratorVoltageGain, 1)
        self.assertIsNone(o.predetectionGain)
        self.assertIsNone(o.predetectionInpNoiseTemp)
        self.assertIsNone(o.predetectionGainVariation)
        self.assertEqual(o.integrationTime, 1)
        self.assertEqual(o.bandwidth, 100e6)

        o = UnbalancedDikeRadiometerSystem.from_json(self.udr_sys2_json) 
        self.assertIsInstance(o, UnbalancedDikeRadiometerSystem)
        self.assertIsNone(o._id) 
        self.assertEqual(o._type, "UnbalancedDikeRadiometerSystem")
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
        self.assertIsNone(o.dickeSwitchOutputNoiseTemperature)
        self.assertEqual(o.referenceTemperature, 300)
        self.assertEqual(o.integratorVoltageGain, 1)
        self.assertEqual(o.predetectionGain, 83)
        self.assertEqual(o.predetectionInpNoiseTemp, 700)
        self.assertEqual(o.predetectionGainVariation, 1995262.314968883)
        self.assertEqual(o.integrationTime, 1)
        self.assertEqual(o.bandwidth, 100e6)

    def test_to_dict(self):
        o = UnbalancedDikeRadiometerSystem.from_json(self.udr_sys1_json)
        self.assertEqual(o.to_dict(), {'tlLoss': 0.5, 'tlPhyTemp': 290.0, 'rfAmpGain': 30.0, 'rfAmpInpNoiseTemp': 200.0, 'rfAmpGainVariation': 10.0, 
                                        'mixerGain,': 23.0, 'mixerInpNoiseTemp': 1200.0, 'mixerGainVariation': 2.0, 'ifAmpGain': 30.0, 'ifAmpInputNoiseTemp': 100.0, 
                                        'ifAmpGainVariation': 10.0, 'dickeSwitchOutputNoiseTemperature':90.0, 'referenceTemperature':300.0, 'integratorVoltageGain': 1.0, 'predetectionGain': None, 'predetectionInpNoiseTemp': None, 
                                        'predetectionGainVariation': None, 'integrationTime': 1.0, 'bandwidth': 100000000.0, '@id': "abc", '@type': 'UNBALANCED_DICKE'}
                        )

        o = UnbalancedDikeRadiometerSystem.from_json(self.udr_sys2_json) 
        self.assertEqual(o.to_dict(), {'tlLoss': None, 'tlPhyTemp': None, 'rfAmpGain': None, 'rfAmpInpNoiseTemp': None, 'rfAmpGainVariation': None, 
                                        'mixerGain,': None, 'mixerInpNoiseTemp': None, 'mixerGainVariation': None, 'ifAmpGain': None, 'ifAmpInputNoiseTemp': None, 
                                        'ifAmpGainVariation': None, 'dickeSwitchOutputNoiseTemperature':None, 'referenceTemperature':300.0, 'integratorVoltageGain': 1.0, 'predetectionGain': 83.0, 'predetectionInpNoiseTemp': 700.0, 
                                        'predetectionGainVariation': 1995262.314968883, 'integrationTime': 1.0, 'bandwidth': 100000000.0, '@id': None, '@type': 'UNBALANCED_DICKE'}
                        )
    
    def test_compute_radiometric_resolution(self):
        antenna = Antenna.from_dict({"radiationEfficiency": 0.8, "phyTemp": 300})

        #################################################### System 1 ####################################################
        ############# Test with T_A  equal to the reference temperature #############
        o = UnbalancedDikeRadiometerSystem.from_json(self.udr_sys1_json) 
        self.assertAlmostEqual(o.compute_radiometric_resolution(td=1.5, antenna=antenna, T_A_q=300), 0.13022711539099396) # note that these is 1s integration time specification
        
        #################################################### System 2 ####################################################
        ############# See Section 7-6, end of Pg. 282. for truth values for the below calculation. #############
        ############# Test with T_A  equal to the reference temperature 
        o = UnbalancedDikeRadiometerSystem.from_json(self.udr_sys2_json) 
        self.assertAlmostEqual(o.compute_radiometric_resolution(td=1.5, antenna=antenna, T_A_q=300), 0.2)
        # Compare with total-power radiometer
        # Initialize a total-power radiometer with the same specifications. Note that however the predetection noise temperature shall be lower 
        # for a total-power radiometer since it does not include the Dicke switch.
        o = TotalPowerRadiometerSystem.from_json(self.udr_sys2_json) 
        self.assertAlmostEqual(o.compute_radiometric_resolution(td=1.5, antenna=antenna, T_A_q=300), 10.000499987500632)     

        ############# Test with T_A not equal to the reference temperature 
        o = UnbalancedDikeRadiometerSystem.from_json(self.udr_sys2_json) # note that these is 1s integration time specification
        antenna = Antenna.from_dict({"radiationEfficiency": 1, "phyTemp": 300}) # setting efficiency to 100% to remove effect of antenna physical temperature
        self.assertAlmostEqual(o.compute_radiometric_resolution(td=1.5, antenna=antenna, T_A_q=0), 3.0049625621627984)
        # Compare with total-power radiometer
        # Initialize a total-power radiometer with the same specifications. Note that however the predetection noise temperature shall be lower 
        # for a total-power radiometer since it does not include the Dicke switch.
        o = TotalPowerRadiometerSystem.from_json(self.udr_sys2_json) 
        self.assertAlmostEqual(o.compute_radiometric_resolution(td=1.5, antenna=antenna, T_A_q=0), 7.000349991250442) 

class TestBalancedDikeRadiometerSystem(unittest.TestCase):   
    @classmethod
    def setUpClass(cls):
        cls.bdr_sys1_json = '{"tlLoss": 0.5,' \
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
                             '"dickeSwitchOutputNoiseTemperature": 90,' \
                             '"integratorVoltageGain": 1,' \
                             '"integrationTime": 1,' \
                             '"bandwidth": 100e6,' \
                             '"@id": "abc"}'
        
        # See Section 7-6, end of Pg. 282.
        cls.bdr_sys2_json = '{"predetectionGain": 83,' \
                             '"predetectionInpNoiseTemp": 700,' \
                             '"predetectionGainVariation": 1995262.314968883,' \
                             '"integrationTime": 1,' \
                             '"bandwidth": 100e6,' \
                             '"integratorVoltageGain": 1 }'

    def test_from_json(self):
        """ Test typical initialization of the balanced Dicke radiometer system.
        """   
        o = BalancedDikeRadiometerSystem.from_json(self.bdr_sys1_json) 
        self.assertIsInstance(o, BalancedDikeRadiometerSystem)
        self.assertEqual(o._id, "abc") 
        self.assertEqual(o._type, "BalancedDikeRadiometerSystem")
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
        self.assertEqual(o.dickeSwitchOutputNoiseTemperature, 90)
        self.assertEqual(o.integratorVoltageGain, 1)
        self.assertIsNone(o.predetectionGain)
        self.assertIsNone(o.predetectionInpNoiseTemp)
        self.assertIsNone(o.predetectionGainVariation)
        self.assertEqual(o.integrationTime, 1)
        self.assertEqual(o.bandwidth, 100e6)

        o = BalancedDikeRadiometerSystem.from_json(self.bdr_sys2_json) 
        self.assertIsInstance(o, BalancedDikeRadiometerSystem)
        self.assertIsNone(o._id) 
        self.assertEqual(o._type, "BalancedDikeRadiometerSystem")
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
        self.assertIsNone(o.dickeSwitchOutputNoiseTemperature)
        self.assertEqual(o.integratorVoltageGain, 1)
        self.assertEqual(o.predetectionGain, 83)
        self.assertEqual(o.predetectionInpNoiseTemp, 700)
        self.assertEqual(o.predetectionGainVariation, 1995262.314968883)
        self.assertEqual(o.integrationTime, 1)
        self.assertEqual(o.bandwidth, 100e6)

    def test_to_dict(self):
        o = BalancedDikeRadiometerSystem.from_json(self.bdr_sys1_json)
        self.assertEqual(o.to_dict(), {'tlLoss': 0.5, 'tlPhyTemp': 290.0, 'rfAmpGain': 30.0, 'rfAmpInpNoiseTemp': 200.0, 'rfAmpGainVariation': 10.0, 
                                        'mixerGain,': 23.0, 'mixerInpNoiseTemp': 1200.0, 'mixerGainVariation': 2.0, 'ifAmpGain': 30.0, 'ifAmpInputNoiseTemp': 100.0, 
                                        'ifAmpGainVariation': 10.0, 'dickeSwitchOutputNoiseTemperature':90.0, 'integratorVoltageGain': 1.0, 'predetectionGain': None, 'predetectionInpNoiseTemp': None, 
                                        'predetectionGainVariation': None, 'integrationTime': 1.0, 'bandwidth': 100000000.0, '@id': "abc", '@type': 'BALANCED_DICKE'}
                        )

        o = BalancedDikeRadiometerSystem.from_json(self.bdr_sys2_json) 
        self.assertEqual(o.to_dict(), {'tlLoss': None, 'tlPhyTemp': None, 'rfAmpGain': None, 'rfAmpInpNoiseTemp': None, 'rfAmpGainVariation': None, 
                                        'mixerGain,': None, 'mixerInpNoiseTemp': None, 'mixerGainVariation': None, 'ifAmpGain': None, 'ifAmpInputNoiseTemp': None, 
                                        'ifAmpGainVariation': None, 'dickeSwitchOutputNoiseTemperature':None, 'integratorVoltageGain': 1.0, 'predetectionGain': 83.0, 'predetectionInpNoiseTemp': 700.0, 
                                        'predetectionGainVariation': 1995262.314968883, 'integrationTime': 1.0, 'bandwidth': 100000000.0, '@id': None, '@type': 'BALANCED_DICKE'}
                        )
    
    def test_compute_radiometric_resolution(self):
        antenna = Antenna.from_dict({"radiationEfficiency": 1, "phyTemp": 300}) # setting efficiency to 100% to remove effect of antenna physical temperature

        o = BalancedDikeRadiometerSystem.from_json(self.bdr_sys1_json) 
        self.assertAlmostEqual(o.compute_radiometric_resolution(td=1.5, antenna=antenna, T_A_q=0), 0.07022711539099395) # note that these is 1s integration time specification
        
        ############# Test with T_A not equal to the reference temperature 
        o = BalancedDikeRadiometerSystem.from_json(self.bdr_sys2_json) # note that these is 1s integration time specification
        
        self.assertAlmostEqual(o.compute_radiometric_resolution(td=1.5, antenna=antenna, T_A_q=0), 0.14)

class TestNoiseAddingRadiometerSystem(unittest.TestCase):   
    @classmethod
    def setUpClass(cls):
        cls.nar_sys1_json = '{"tlLoss": 0.5,' \
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
                             '"excessNoiseTemperature": 1000,' \
                             '"integratorVoltageGain": 1,' \
                             '"integrationTime": 1,' \
                             '"bandwidth": 100e6,' \
                             '"@id": "abc"}'
        
        # See Section 7-6, end of Pg. 282.
        cls.nar_sys2_json = '{"predetectionGain": 83,' \
                             '"predetectionInpNoiseTemp": 700,' \
                             '"predetectionGainVariation": 1995262.314968883,' \
                             '"excessNoiseTemperature": 10000,' \
                             '"integrationTime": 1,' \
                             '"bandwidth": 100e6,' \
                             '"integratorVoltageGain": 1 }'

    def test_from_json(self):
        """ Test typical initialization of the noise-adding radiometer system.
        """   
        o = NoiseAddingRadiometerSystem.from_json(self.nar_sys1_json) 
        self.assertIsInstance(o, NoiseAddingRadiometerSystem)
        self.assertEqual(o._id, "abc") 
        self.assertEqual(o._type, "NoiseAddingRadiometerSystem")
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
        self.assertEqual(o.excessNoiseTemperature, 1000)
        self.assertEqual(o.integratorVoltageGain, 1)
        self.assertIsNone(o.predetectionGain)
        self.assertIsNone(o.predetectionInpNoiseTemp)
        self.assertIsNone(o.predetectionGainVariation)
        self.assertEqual(o.integrationTime, 1)
        self.assertEqual(o.bandwidth, 100e6)

        o = NoiseAddingRadiometerSystem.from_json(self.nar_sys2_json) 
        self.assertIsInstance(o, NoiseAddingRadiometerSystem)
        self.assertIsNone(o._id) 
        self.assertEqual(o._type, "NoiseAddingRadiometerSystem")
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
        self.assertEqual(o.excessNoiseTemperature, 10000)
        self.assertEqual(o.integratorVoltageGain, 1)
        self.assertEqual(o.predetectionGain, 83)
        self.assertEqual(o.predetectionInpNoiseTemp, 700)
        self.assertEqual(o.predetectionGainVariation, 1995262.314968883)
        self.assertEqual(o.integrationTime, 1)
        self.assertEqual(o.bandwidth, 100e6)

    def test_to_dict(self):
        o = NoiseAddingRadiometerSystem.from_json(self.nar_sys1_json)
        self.assertEqual(o.to_dict(), {'tlLoss': 0.5, 'tlPhyTemp': 290.0, 'rfAmpGain': 30.0, 'rfAmpInpNoiseTemp': 200.0, 'rfAmpGainVariation': 10.0, 
                                        'mixerGain,': 23.0, 'mixerInpNoiseTemp': 1200.0, 'mixerGainVariation': 2.0, 'ifAmpGain': 30.0, 'ifAmpInputNoiseTemp': 100.0, 
                                        'ifAmpGainVariation': 10.0, 'excessNoiseTemperature':1000.0, 'integratorVoltageGain': 1.0, 'predetectionGain': None, 'predetectionInpNoiseTemp': None, 
                                        'predetectionGainVariation': None, 'integrationTime': 1.0, 'bandwidth': 100000000.0, '@id': "abc", '@type': 'NOISE_ADDING'}
                        )

        o = NoiseAddingRadiometerSystem.from_json(self.nar_sys2_json) 
        self.assertEqual(o.to_dict(), {'tlLoss': None, 'tlPhyTemp': None, 'rfAmpGain': None, 'rfAmpInpNoiseTemp': None, 'rfAmpGainVariation': None, 
                                        'mixerGain,': None, 'mixerInpNoiseTemp': None, 'mixerGainVariation': None, 'ifAmpGain': None, 'ifAmpInputNoiseTemp': None, 
                                        'ifAmpGainVariation': None, 'excessNoiseTemperature':10000.0, 'integratorVoltageGain': 1.0, 'predetectionGain': 83.0, 'predetectionInpNoiseTemp': 700.0, 
                                        'predetectionGainVariation': 1995262.314968883, 'integrationTime': 1.0, 'bandwidth': 100000000.0, '@id': None, '@type': 'NOISE_ADDING'}
                        )

    def test_compute_radiometric_resolution(self):
        antenna = Antenna.from_dict({"radiationEfficiency": 0.8, "phyTemp": 300})

        o = NoiseAddingRadiometerSystem.from_json(self.nar_sys1_json) 
        self.assertAlmostEqual(o.compute_radiometric_resolution(td=1.5, antenna=antenna, T_A_q=300), 0.23817636968082867) # note that these is 1s integration time specification

        o = NoiseAddingRadiometerSystem.from_json(self.nar_sys2_json) 
        self.assertAlmostEqual(o.compute_radiometric_resolution(td=1.5, antenna=antenna, T_A_q=300), 0.24) # note that these is 1s integration time specification


class TestFixedScan(unittest.TestCase):   
    def test_from_json(self):
        """ Test typical initialization of the FixedScan object
        """   
        o = FixedScan.from_json('{"@id": 123}') 
        self.assertIsInstance(o, FixedScan)
        self.assertEqual(o._id, 123) 
        self.assertEqual(o._type, "FixedScan")

        o = FixedScan.from_json('{"@id": "abc"}') 
        self.assertIsInstance(o, FixedScan)
        self.assertEqual(o._id, "abc") 
        self.assertEqual(o._type, "FixedScan")

        o = FixedScan.from_json('{}') 
        self.assertIsInstance(o, FixedScan)
        self.assertIsNone(o._id) 
        self.assertEqual(o._type, "FixedScan")
    
    def test_to_dict(self):
        o = FixedScan.from_json('{"@id": 123}')
        self.assertEqual(o.to_dict(), {'@id': 123, '@type': 'FIXED'})

        o = FixedScan.from_json('{"@id": "abc"}')
        self.assertEqual(o.to_dict(), {'@id': "abc", '@type': 'FIXED'})

        o = FixedScan.from_json('{}')
        self.assertEqual(o.to_dict(), {'@id': None, '@type': 'FIXED'})
    
    def test_compute_instru_field_of_view(self):
        o = FixedScan.from_json('{"@id": "abc"}')
        instru_orientation = Orientation.from_dict({"referenceFrame": "SC_BODY_FIXED", "convention": "SIDE_LOOK","sideLookAngle":10})

        antenna_fov_sph_geom = SphericalGeometry.from_dict({"shape": "CIRCULAR", "diameter": 30})
        instru_fov_sph_geom = antenna_fov_sph_geom
        self.assertEqual(o.compute_instru_field_of_view(antenna_fov_sph_geom=antenna_fov_sph_geom, instru_orientation=instru_orientation), ViewGeometry(orien=instru_orientation, sph_geom=instru_fov_sph_geom))

        antenna_fov_sph_geom = SphericalGeometry.from_dict({"shape": "RECTANGULAR", "angleHeight": 10, "angleWidth": 20})
        instru_fov_sph_geom = antenna_fov_sph_geom
        self.assertEqual(o.compute_instru_field_of_view(antenna_fov_sph_geom=antenna_fov_sph_geom, instru_orientation=instru_orientation), ViewGeometry(orien=instru_orientation, sph_geom=instru_fov_sph_geom))

    def test_compute_dwell_time_per_ground_pixel(self):
        o = FixedScan.from_json('{"@id": 123}') 
        self.assertAlmostEqual(o.compute_dwell_time_per_ground_pixel(res_AT_m=1000, sat_speed_kmps=7.8), 0.1282051282051282)
    
    def test_compute_swath_width(self):
        o = FixedScan.from_json('{"@id": 123}')
        antenna_fov_sph_geom = SphericalGeometry.from_dict({"shape": "CIRCULAR", "diameter": 30})
        # using approximate swath formula as the truth data                                      
        self.assertAlmostEqual(o.compute_swath_width(alt_km=500, instru_look_angle_deg=0, antenna_fov_sph_geom=antenna_fov_sph_geom), 30*np.pi/180*500, delta=25) 
        self.assertAlmostEqual(o.compute_swath_width(alt_km=700, instru_look_angle_deg=0, antenna_fov_sph_geom=antenna_fov_sph_geom), 30*np.pi/180*700, delta=25)
        self.assertAlmostEqual(o.compute_swath_width(alt_km=500, instru_look_angle_deg=15, antenna_fov_sph_geom=antenna_fov_sph_geom), 30*np.pi/180*(500/np.cos(np.deg2rad(15))), delta=25)

class TestCrossTrackScan(unittest.TestCase):   
    def test_from_json(self):
        """ Test typical initialization of the CrossTrackScan object
        """   
        o = CrossTrackScan.from_json('{"@id": 123, "scanWidth": 120, "interScanOverheadTime": 1e-3}') 
        self.assertIsInstance(o, CrossTrackScan)
        self.assertEqual(o._id, 123) 
        self.assertEqual(o._type, "CrossTrackScan")
        self.assertEqual(o.scanWidth, 120) 
        self.assertEqual(o.interScanOverheadTime, 1e-3) 
    
    def test_to_dict(self):
        o = CrossTrackScan.from_json('{"@id": 123, "scanWidth": 120, "interScanOverheadTime": 1e-3}') 
        self.assertEqual(o.to_dict(), {'@id': 123, '@type': 'CROSS_TRACK', "scanWidth": 120.0, "interScanOverheadTime": 0.001})
    
    def test_compute_instru_field_of_view(self):
        o = CrossTrackScan.from_json('{"@id": 123, "scanWidth": 120, "interScanOverheadTime": 1e-3}') 
        instru_orientation = Orientation.from_dict({"referenceFrame": "SC_BODY_FIXED", "convention": "SIDE_LOOK","sideLookAngle":10})

        antenna_fov_sph_geom = SphericalGeometry.from_dict({"shape": "CIRCULAR", "diameter": 30})
        instru_fov_sph_geom = SphericalGeometry.from_dict({"shape": "RECTANGULAR", "angleHeight": 30, "angleWidth": 150})
        self.assertEqual(o.compute_instru_field_of_view(antenna_fov_sph_geom=antenna_fov_sph_geom, instru_orientation=instru_orientation), ViewGeometry(orien=instru_orientation, sph_geom=instru_fov_sph_geom))

        antenna_fov_sph_geom = SphericalGeometry.from_dict({"shape": "RECTANGULAR", "angleHeight": 16, "angleWidth": 50})
        instru_fov_sph_geom = SphericalGeometry.from_dict({"shape": "RECTANGULAR", "angleHeight": 16, "angleWidth": 170})
        ddiff = DeepDiff(o.compute_instru_field_of_view(antenna_fov_sph_geom=antenna_fov_sph_geom, instru_orientation=instru_orientation), 
                         ViewGeometry(orien=instru_orientation, sph_geom=instru_fov_sph_geom), 
                         ignore_numeric_type_changes=True)
        self.assertEqual(ddiff, {})

    def test_compute_dwell_time_per_ground_pixel(self):
        o = CrossTrackScan.from_json('{"@id": 123, "scanWidth": 120, "interScanOverheadTime": 1e-3}') 
        self.assertAlmostEqual(o.compute_dwell_time_per_ground_pixel(res_AT_m=5000, sat_speed_kmps=7.8, iFOV_CT_deg=4), 0.021334188034188035)        
        self.assertAlmostEqual(o.compute_dwell_time_per_ground_pixel(res_AT_m=10000, sat_speed_kmps=7.8, iFOV_CT_deg=4), 2*0.021334188034188035, places=4) # dwell time should be around doubled in case of double along-track pixel resolution
        self.assertAlmostEqual(o.compute_dwell_time_per_ground_pixel(res_AT_m=5000, sat_speed_kmps=7.8, iFOV_CT_deg=8), 2*0.021334188034188035, places=4) # dwell time should be around doubled in case of cross-track iFOV

        o = CrossTrackScan.from_json('{"@id": 123, "scanWidth": 120, "interScanOverheadTime": 10e-3}') 
        self.assertAlmostEqual(o.compute_dwell_time_per_ground_pixel(res_AT_m=5000, sat_speed_kmps=7.8, iFOV_CT_deg=4), 0.021034188034188037) 

    def test_compute_swath_width(self):
        o = CrossTrackScan.from_json('{"@id": 123, "scanWidth": 20, "interScanOverheadTime": 1e-3}') 
        antenna_fov_sph_geom = SphericalGeometry.from_dict({"shape": "CIRCULAR", "diameter": 1})

        # using approximate swath formula as the truth data                                      
        self.assertAlmostEqual(o.compute_swath_width(alt_km=500, instru_look_angle_deg=0, antenna_fov_sph_geom=antenna_fov_sph_geom), 20*np.pi/180*500, delta=25) 
        self.assertAlmostEqual(o.compute_swath_width(alt_km=700, instru_look_angle_deg=0, antenna_fov_sph_geom=antenna_fov_sph_geom), 20*np.pi/180*700, delta=25) 

        o = CrossTrackScan.from_json('{"@id": 123, "scanWidth": 60, "interScanOverheadTime": 1e-3}') 
        self.assertAlmostEqual(o.compute_swath_width(alt_km=500, instru_look_angle_deg=0, antenna_fov_sph_geom=antenna_fov_sph_geom), 60*np.pi/180*500, delta=75) 

class TestConicalScan(unittest.TestCase):   
    def test_from_json(self):
        """ Test typical initialization of the ConicalScan object
        """   
        o = ConicalScan.from_json('{"@id": "abc", "offNadirAngle": 30, "clockAngleRange": 60, "interScanOverheadTime": 1e-3}') 
        self.assertIsInstance(o, ConicalScan)
        self.assertEqual(o._id, "abc") 
        self.assertEqual(o._type, "ConicalScan")
        self.assertEqual(o.offNadirAngle, 30)
        self.assertEqual(o.clockAngleRange, 60)
        self.assertEqual(o.interScanOverheadTime, 1e-3) 
    
    def test_to_dict(self):
        o = ConicalScan.from_json('{"@id": "abc", "offNadirAngle": 30, "clockAngleRange": 60, "interScanOverheadTime": 1e-3}') 
        self.assertEqual(o.to_dict(), {'@id': "abc", '@type': 'CONICAL', "offNadirAngle": 30.0, "clockAngleRange": 60.0, "interScanOverheadTime": 0.001})
    
    def test_compute_instru_field_of_view(self):
        o = ConicalScan.from_json('{"@id": "abc", "offNadirAngle": 30, "clockAngleRange": 60, "interScanOverheadTime": 1e-3}')        
        instru_orientation = Orientation.from_dict({"referenceFrame": "SC_BODY_FIXED", "convention": "SIDE_LOOK","sideLookAngle":10})

        antenna_fov_sph_geom = SphericalGeometry.from_dict({"shape": "CIRCULAR", "diameter": 30})
        with self.assertRaises(NotImplementedError):
            o.compute_instru_field_of_view(antenna_fov_sph_geom=antenna_fov_sph_geom, instru_orientation=instru_orientation)

    def test_compute_dwell_time_per_ground_pixel(self):

        # results are the same as that of the CrossTrackScan
        o = ConicalScan.from_json('{"@id": "abc", "offNadirAngle": 30, "clockAngleRange": 120, "interScanOverheadTime": 1e-3}') 
        self.assertAlmostEqual(o.compute_dwell_time_per_ground_pixel(res_AT_m=5000, sat_speed_kmps=7.8, iFOV_CT_deg=4), 0.021334188034188035)
        self.assertAlmostEqual(o.compute_dwell_time_per_ground_pixel(res_AT_m=10000, sat_speed_kmps=7.8, iFOV_CT_deg=4), 2*0.021334188034188035, places=4) # dwell time should be around doubled in case of double along-track pixel resolution
        self.assertAlmostEqual(o.compute_dwell_time_per_ground_pixel(res_AT_m=5000, sat_speed_kmps=7.8, iFOV_CT_deg=8), 2*0.021334188034188035, places=4) # dwell time should be around doubled in case of cross-track iFOV

        o = CrossTrackScan.from_json('{"scanWidth": 120, "interScanOverheadTime": 10e-3}') 
        self.assertAlmostEqual(o.compute_dwell_time_per_ground_pixel(res_AT_m=5000, sat_speed_kmps=7.8, iFOV_CT_deg=4), 0.021034188034188037) 

    def test_compute_swath_width(self):
        o = ConicalScan.from_json('{"@id": "abc", "offNadirAngle": 30, "clockAngleRange": 120, "interScanOverheadTime": 1e-3}') 

        # using approximate swath formula as the truth data                                      
        self.assertAlmostEqual(o.compute_swath_width(alt_km=500, instru_look_angle_deg=0), 612.7169711748869) 
        self.assertAlmostEqual(o.compute_swath_width(alt_km=700, instru_look_angle_deg=0), 862.5336432436297) 

        with self.assertRaises(Exception):
            o.compute_swath_width(alt_km=500, instru_look_angle_deg=30) # instrument look angle is not 0 degrees

class TestRadiometerModel(unittest.TestCase):   
    
    @classmethod
    def setUpClass(cls):

        cls.radio1_json = '{"@type": "Radiometer", "name": "ray1", "mass": 50, "volume": 3, "power": 10,' \
                          ' "orientation": {"referenceFrame": "SC_BODY_FIXED", "convention": "REF_FRAME_ALIGNED"},' \
                          ' "bitsPerPixel": 16,' \
                          ' "operatingFrequency": 1.25e9,' \
                          ' "antenna": {"shape": "CIRCULAR", "diameter": 1, "apertureExcitationProfile": "UNIFORM",' \
                          '             "radiationEfficiency": 0.8, "phyTemp": 300},' \
                          ' "system": {"tlLoss": 0.5, "tlPhyTemp": 290, ' \
                          '            "rfAmpGain": 30, "rfAmpInpNoiseTemp": 200, ' \
                          '            "rfAmpGainVariation": 10, "mixerGain": 23, "mixerInpNoiseTemp": 1200,' \
                          '            "mixerGainVariation": 2, "ifAmpGain": 30, "ifAmpInputNoiseTemp": 100,' \
                          '            "ifAmpGainVariation": 10, "integratorVoltageGain": 1, "integrationTime": 100e-3,' \
                          '            "bandwidth": 10e6, "@type": "TOTAL_POWER"},' \
                          ' "scan": {"@type": "FIXED"},' \
                          ' "targetBrightnessTemp": 345' \
                          '}'

        cls.radio2_json = '{"@type": "Radiometer", "name": "ray2", "mass": 50, ' \
                          ' "operatingFrequency": 1.25e9,' \
                          ' "antenna": {"shape": "RECTANGULAR", "height": 1, "width": 1, "apertureExcitationProfile": "UNIFORM",' \
                          '             "radiationEfficiency": 0.75, "phyTemp": 300},' \
                          ' "system": { "predetectionGain": 83, "predetectionInpNoiseTemp": 700, ' \
                          '             "predetectionGainVariation": 1995262.314968883, "integrationTime": 1, ' \
                          '             "bandwidth": 100e6, "referenceTemperature": 300, "integratorVoltageGain": 1,' \
                          '             "@type": "UNBALANCED_DICKE"},' \
                          ' "scan": {"@type": "CROSS_TRACK", "scanWidth": 120, "interScanOverheadTime": 1e-3},' \
                          ' "targetBrightnessTemp": 301' \
                          '}'

        cls.radio3_json = '{"@type": "Radiometer", "@id": "ray3",' \
                          ' "orientation": {"referenceFrame": "SC_BODY_FIXED", "convention": "REF_FRAME_ALIGNED"},' \
                          ' "bitsPerPixel": 16,' \
                          ' "operatingFrequency": 1.25e9,' \
                          ' "antenna": {"shape": "CIRCULAR", "diameter": 3.5, "apertureExcitationProfile": "UNIFORM",' \
                          '             "radiationEfficiency": 1, "phyTemp": 300},' \
                          ' "system": { "tlLoss": 0.5, "tlPhyTemp": 290, "rfAmpGain": 30, "rfAmpInpNoiseTemp": 200,' \
                          '             "rfAmpGainVariation": 10, "mixerGain": 23, "mixerInpNoiseTemp": 1200, "mixerGainVariation": 2,' \
                          '             "ifAmpGain": 30, "ifAmpInputNoiseTemp": 100, "ifAmpGainVariation": 10, "dickeSwitchOutputNoiseTemperature": 90,' \
                          '             "integratorVoltageGain": 1, "integrationTime": 1, "bandwidth": 100e6, "@type": "BALANCED_DICKE"},' \
                          ' "scan": {"@type": "CONICAL", "offNadirAngle": 30, "clockAngleRange": 60, "interScanOverheadTime": 1e-3},' \
                          ' "targetBrightnessTemp": 295' \
                          '}'
        
        cls.radio4_json = '{"@type": "Radiometer", "@id": "ray4",' \
                          ' "orientation": {"referenceFrame": "SC_BODY_FIXED", "convention": "SIDE_LOOK", "sideLookAngle":-30},' \
                          ' "operatingFrequency": 1.25e9,' \
                          ' "antenna": {"shape": "CIRCULAR", "diameter": 1, "apertureExcitationProfile": "UNIFORM",' \
                          '             "radiationEfficiency": 1, "phyTemp": 300},' \
                          ' "system": { "tlLoss": 0.5, "tlPhyTemp": 290, "rfAmpGain": 30, "rfAmpInpNoiseTemp": 200,' \
                          '             "rfAmpGainVariation": 10, "mixerGain": 23, "mixerInpNoiseTemp": 1200, "mixerGainVariation": 2,' \
                          '             "ifAmpGain": 30, "ifAmpInputNoiseTemp": 100, "ifAmpGainVariation": 10, "excessNoiseTemperature": 1000,' \
                          '             "integratorVoltageGain": 1, "integrationTime": 1, "bandwidth": 100e6, "@type": "NOISE_ADDING"},' \
                          ' "scan": {"@type": "FIXED"},' \
                          ' "targetBrightnessTemp": 295' \
                          '}'


    def test_from_json(self):
        """ Test typical initializations of the RadiometerModel object
        """   
        o = RadiometerModel.from_json(self.radio1_json)
        self.assertIsInstance(o, RadiometerModel)
        self.assertIsNotNone(o._id)
        self.assertEqual(o.name, "ray1")
        self.assertEqual(o.mass, 50)
        self.assertEqual(o.volume, 3)
        self.assertEqual(o.power, 10)
        self.assertEqual(o.orientation, Orientation.from_dict({"referenceFrame": "SC_BODY_FIXED", "convention": "REF_FRAME_ALIGNED"}))
        self.assertEqual(o.bitsPerPixel, 16)
        self.assertEqual(o.operatingFrequency, 1.25e9)
        self.assertEqual(o.antenna, Antenna.from_dict({"shape": "CIRCULAR", "diameter": 1, "apertureExcitationProfile": "UNIFORM", "radiationEfficiency": 0.8, "phyTemp": 300}))
        self.assertEqual(o.system, TotalPowerRadiometerSystem.from_dict({"tlLoss": 0.5, "tlPhyTemp": 290, "rfAmpGain": 30, "rfAmpInpNoiseTemp": 200, "rfAmpGainVariation": 10, "mixerGain": 23, "mixerInpNoiseTemp": 1200, "mixerGainVariation": 2, 
                                                                         "ifAmpGain": 30, "ifAmpInputNoiseTemp": 100, "ifAmpGainVariation": 10, 
                                                                         "integratorVoltageGain": 1, "integrationTime": 100e-3, "bandwidth": 10e6}))
        self.assertEqual(o.scan, FixedScan.from_dict({}))
        self.assertEqual(o.targetBrightnessTemp, 345)

        o = RadiometerModel.from_json(self.radio2_json)
        self.assertIsInstance(o, RadiometerModel)
        self.assertIsNotNone(o._id)
        self.assertEqual(o.name, "ray2")
        self.assertEqual(o.mass, 50)
        self.assertIsNone(o.volume)
        self.assertIsNone(o.power)
        self.assertEqual(o.orientation, Orientation.from_dict({"referenceFrame": "SC_BODY_FIXED", "convention": "REF_FRAME_ALIGNED"}))
        self.assertIsNone(o.bitsPerPixel)
        self.assertEqual(o.operatingFrequency, 1.25e9)
        self.assertEqual(o.antenna, Antenna.from_dict({"shape": "RECTANGULAR", "height": 1, "width":1, "apertureExcitationProfile": "UNIFORM", "radiationEfficiency": 0.75, "phyTemp": 300}))
        self.assertEqual(o.system, UnbalancedDikeRadiometerSystem.from_dict({ "predetectionGain": 83, "predetectionInpNoiseTemp": 700,
                                                                          "predetectionGainVariation": 1995262.314968883, "integrationTime": 1,
                                                                          "bandwidth": 100e6, "referenceTemperature": 300, "integratorVoltageGain": 1,
                                                                          "@type": "UNBALANCED_DICKE"}))
        self.assertEqual(o.scan, CrossTrackScan.from_dict({"scanWidth": 120, "interScanOverheadTime": 1e-3}))
        self.assertEqual(o.targetBrightnessTemp, 301)         

        o = RadiometerModel.from_json(self.radio3_json)
        self.assertIsInstance(o, RadiometerModel)
        self.assertEqual(o._id, "ray3")
        self.assertIsNone(o.name)
        self.assertIsNone(o.mass)
        self.assertIsNone(o.volume)
        self.assertIsNone(o.power)
        self.assertEqual(o.orientation, Orientation.from_dict({"referenceFrame": "SC_BODY_FIXED", "convention": "REF_FRAME_ALIGNED"}))
        self.assertEqual(o.bitsPerPixel, 16)
        self.assertEqual(o.operatingFrequency, 1.25e9)
        self.assertEqual(o.antenna, Antenna.from_dict({"shape": "CIRCULAR", "diameter": 3.5, "apertureExcitationProfile": "UNIFORM", "radiationEfficiency": 1, "phyTemp": 300}))
        self.assertEqual(o.system, BalancedDikeRadiometerSystem.from_dict({ "tlLoss": 0.5, "tlPhyTemp": 290, "rfAmpGain": 30, "rfAmpInpNoiseTemp": 200,
                                                                              "rfAmpGainVariation": 10, "mixerGain": 23, "mixerInpNoiseTemp": 1200, "mixerGainVariation": 2,
                                                                              "ifAmpGain": 30, "ifAmpInputNoiseTemp": 100, "ifAmpGainVariation": 10, "dickeSwitchOutputNoiseTemperature": 90,
                                                                              "integratorVoltageGain": 1, "integrationTime": 1, "bandwidth": 100e6, "@type": "BALANCED_DICKE"}))
        self.assertEqual(o.scan, ConicalScan.from_dict({"offNadirAngle": 30, "clockAngleRange": 60, "interScanOverheadTime": 1e-3}))
        self.assertEqual(o.targetBrightnessTemp, 295)

        o = RadiometerModel.from_json(self.radio4_json)
        self.assertIsInstance(o, RadiometerModel)
        self.assertEqual(o._id, "ray4")
        self.assertIsNone(o.name)
        self.assertIsNone(o.mass)
        self.assertIsNone(o.volume)
        self.assertIsNone(o.power)
        self.assertEqual(o.orientation, Orientation.from_dict({"referenceFrame": "SC_BODY_FIXED", "convention": "SIDE_LOOK", "sideLookAngle":-30}))
        self.assertIsNone(o.bitsPerPixel)
        self.assertEqual(o.operatingFrequency, 1.25e9)
        self.assertEqual(o.antenna, Antenna.from_dict({"shape": "CIRCULAR", "diameter": 1, "apertureExcitationProfile": "UNIFORM", "radiationEfficiency": 1, "phyTemp": 300}))
        self.assertEqual(o.system, NoiseAddingRadiometerSystem.from_dict({ "tlLoss": 0.5, "tlPhyTemp": 290, "rfAmpGain": 30, "rfAmpInpNoiseTemp": 200,
                                                                              "rfAmpGainVariation": 10, "mixerGain": 23, "mixerInpNoiseTemp": 1200, "mixerGainVariation": 2,
                                                                              "ifAmpGain": 30, "ifAmpInputNoiseTemp": 100, "ifAmpGainVariation": 10, "excessNoiseTemperature": 1000,
                                                                              "integratorVoltageGain": 1, "integrationTime": 1, "bandwidth": 100e6, "@type": "NOISE_ADDING"}))
        self.assertEqual(o.scan, FixedScan.from_dict({}))
        self.assertEqual(o.targetBrightnessTemp, 295)
    
    def test_to_dict(self):
        o = RadiometerModel.from_json(self.radio1_json)
        _id = o._id # id is generated randomly
        self.assertEqual(o.to_dict(), {'@type': 'Radiometer', 'name': 'ray1', 'mass': 50.0, 'volume': 3.0, 'power': 10.0, 
                                       'orientation': {'referenceFrame': 'SC_BODY_FIXED', 'convention': 'EULER', 'eulerAngle1': 0.0, 'eulerAngle2': 0.0, 'eulerAngle3': 0.0, 'eulerSeq1': 1, 'eulerSeq2': 2, 'eulerSeq3': 3, '@id': None}, 'fieldOfViewGeometry': {'shape': 'CIRCULAR', 'diameter': 13.741474058602394, '@id': None}, 
                                       'sceneFieldOfViewGeometry': {'shape': 'CIRCULAR', 'diameter': 13.741474058602394, '@id': None}, 
                                       'maneuver': None, 'pointingOption': None, 'dataRate': None, 'bitsPerPixel': 16, 
                                       'antenna': {'shape': 'CIRCULAR', 'apertureExcitationProfile': 'UNIFORM', 'diameter': 1.0, 'height': None, 'width': None, 'apertureEfficiency': None, 'radiationEfficiency': 0.8, 'phyTemp': 300.0, '@id': None}, 
                                       'operatingFrequency': 1250000000.0, 
                                       'system': {'tlLoss': 0.5, 'tlPhyTemp': 290.0, 'rfAmpGain': 30.0, 'rfAmpInpNoiseTemp': 200.0, 'rfAmpGainVariation': 10.0, 
                                                  'mixerGain,': 23.0, 'mixerInpNoiseTemp': 1200.0, 'mixerGainVariation': 2.0, 
                                                  'ifAmpGain': 30.0, 'ifAmpInputNoiseTemp': 100.0, 'ifAmpGainVariation': 10.0, 
                                                  'integratorVoltageGain': 1.0, 'predetectionGain': None, 'predetectionInpNoiseTemp': None, 
                                                  'predetectionGainVariation': None, 'integrationTime': 0.1, 'bandwidth': 10000000.0, '@id': None, '@type': 'TOTAL_POWER'}, 
                                        'scan': {'@id': None, '@type': 'FIXED'}, 'targetBrightnessTemp': 345.0, '@id': _id})

        o = RadiometerModel.from_json(self.radio2_json)
        _id = o._id # id is generated randomly

        ddiff = DeepDiff(o.to_dict(), {'@type': 'Radiometer', 'name': 'ray2', 'mass': 50.0, 'volume': None, 'power': None, 
                                       'orientation': {'referenceFrame': 'SC_BODY_FIXED', 'convention': 'EULER', 'eulerAngle1': 0.0, 'eulerAngle2': 0.0, 'eulerAngle3': 0.0, 'eulerSeq1': 1, 'eulerSeq2': 2, 'eulerSeq3': 3, '@id': None}, 
                                       'fieldOfViewGeometry': {'shape': 'RECTANGULAR', 'angleHeight': 12.092497171570107, 'angleWidth': 12.092497171570086, '@id': None}, 
                                       'sceneFieldOfViewGeometry': {'shape': 'RECTANGULAR', 'angleHeight': 12.092497171570107, 'angleWidth': 12.092497171570086, '@id': None}, 
                                       'maneuver': None, 'pointingOption': None, 'dataRate': None, 'bitsPerPixel': None, 
                                       'antenna': {'shape': 'RECTANGULAR', 'apertureExcitationProfile': 'UNIFORM', 'diameter': None, 'height': 1.0, 'width': 1.0, 'apertureEfficiency': None, 'radiationEfficiency': 0.75, 'phyTemp': 300.0, '@id': None}, 
                                       'operatingFrequency': 1250000000.0, 
                                       'system': {'tlLoss': None, 'tlPhyTemp': None, 'rfAmpGain': None, 'rfAmpInpNoiseTemp': None, 'rfAmpGainVariation': None, 
                                                  'mixerGain,': None, 'mixerInpNoiseTemp': None, 'mixerGainVariation': None, 
                                                  'ifAmpGain': None, 'ifAmpInputNoiseTemp': None, 'ifAmpGainVariation': None, 
                                                  'dickeSwitchOutputNoiseTemperature': None, 'referenceTemperature': 300.0, 
                                                  'integratorVoltageGain': 1.0, 'predetectionGain': 83.0, 'predetectionInpNoiseTemp': 700.0, 'predetectionGainVariation': 1995262.314968883, 
                                                  'integrationTime': 1.0, 'bandwidth': 100000000.0, '@id': None, '@type': 'UNBALANCED_DICKE'}, 
                                       'scan': {'scanWidth': 120.0, 'interScanOverheadTime': 0.001, '@id': None, '@type': 'CROSS_TRACK'}, 
                                       'targetBrightnessTemp': 301.0, '@id': _id}, ignore_numeric_type_changes=True)
        self.assertEqual(ddiff, {})

        o = RadiometerModel.from_json(self.radio3_json)
        self.assertEqual(o.to_dict(), {'@type': 'Radiometer', 'name': None, 'mass': None, 'volume': None, 'power': None, 
                            'orientation': {'referenceFrame': 'SC_BODY_FIXED', 'convention': 'EULER', 'eulerAngle1': 0.0, 'eulerAngle2': 0.0, 'eulerAngle3': 0.0, 'eulerSeq1': 1, 'eulerSeq2': 2, 'eulerSeq3': 3, '@id': None}, 
                            'fieldOfViewGeometry': {'shape': 'CIRCULAR', 'diameter': 3.9261354453149697, '@id': None}, 
                            'sceneFieldOfViewGeometry': {'shape': 'CIRCULAR', 'diameter': 3.9261354453149697, '@id': None}, 
                            'maneuver': None, 'pointingOption': None, 'dataRate': None, 'bitsPerPixel': 16, 
                            'antenna': {'shape': 'CIRCULAR', 'apertureExcitationProfile': 'UNIFORM', 'diameter': 3.5, 'height': None, 'width': None, 'apertureEfficiency': None, 'radiationEfficiency': 1.0, 'phyTemp': 300.0, '@id': None}, 
                            'operatingFrequency': 1250000000.0, 
                            'system': {'tlLoss': 0.5, 'tlPhyTemp': 290.0, 'rfAmpGain': 30.0, 'rfAmpInpNoiseTemp': 200.0, 'rfAmpGainVariation': 10.0, 
                                       'mixerGain,': 23.0, 'mixerInpNoiseTemp': 1200.0, 'mixerGainVariation': 2.0, 
                                       'ifAmpGain': 30.0, 'ifAmpInputNoiseTemp': 100.0, 'ifAmpGainVariation': 10.0, 
                                       'dickeSwitchOutputNoiseTemperature': 90.0, 'integratorVoltageGain': 1.0, 
                                       'predetectionGain': None, 'predetectionInpNoiseTemp': None, 'predetectionGainVariation': None, 
                                       'integrationTime': 1.0, 'bandwidth': 100000000.0, '@id': None, '@type': 'BALANCED_DICKE'}, 
                            'scan': {'offNadirAngle': 30.0, 'clockAngleRange': 60.0, 'interScanOverheadTime': 0.001, '@id': None, '@type': 'CONICAL'}, 
                            'targetBrightnessTemp': 295.0, '@id': 'ray3'})
        
        o = RadiometerModel.from_json(self.radio4_json)
        self.assertEqual(o.to_dict(), {'@type': 'Radiometer', 'name': None, 'mass': None, 'volume': None, 'power': None, 
                                       'orientation': {'referenceFrame': 'SC_BODY_FIXED', 'convention': 'EULER', 'eulerAngle1': 0.0, 'eulerAngle2': 330.0, 'eulerAngle3': 0.0, 'eulerSeq1': 1, 'eulerSeq2': 2, 'eulerSeq3': 3, '@id': None}, 
                                       'fieldOfViewGeometry': {'shape': 'CIRCULAR', 'diameter': 13.741474058602394, '@id': None}, 
                                       'sceneFieldOfViewGeometry': {'shape': 'CIRCULAR', 'diameter': 13.741474058602394, '@id': None}, 
                                       'maneuver': None, 'pointingOption': None, 'dataRate': None, 'bitsPerPixel': None, 
                                       'antenna': {'shape': 'CIRCULAR', 'apertureExcitationProfile': 'UNIFORM', 'diameter': 1.0, 'height': None, 'width': None, 'apertureEfficiency': None, 'radiationEfficiency': 1.0, 'phyTemp': 300.0, '@id': None}, 
                                       'operatingFrequency': 1250000000.0, 
                                       'system': {'tlLoss': 0.5, 'tlPhyTemp': 290.0, 'rfAmpGain': 30.0, 'rfAmpInpNoiseTemp': 200.0, 'rfAmpGainVariation': 10.0, 
                                                  'mixerGain,': 23.0, 'mixerInpNoiseTemp': 1200.0, 'mixerGainVariation': 2.0, 
                                                  'ifAmpGain': 30.0, 'ifAmpInputNoiseTemp': 100.0, 'ifAmpGainVariation': 10.0, 
                                                  'excessNoiseTemperature': 1000.0, 'integratorVoltageGain': 1.0, 'predetectionGain': None, 'predetectionInpNoiseTemp': None, 'predetectionGainVariation': None, 
                                                  'integrationTime': 1.0, 'bandwidth': 100000000.0, '@id': None, '@type': 'NOISE_ADDING'}, 
                                        'scan': {'@id': None, '@type': 'FIXED'}, 'targetBrightnessTemp': 295.0, '@id': 'ray4'})
    
    def test_calc_data_metrics_1(self):
        """ ``instru_look_angle_from_target_inc_angle`` flag is set to False."""

        epoch_JDUT1 =  2458543.06088 # 2019 Feb 28 13:27:40 is time at which the ECEF and ECI frames approximately align, hence ECEF to ECI rotation is identity. See <https://www.celnav.de/longterm.htm> online calculator of GMST.
        SpacecraftOrbitState = {'time [JDUT1]':epoch_JDUT1, 'x [km]': 6878.137, 'y [km]': 0, 'z [km]': 0, 'vx [km/s]': 0, 'vy [km/s]': 7.6126, 'vz [km/s]': 0} # altitude 500 km
        TargetCoords = {'lat [deg]': 0, 'lon [deg]': 0}

        o = RadiometerModel.from_json(self.radio1_json)
        data_metrics = o.calc_data_metrics(sc_orbit_state=SpacecraftOrbitState, target_coords=TargetCoords, instru_look_angle_from_target_inc_angle=False)
        
        self.assertEqual(data_metrics, {'ground pixel along-track resolution [m]': 119917.0, 'ground pixel cross-track resolution [m]': 119917.02, 
                                        'swath-width [m]': 120565.56, 'sensitivity [K]': 17.94, 'incidence angle [deg]': 0.03, 'beam efficiency': np.nan})

        o = RadiometerModel.from_json(self.radio2_json)
        TargetCoords = {'lat [deg]': 5, 'lon [deg]': 0} # target pixel somewhere on the cross-track direction.
        data_metrics = o.calc_data_metrics(sc_orbit_state=SpacecraftOrbitState, target_coords=TargetCoords, instru_look_angle_from_target_inc_angle=False)
        
        self.assertEqual(data_metrics, {'ground pixel along-track resolution [m]': 161269.89, 'ground pixel cross-track resolution [m]': 260072.11, 'swath-width [m]': 3159185.93, 
                                        'sensitivity [K]': 0.2, 'incidence angle [deg]': 51.68, 'beam efficiency': 0.24})

        o = RadiometerModel.from_json(self.radio3_json)
        TargetCoords = {'lat [deg]': 0, 'lon [deg]': 2.62}
        data_metrics = o.calc_data_metrics(sc_orbit_state=SpacecraftOrbitState, target_coords=TargetCoords, instru_look_angle_from_target_inc_angle=False)
        
        # calculated pixels resolutions are not accurate since the imaged pixel is not at side-looking geometry      
        self.assertEqual(data_metrics, {'ground pixel along-track resolution [m]': 40047.33, 'ground pixel cross-track resolution [m]': 47491.27, 
                                        'swath-width [m]': 306358.49, 'sensitivity [K]': 0.21, 'incidence angle [deg]': 32.51, 'beam efficiency': np.nan})
        
        o = RadiometerModel.from_json(self.radio4_json)
        TargetCoords = {'lat [deg]': -2.62, 'lon [deg]': 0}
        data_metrics = o.calc_data_metrics(sc_orbit_state=SpacecraftOrbitState, target_coords=TargetCoords, instru_look_angle_from_target_inc_angle=False)
        self.assertEqual(data_metrics, {'ground pixel along-track resolution [m]': 140198.56, 'ground pixel cross-track resolution [m]': 166301.75, 
                                        'swath-width [m]': 168745.48, 'sensitivity [K]': 0.23, 'incidence angle [deg]': 32.54, 'beam efficiency': np.nan})

                                
    def test_calc_data_metrics_2(self):
        """ ``instru_look_angle_from_target_inc_angle`` flag is set to True."""

        epoch_JDUT1 =  2458543.06088 # 2019 Feb 28 13:27:40 is time at which the ECEF and ECI frames approximately align, hence ECEF to ECI rotation is identity. See <https://www.celnav.de/longterm.htm> online calculator of GMST.
        SpacecraftOrbitState = {'time [JDUT1]':epoch_JDUT1, 'x [km]': 6878.137, 'y [km]': 0, 'z [km]': 0, 'vx [km/s]': 0, 'vy [km/s]': 7.6126, 'vz [km/s]': 0} # altitude 500 km
        TargetCoords = {'lat [deg]': 0, 'lon [deg]': 0.5}

        o = RadiometerModel.from_json(self.radio1_json)
        data_metrics = o.calc_data_metrics(sc_orbit_state=SpacecraftOrbitState, target_coords=TargetCoords, instru_look_angle_from_target_inc_angle=True)
        self.assertEqual(data_metrics, {'ground pixel along-track resolution [m]': 120708.29, 'ground pixel cross-track resolution [m]': 121567.92, 
                                        'swath-width [m]': 122255.0, 'sensitivity [K]': 17.94, 'incidence angle [deg]': 6.82, 'beam efficiency': np.nan})

        

