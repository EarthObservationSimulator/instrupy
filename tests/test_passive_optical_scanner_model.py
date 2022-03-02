"""Unit tests for instrupy.passive_optical_sensor_model

Tests:

* test_from_json_basic
* test_calculate_integration_time
* test_calculate_number_of_signal_electrons
* test_radianceWithEarthAsReflector
* test_radianceWithEarthAsBlackBodyRadiator
* test_planck_photon_integral
* test_calc_data_metrics_smad_truth: Test with SMAD 3rd ed truth data (firesat example).

For the below tests refer to the following article more context:
V. Ravindra and S. Nag, "Instrument Data Metrics Evaluator for Tradespace Analysis of Earth Observing Constellations," 2020 IEEE Aerospace Conference, Big Sky, MT, USA, 2020.

* test_calc_data_metrics_TIRSBand1_precomputed: Model instrument with TIRS Band 1 specs (IR, pushbroom), Landsat-8 orbit and test with the results as computed on 4 April 2021.
* test_calc_data_metrics_OLIBlueBand_precomputed: Model instrument with OLI Blue band specs (Optical, pushbroom), Landsat-8 orbit and test with the results as computed on 4 April 2021.
* test_calc_data_metrics_MODISBand10_precomputed: Model instrument with MODIS Band 10 specs (Optical, whiskbroom), Aqua orbit and test with the results as computed on 4 April 2021.
* test_calc_data_metrics_MODISBand1_precomputed: Model instrument with MODIS Band 1 specs (Optical, whiskbroom), Aqua orbit and test with the results as computed on 4 April 2021.
* test_calc_data_metrics_CCAMBlueBand_precomputed: Model instrument with CCAM Blue Band specs (Matrix, optical), Aqua orbit and test with the results as computed on 4 April 2021.

.. note:: The results using the LOWTRAN7 model are computed at resolution wav_step_percm = 5.

"""

import unittest
import json
import numpy as np
import sys, os
import math

from instrupy.passive_optical_scanner_model import ScanTech, PassiveOpticalScannerModel, AtmosphericLossModel
from instrupy.util import Orientation, SphericalGeometry, ViewGeometry, Maneuver

firesat_json =  '{"@type": "Passive Optical Scanner",' \
                '"name": "FireSat",' \
                '"mass": 28,' \
                '"volume": 0.12,' \
                '"power": 32,' \
                '"fieldOfViewGeometry": {' \
                '   "shape": "RECTanGULAR",' \
                '   "angleHeight": 0.628,' \
                '   "angleWidth": 115.8' \
                ' },' \
                '"scanTechnique": "WhiskBROOM",' \
                '"orientation": {' \
                '   "referenceFrame": "SC_BODY_FIXED",' \
                '   "convention": "SIDE_loOK",' \
                '   "sideLookAngle": 0' \
                ' },' \
                '"dataRate": 85,' \
                '"numberDetectorRows": 256,' \
                '"numberDetectorCols": 1,' \
                '"detectorWidth": 30e-6,' \
                '"focalLength": 0.7,' \
                '"operatingWavelength": 4.2e-6,' \
                '"bandwidth": 1.9e-6,' \
                '"quantumEff": 0.5,' \
                '"targetBlackBodyTemp": 290,' \
                '"bitsPerPixel": 8,' \
                '"opticsSysEff": 0.75,' \
                '"numOfReadOutE": 25,' \
                '"apertureDia": 0.26,' \
                '"Fnum": 2.7,' \
                '"atmosLossModel": "LOWTRAN7"}'
class TestPassiveOpticalScannerModel(unittest.TestCase):   

    def assertNearlyZeroErrorFraction(self,a,b,fraction=0.01,msg=None):
        if abs(a-b) > abs(fraction*a):
            if msg is None:
                self.fail("The given numbers %s and %s are not near each other."%(a,b))
            else:
                self.fail(msg)

    def test_from_json_basic(self):        
        # Test: Typical case      
        firesat = PassiveOpticalScannerModel.from_json(firesat_json)
        self.assertEqual(firesat._type, "Passive Optical Scanner")
        self.assertIsNotNone(firesat._id) # default random-id assigned
        self.assertIsInstance(firesat.name, str)
        self.assertEqual(firesat.name, "FireSat")
        self.assertIsInstance(firesat.mass, float)
        self.assertEqual(firesat.mass, 28)
        self.assertIsInstance(firesat.volume, float)
        self.assertEqual(firesat.volume, 0.12)
        self.assertIsInstance(firesat.power, float)
        self.assertEqual(firesat.power, 32)
        self.assertIsInstance(firesat.dataRate, float)
        self.assertEqual(firesat.dataRate, 85)        
        self.assertIsInstance(firesat.orientation, Orientation)
        self.assertEqual(firesat.orientation, Orientation.from_dict({'referenceFrame': 'SC_BODY_FIXED', 'convention': 'EULER', 'eulerSeq1':1, 'eulerSeq2':2, 'eulerSeq3':3, 'eulerAngle1':0, 'eulerAngle2':0, 'eulerAngle3':0}))
        self.assertIsInstance(firesat.fieldOfView, ViewGeometry)
        self.assertEqual(firesat.fieldOfView, ViewGeometry.from_dict({'orientation': {'referenceFrame': 'SC_BODY_FIXED', 'convention': 'EULER', 'eulerSeq1':1, 'eulerSeq2':2, 'eulerSeq3':3, 'eulerAngle1':0, 'eulerAngle2':0, 'eulerAngle3':0}, 
                                                                           'sphericalGeometry':{'shape':'rectangular', 'angleHeight':0.628, 'angleWidth':115.8}}))
        self.assertIsInstance(firesat.sceneFieldOfView, ViewGeometry)
        self.assertEqual(firesat.sceneFieldOfView, firesat.fieldOfView)
        self.assertIsNone(firesat.maneuver)
        self.assertIsNone(firesat.fieldOfRegard)
        self.assertIsNone(firesat.pointingOption)
        self.assertIsInstance(firesat.scanTechnique, str)
        self.assertEqual(firesat.scanTechnique, "WHISKBROOM")
        self.assertIsInstance(firesat.detectorWidth, float)        
        self.assertEqual(firesat.detectorWidth, 30e-6)
        self.assertIsInstance(firesat.focalLength, float)
        self.assertEqual(firesat.focalLength, 0.7)
        self.assertIsInstance(firesat.operatingWavelength, float)
        self.assertEqual(firesat.operatingWavelength, 4.2e-6)
        self.assertIsInstance(firesat.bandwidth, float)
        self.assertEqual(firesat.bandwidth, 1.9e-6)
        self.assertIsInstance(firesat.quantumEff, float)
        self.assertEqual(firesat.quantumEff, 0.5)
        self.assertIsInstance(firesat.targetBlackBodyTemp, float)
        self.assertEqual(firesat.targetBlackBodyTemp, 290)
        self.assertIsInstance(firesat.bitsPerPixel, int)
        self.assertEqual(firesat.bitsPerPixel, 8)
        self.assertIsInstance(firesat.opticsSysEff, float)
        self.assertEqual(firesat.opticsSysEff, 0.75)
        self.assertIsInstance(firesat.numOfReadOutE, float)
        self.assertEqual(firesat.numOfReadOutE, 25)
        self.assertIsInstance(firesat.apertureDia, float)
        self.assertEqual(firesat.apertureDia, 0.26)
        self.assertIsInstance(firesat.Fnum, float)
        self.assertEqual(firesat.Fnum,  2.7)
        self.assertIsInstance(firesat.atmosLossModel, AtmosphericLossModel)
        self.assertTrue(firesat.atmosLossModel, AtmosphericLossModel.LOWTRAN7)        
        self.assertIsInstance(firesat, PassiveOpticalScannerModel)

        # Test the setting of default parameters
        o = PassiveOpticalScannerModel.from_json('{"@type": "Passive Optical Scanner",'
                                                    '"fieldOfViewGeometry": {'
                                                    '   "shape": "RECTanGULAR",'
                                                    '   "angleHeight": 0.628,'
                                                    '   "angleWidth": 115.8'
                                                    ' },'
                                                    '"scanTechnique": "WhiskBROOM",'
                                                    '"numberDetectorRows": 256,'
                                                    '"numberDetectorCols": 1,'
                                                    '"detectorWidth": 30e-6,'
                                                    '"focalLength": 0.7,'
                                                    '"operatingWavelength": 4.2e-6,'
                                                    '"bandwidth": 1.9e-6,'
                                                    '"quantumEff": 0.5,'
                                                    '"opticsSysEff": 0.75,'
                                                    '"numOfReadOutE": 25,'
                                                    '"apertureDia": 0.26,'
                                                    '"Fnum": 2.7}')

        self.assertIsNotNone(o._id)
        self.assertEqual(o.orientation, Orientation.from_dict({'referenceFrame': 'SC_BODY_FIXED', 'convention': 'EULER', 'eulerSeq1':1, 'eulerSeq2':2, 'eulerSeq3':3, 'eulerAngle1':0, 'eulerAngle2':0, 'eulerAngle3':0}))
        self.assertEqual(o.fieldOfView, o.sceneFieldOfView)
        self.assertEqual(firesat.targetBlackBodyTemp, 290)
        self.assertIsNone(o.atmosLossModel)

        # Test of an un-supported field-of-view specification. 
        with self.assertRaises(Exception):
            o = PassiveOpticalScannerModel.from_json('{"@type": "Passive Optical Scanner",'
                                                       '"name": "FireSat",'
                                                        '"mass": 28,'
                                                        '"volume": 0.12,'
                                                        '"power": 32,'
                                                        '"fieldOfViewGeometry": {'
                                                        '   "shape": "CIRCULAR",'
                                                        '   "diameter": 60,'
                                                        ' },'
                                                        '"scanTechnique": "WhiskBROOM",'
                                                        '"orientation": {'
                                                        '   "referenceFrame": "SC_BODY_FIXED",'
                                                        '   "convention": "SIDE_loOK",'
                                                        '   "sideLookAngle": 0'
                                                        ' },'
                                                        '"dataRate": 85,'
                                                        '"numberDetectorRows": 256,'
                                                        '"numberDetectorCols": 1,'
                                                        '"detectorWidth": 30e-6,'
                                                        '"focalLength": 0.7,'
                                                        '"operatingWavelength": 4.2e-6,'
                                                        '"bandwidth": 1.9e-6,'
                                                        '"quantumEff": 0.5,'
                                                        '"targetBlackBodyTemp": 290,'
                                                        '"bitsPerPixel": 8,'
                                                        '"opticsSysEff": 0.75,'
                                                        '"numOfReadOutE": 25,'
                                                        '"apertureDia": 0.26,'
                                                        '"Fnum": 2.7,'
                                                        '"atmosLossModel": "LOWTRAN7"}')

        # Test of an improper scanning technique specification
        with self.assertRaises(Exception):
            o = PassiveOpticalScannerModel.from_json('{"@type": "Passive Optical Scanner",'
                                                       '"name": "FireSat",'
                                                        '"mass": 28,'
                                                        '"volume": 0.12,'
                                                        '"power": 32,'
                                                        '"fieldOfViewGeometry": {'
                                                        '   "shape": "RECTanGULAR",'
                                                        '   "angleHeight": 0.628,'
                                                        '   "angleWidth": 115.8'
                                                        ' },'
                                                        '"scanTechnique": "Broombroom",'
                                                        '"orientation": {'
                                                        '   "referenceFrame": "SC_BODY_FIXED",'
                                                        '   "convention": "SIDE_loOK",'
                                                        '   "sideLookAngle": 0'
                                                        ' },'
                                                        '"dataRate": 85,'
                                                        '"numberDetectorRows": 256,'
                                                        '"numberDetectorCols": 1,'
                                                        '"detectorWidth": 30e-6,'
                                                        '"focalLength": 0.7,'
                                                        '"operatingWavelength": 4.2e-6,'
                                                        '"bandwidth": 1.9e-6,'
                                                        '"quantumEff": 0.5,'
                                                        '"targetBlackBodyTemp": 290,'
                                                        '"bitsPerPixel": 8,'
                                                        '"opticsSysEff": 0.75,'
                                                        '"numOfReadOutE": 25,'
                                                        '"apertureDia": 0.26,'
                                                        '"Fnum": 2.7,'
                                                        '"atmosLossModel": "LOWTRAN7"}')


        # Test of an PUSHBROOM scanning technique specification and more than one :code:`numberDetectorRows` specification.
        with self.assertRaises(Exception):
            o = PassiveOpticalScannerModel.from_json('{"@type": "Passive Optical Scanner",'
                                                       '"name": "FireSat",'
                                                        '"mass": 28,'
                                                        '"volume": 0.12,'
                                                        '"power": 32,'
                                                        '"fieldOfViewGeometry": {'
                                                        '   "shape": "RECTanGULAR",'
                                                        '   "angleHeight": 0.628,'
                                                        '   "angleWidth": 115.8'
                                                        ' },'
                                                        '"scanTechnique": "Pushbroom",'
                                                        '"orientation": {'
                                                        '   "referenceFrame": "SC_BODY_FIXED",'
                                                        '   "convention": "SIDE_loOK",'
                                                        '   "sideLookAngle": 0'
                                                        ' },'
                                                        '"dataRate": 85,'
                                                        '"numberDetectorRows": 10,'
                                                        '"numberDetectorCols": 10,'
                                                        '"detectorWidth": 30e-6,'
                                                        '"focalLength": 0.7,'
                                                        '"operatingWavelength": 4.2e-6,'
                                                        '"bandwidth": 1.9e-6,'
                                                        '"quantumEff": 0.5,'
                                                        '"targetBlackBodyTemp": 290,'
                                                        '"bitsPerPixel": 8,'
                                                        '"opticsSysEff": 0.75,'
                                                        '"numOfReadOutE": 25,'
                                                        '"apertureDia": 0.26,'
                                                        '"Fnum": 2.7,'
                                                        '"atmosLossModel": "LOWTRAN7"}')

        # Test of an WHISKBROOM scanning technique specification and more than one :code:`numberDetectorCols` specification.
        with self.assertRaises(Exception):
            o = PassiveOpticalScannerModel.from_json('{"@type": "Passive Optical Scanner",'
                                                       '"name": "FireSat",'
                                                        '"mass": 28,'
                                                        '"volume": 0.12,'
                                                        '"power": 32,'
                                                        '"fieldOfViewGeometry": {'
                                                        '   "shape": "RECTanGULAR",'
                                                        '   "angleHeight": 0.628,'
                                                        '   "angleWidth": 115.8'
                                                        ' },'
                                                        '"scanTechnique": "WhiskBROOM",'
                                                        '"orientation": {'
                                                        '   "referenceFrame": "SC_BODY_FIXED",'
                                                        '   "convention": "SIDE_loOK",'
                                                        '   "sideLookAngle": 0'
                                                        ' },'
                                                        '"dataRate": 85,'
                                                        '"numberDetectorRows": 10,'
                                                        '"numberDetectorCols": 10,'
                                                        '"detectorWidth": 30e-6,'
                                                        '"focalLength": 0.7,'
                                                        '"operatingWavelength": 4.2e-6,'
                                                        '"bandwidth": 1.9e-6,'
                                                        '"quantumEff": 0.5,'
                                                        '"targetBlackBodyTemp": 290,'
                                                        '"bitsPerPixel": 8,'
                                                        '"opticsSysEff": 0.75,'
                                                        '"numOfReadOutE": 25,'
                                                        '"apertureDia": 0.26,'
                                                        '"Fnum": 2.7,'
                                                        '"atmosLossModel": "LOWTRAN7"}')
    
    def test_planck_photon_integral(self):
        # Test trivial case with 0 wavelength
        self.assertAlmostEqual(PassiveOpticalScannerModel.planck_photon_integral(0,290), 0)

        """ Tests using online calculator from <https://www.opticsthewebsite.com/OpticsCalculators.aspx> as truth data.
            Note that the online calculator requires minimum wavelength to be set as 1e-9 um which is nearly 0 wavelength.

        """
        self.assertNearlyZeroErrorFraction(PassiveOpticalScannerModel.planck_photon_integral(12e-6,1500), 1.46801e+20 * 1e4)
        self.assertNearlyZeroErrorFraction(PassiveOpticalScannerModel.planck_photon_integral(2e-6,500), 3.36875e+14 * 1e4) 
        self.assertNearlyZeroErrorFraction(PassiveOpticalScannerModel.planck_photon_integral(500e-6,45), 4.10754e+15 * 1e4) 
        self.assertNearlyZeroErrorFraction(PassiveOpticalScannerModel.planck_photon_integral(1e9,45), 4.40891e+15 * 1e4)  # specifying 1e9 m as wavelength is to get approximately the radiance over entire spectrum
    
    def test_radianceWithEarthAsBlackBodyRadiator(self):
        """ Tests using online calculator from <https://www.opticsthewebsite.com/OpticsCalculators.aspx> as truth data for in-band radiance calculation.
            Note that the online calculator requires minimum wavelength to be set as 1e-9 um which is nearly 0 wavelength.

        """
        self.assertNearlyZeroErrorFraction(PassiveOpticalScannerModel.radiance_with_earth_as_bb_radiator(16.25e-6, 12.5e-6, 290, 0, atmos_loss_model=None),  5.08113e+17 * 1e4) # 10 um to 22.5 um at 290 K, 0 deg incidence angle
        self.assertNearlyZeroErrorFraction(PassiveOpticalScannerModel.radiance_with_earth_as_bb_radiator(6e-6, 4e-6, 180.5, 0, atmos_loss_model=None),  6.74299e+14 * 1e4) 
        self.assertNearlyZeroErrorFraction(PassiveOpticalScannerModel.radiance_with_earth_as_bb_radiator(0.5e-6, 0.8e-6, 270, 0, atmos_loss_model=None),  2.74990e-5 * 1e4) 
        self.assertNearlyZeroErrorFraction(PassiveOpticalScannerModel.radiance_with_earth_as_bb_radiator(140e-6, 40e-6, 330, 0, atmos_loss_model=None),  1.77269e+16 * 1e4)

        # Tests with 10 um to 22.5 um at 290 K at different incidence angles
        self.assertNearlyZeroErrorFraction(PassiveOpticalScannerModel.radiance_with_earth_as_bb_radiator(16.25e-6, 12.5e-6, 290, 0.1, atmos_loss_model=None),  5.08113e+17 * 1e4 * np.cos(0.1)) 
        self.assertNearlyZeroErrorFraction(PassiveOpticalScannerModel.radiance_with_earth_as_bb_radiator(16.25e-6, 12.5e-6, 290, -0.1, atmos_loss_model=None),  5.08113e+17 * 1e4 * np.cos(0.1))
        self.assertNearlyZeroErrorFraction(PassiveOpticalScannerModel.radiance_with_earth_as_bb_radiator(16.25e-6, 12.5e-6, 290, 0.5 + 2*3.141, atmos_loss_model=None),  5.08113e+17 * 1e4 * np.cos(0.5))
        self.assertNearlyZeroErrorFraction(PassiveOpticalScannerModel.radiance_with_earth_as_bb_radiator(16.25e-6, 12.5e-6, 290, -0.5 - 8*3.141, atmos_loss_model=None),  5.08113e+17 * 1e4 * np.cos(0.5))

        # Tests with unrealistic observation incidence angles
        with self.assertRaises(Exception):
            PassiveOpticalScannerModel.radiance_with_earth_as_bb_radiator(16.25e-6, 12.5e-6, 290, np.deg2rad(91), atmos_loss_model=None)
        with self.assertRaises(Exception):
            PassiveOpticalScannerModel.radiance_with_earth_as_bb_radiator(16.25e-6, 12.5e-6, 290, np.deg2rad(-91), atmos_loss_model=None)
        with self.assertRaises(Exception):
            PassiveOpticalScannerModel.radiance_with_earth_as_bb_radiator(16.25e-6, 12.5e-6, 290, np.deg2rad(-91 + 360*8), atmos_loss_model=None)
    
    def test_radianceWithEarthAsReflector(self):        
        """ Initialize parameters which ensure that the target (ground-pixel) and observer (satelltie) and the Sun have a LOS geometry. 
            This is verified using GMAT (Orbit View animation).

        """
        tObs_JDUT1 = 2451623.999630 # spring equinox day
        obs_pos_km = [6577.848345501363, -9.521529479781905e-013, 2394.141003279681] # satelltie is in the XZ plane
        tar_pos_km = [6577.848345501363 - 622, -9.521529479781905e-013, 2394.141003279681] # satellite altitude is 622 km and target is in the XZ plane
        obs_area_m2 = 1

        """ Test: Reflected energy far-outside visible wavelengths must be near 0. 
        """
        self.assertAlmostEqual(PassiveOpticalScannerModel.radiance_with_earth_as_reflector(0.5e-9, 0.2e-9, tObs_JDUT1, obs_pos_km, tar_pos_km, obs_area_m2, atmos_loss_model='LOWTRAN7'), 0)
        self.assertAlmostEqual(PassiveOpticalScannerModel.radiance_with_earth_as_reflector(1, 1e-2, tObs_JDUT1, obs_pos_km, tar_pos_km, obs_area_m2, atmos_loss_model='LOWTRAN7'), 0, places = 3)

        """ Test: Reflected energy for visible wavelengths must be greater than that of other wavelengths, keeping bandwidth same.
        """        
        self.assertGreater(PassiveOpticalScannerModel.radiance_with_earth_as_reflector(0.5e-6, 0.2e-6, tObs_JDUT1, obs_pos_km, tar_pos_km, obs_area_m2, atmos_loss_model='LOWTRAN7'), PassiveOpticalScannerModel.radiance_with_earth_as_reflector(6e-6, 0.2e-6, tObs_JDUT1, obs_pos_km, tar_pos_km, obs_area_m2, atmos_loss_model='LOWTRAN7')) # longer wavelengths than visible
        self.assertGreater(PassiveOpticalScannerModel.radiance_with_earth_as_reflector(0.5e-6, 0.2e-6, tObs_JDUT1, obs_pos_km, tar_pos_km, obs_area_m2, atmos_loss_model='LOWTRAN7'), PassiveOpticalScannerModel.radiance_with_earth_as_reflector(0.25e-6, 0.2e-6, tObs_JDUT1, obs_pos_km, tar_pos_km, obs_area_m2, atmos_loss_model='LOWTRAN7')) # shorter wavelengths than visible

        """
            Test with another observer position, one where the observer sees the target pixel (of fixed area) at a larger angle.
        """
        opWav_m = 0.5e-6
        bw_m = 0.2e-6
        obs2_pos_km = [6893.654271085462, -9.186593534864809e-013, 1215.537243668513]
        self.assertGreater(PassiveOpticalScannerModel.radiance_with_earth_as_reflector(opWav_m, bw_m, tObs_JDUT1, obs_pos_km, tar_pos_km, obs_area_m2, atmos_loss_model='LOWTRAN7'), PassiveOpticalScannerModel.radiance_with_earth_as_reflector(opWav_m, bw_m, tObs_JDUT1, obs2_pos_km, tar_pos_km, obs_area_m2, atmos_loss_model='LOWTRAN7')) # longer wavelengths than visible
    
    def test_calculate_number_of_signal_electrons(self):        
        # Test: Truth data from SMAD 3rd edition, Table 9-15.
        opWav_m = 4.2e-6
        bw_m = 1.9e-6
        bbT_K = 290
        apDia_m = 0.2626
        opTrns = 0.75
        QE = 0.5   
        tObs_JDUT1 = 2451623.999630 #  date corresponds to satellite over target at day-time.
        obs_pos_km = [6378+700, 0, 0] # satellite on the X-axis at altitude 700 km
        tar_pos_km = [6378, 0, 0] 
        pixelArea_m2 = 30.0519 * 30.0519
        Ti_s = 24.1827e-6
        # The InstruPy computed value must be greater than truth value since SMAD does not consider the energy reflected off Sun, and the date corresponds to satellite over target at day-time.
        self.assertGreater(PassiveOpticalScannerModel.calculate_number_of_signal_electrons(opWav_m, bw_m, bbT_K, apDia_m, opTrns, QE, tObs_JDUT1, obs_pos_km, tar_pos_km, pixelArea_m2, Ti_s, atmos_loss_model='LOWTRAN7'), 8286.104444633884)        
    
    def test_calculate_integration_time(self):
        # Test: PUSHBROOM scanning
        self.assertAlmostEqual(PassiveOpticalScannerModel.calculate_integration_time(ScanTech.PUSHBROOM, 1, 1, 12.5, 1), 12.5)
        self.assertAlmostEqual(PassiveOpticalScannerModel.calculate_integration_time(ScanTech.PUSHBROOM, 1, 150, 12.5, 0.1), 12.5)
        self.assertAlmostEqual(PassiveOpticalScannerModel.calculate_integration_time(ScanTech.PUSHBROOM, 1, 150, 12.5, 0.1, angle_width_deg=30), 12.5)
        # check max exposure time functionality
        self.assertAlmostEqual(PassiveOpticalScannerModel.calculate_integration_time(ScanTech.PUSHBROOM, 1, 150, 12.5, 0.1, max_det_exp_time=5), 5)
        
        with self.assertRaises(Exception): # Exception expected if number of detector rows in along-track direction is not 1.
            PassiveOpticalScannerModel.calculate_integration_time(ScanTech.PUSHBROOM, 10, 150, 12.5, 1)

        # Test: Whiskroom scanning
        self.assertAlmostEqual(PassiveOpticalScannerModel.calculate_integration_time(ScanTech.WHISKBROOM, 1, 1, 12.5, 1, max_det_exp_time=None, angle_width_deg=30), 12.5/30)
        self.assertAlmostEqual(PassiveOpticalScannerModel.calculate_integration_time(ScanTech.WHISKBROOM, 1, 1, 12.5, 0.1, max_det_exp_time=None, angle_width_deg=30), 12.5/300)
        self.assertAlmostEqual(PassiveOpticalScannerModel.calculate_integration_time(ScanTech.WHISKBROOM, 20, 1, 12.5, 0.1, max_det_exp_time=None, angle_width_deg=30), 12.5/300)
        # check max exposure time functionality
        self.assertAlmostEqual(PassiveOpticalScannerModel.calculate_integration_time(ScanTech.WHISKBROOM, 20, 1, 12.5, 0.1, max_det_exp_time=(12.5/3000), angle_width_deg=30), 12.5/3000)

        with self.assertRaises(Exception): # Exception expected if number of detector columns in cross-track direction is not 1.
            PassiveOpticalScannerModel.calculate_integration_time(ScanTech.WHISKBROOM, 10, 150, 12.5, 1, max_det_exp_time=None, angle_width_deg=30)

        with self.assertRaises(Exception): # Exception expected if cross-track-fov is not specified
            PassiveOpticalScannerModel.calculate_integration_time(ScanTech.WHISKBROOM, 20, 1, 12.5, 0.1, max_det_exp_time=None, angle_width_deg=None)

        # Test: MATRIX_IMAGER scanning
        self.assertAlmostEqual(PassiveOpticalScannerModel.calculate_integration_time(ScanTech.MATRIX_IMAGER, 1, 1, 12.5, 1), 12.5)
        self.assertAlmostEqual(PassiveOpticalScannerModel.calculate_integration_time(ScanTech.MATRIX_IMAGER, 1, 1, 12.5, 0.1), 12.5)
        self.assertAlmostEqual(PassiveOpticalScannerModel.calculate_integration_time(ScanTech.MATRIX_IMAGER, 20, 1, 12.5, 0.1), 12.5)
        self.assertAlmostEqual(PassiveOpticalScannerModel.calculate_integration_time(ScanTech.MATRIX_IMAGER, 1, 20, 12.5, 0.1), 12.5)
        self.assertAlmostEqual(PassiveOpticalScannerModel.calculate_integration_time(ScanTech.MATRIX_IMAGER, 20, 20, 12.5, 0.1, 30), 12.5)
        # check max exposure time functionality
        self.assertAlmostEqual(PassiveOpticalScannerModel.calculate_integration_time(ScanTech.MATRIX_IMAGER, 1, 20, 12.5, 0.1, max_det_exp_time=1.2), 1.2)
        
    def test_calc_data_metrics_smad_truth(self):        
        """ Test: Truth data from SMAD 3rd edition, Table 9-15. Note that the  inputs are made so that they are consistent with the example.
             Further note that SMAD does not consider the energy from Sun reflected off Sun in their calculations.
        """
        firesat = PassiveOpticalScannerModel.from_json(firesat_json)
        epoch_JDUT1 = 2451623.999630
        sc_orbit_state = {'time [JDUT1]':epoch_JDUT1, 'x [km]': 7078.137, 'y [km]': 0, 'z [km]': 0, 'vx [km/s]': 0, 'vy [km/s]': 7.5, 'vz [km/s]': 0} # equatorial orbit, altitude about 700 km
        target_coords = {'lat [deg]': 0, 'lon [deg]': 0} # lat = 0, lon = 0 corresponds to [6378, 0, 0] km in ECI for observer position, check using Matlab function: eci2lla([6378, 0, 0] ,[2000 3 20 11 59 28.000])
        obsv_metrics = firesat.calc_data_metrics(sc_orbit_state, target_coords)
        self.assertAlmostEqual(obsv_metrics["ground pixel along-track resolution [m]"], 30, delta=3)
        self.assertAlmostEqual(obsv_metrics["ground pixel cross-track resolution [m]"], 30, delta=3)
        # A (positive) deviation is expected since SMAD does not consider the energy from Sun reflected off Sun
        self.assertGreater(obsv_metrics["SNR"], 88)
        self.assertGreater(obsv_metrics["dynamic range"], 332.9)
        self.assertAlmostEqual(obsv_metrics["noise-equivalent delta T [K]"], 0.3 , delta=0.01)
        self.assertEqual(obsv_metrics, {'ground pixel along-track resolution [m]': 31.34, 'ground pixel cross-track resolution [m]': 32.91, 'SNR': 148.15, 'dynamic range': 902.27, 'noise-equivalent delta T [K]': 0.30609})


        # @TODO, make test with satellite in night region, thus no reflected energy of Sun. Result should match SMAD.
    
    def test_calc_data_metrics_TIRSBand1_precomputed(self):
        """ Model instrument with TIRS Band 1 specs (IR, pushbroom) and test with the results as computed on 4 April 2021.

        Refer to the following article more context:
        V. Ravindra and S. Nag, "Instrument Data Metrics Evaluator for Tradespace Analysis of Earth Observing Constellations," 
        2020 IEEE Aerospace Conference, Big Sky, MT, USA, 2020.

        .. note:: The results using the LOWTRAN7 model are computed at resolution wav_step_percm = 5.

        """
        landsat_tirs_band1_dict = {
                                    "@type": "Passive Optical Scanner",
                                    "name": "Landsat 8 TIRS Band1",
                                    "mass": 236,
                                    "volume": 0.261, 
                                    "power": 380, 
                                    "fieldOfViewGeometry": {
                                        "shape": "RECTANGULAR",
                                        "angleHeight": 0.0081,
                                        "angleWidth": 15
                                    },
                                    "scanTechnique": "PUSHBROOM",
                                    "orientation": {
                                        "referenceFrame": "SC_BODY_FIXED",
                                        "convention": "REF_FRAME_ALIGNED"
                                    },
                                    "dataRate":  384,
                                    "numberDetectorRows": 1,
                                    "numberDetectorCols": 1850,
                                    "detectorWidth": 25e-6,
                                    "focalLength": 0.178,
                                    "operatingWavelength": 10.9e-6,
                                    "bandwidth": 0.6e-6,
                                    "quantumEff": 0.025,
                                    "targetBlackBodyTemp": 290,
                                    "bitsPerPixel": 12,
                                    "opticsSysEff": 0.60 ,
                                    "numOfReadOutE":  20,
                                    "apertureDia":  0.1085366,
                                    "Fnum":  1.64,
                                    "maxDetectorExposureTime": 3.49e-3,
                                    "atmosLossModel": "LOWTRAN7",
                                    "_comments": ["Above is Total payload data-rate not just off the TIRS.",
                                                "numReadOutE is guessed."]
                                   }
        landsat_tirs_band1 = PassiveOpticalScannerModel.from_dict(landsat_tirs_band1_dict)
        # landsat 8 orbit at 10 Apr 2021 14:24:17.819 UTC            
        sc_orbit_state = {'time [JDUT1]':2459315.100208333,  'x [km]': -7012.215259847972,    'y [km]': 981.6284579029395,    'z [km]': 16.62328546479549, 
                                                            'vx [km/s]': 0.1664588472531363, 'vy [km/s]': 1.055747095699285, 'vz [km/s]': 7.426472416008381 }
        target_coords = {'lat [deg]': 0.01942147899019397 , 'lon [deg]': 117.1899962481559} # nadir position of satellite
        obsv_metrics = landsat_tirs_band1.calc_data_metrics(sc_orbit_state, target_coords)
        self.assertEqual(obsv_metrics, {'ground pixel along-track resolution [m]': 98.78, 'ground pixel cross-track resolution [m]': 98.92, 'SNR': 1507.48, 'dynamic range': 113645.23, 'noise-equivalent delta T [K]': 0.04162})
        # disable LOWTRAN7 atmospheric loss model and evaluate results
        landsat_tirs_band1_dict["atmosLossModel"] = None
        landsat_tirs_band1 = PassiveOpticalScannerModel.from_dict(landsat_tirs_band1_dict)
        obsv_metrics = landsat_tirs_band1.calc_data_metrics(sc_orbit_state, target_coords)
        self.assertEqual(obsv_metrics, {'ground pixel along-track resolution [m]': 98.78, 'ground pixel cross-track resolution [m]': 98.92, 'SNR': 1610.69, 'dynamic range': 129736.13, 'noise-equivalent delta T [K]': 0.03895})
 
    def test_calc_data_metrics_OLIBlueBand_precomputed(self):
        """  Model instrument with OLI Blue band specs (Optical, pushbroom) and test with the results as computed on 4 April 2021.

        Refer to the following article more context:
        V. Ravindra and S. Nag, "Instrument Data Metrics Evaluator for Tradespace Analysis of Earth Observing Constellations," 
        2020 IEEE Aerospace Conference, Big Sky, MT, USA, 2020.

        .. note:: The results using the LOWTRAN7 model are computed at resolution wav_step_percm = 5.

        """
        landsat_oli_blue_dict = {   "@type": "Passive Optical Scanner",
                                    "name": "Landsat 8 OLI Blue band",
                                    "mass": 1,
                                    "volume": 1, 
                                    "power": 1, 
                                    "fieldOfViewGeometry": {
                                        "shape": "RECTANGULAR",
                                        "angleHeight": 0.00244080020725731,
                                        "angleWidth": 15
                                    },
                                    "scanTechnique": "PUSHBROOM",
                                    "dataRate":  384,
                                    "numberDetectorRows": 1,
                                    "numberDetectorCols": 6146,
                                    "detectorWidth": 36e-6,
                                    "focalLength": 845.1e-3,
                                    "operatingWavelength": 482e-9,
                                    "bandwidth": 65e-9,
                                    "quantumEff":  0.85,
                                    "targetBlackBodyTemp": 290,
                                    "bitsPerPixel": 12,
                                    "opticsSysEff": 0.90 ,
                                    "numOfReadOutE":  8,
                                    "apertureDia": 0.1320,
                                    "Fnum":  6.4,
                                    "maxDetectorExposureTime": 3.6e-3,
                                    "atmosLossModel": "LOWTRAN7",
                                    "_comments": ["Above is Total payload data-rate not just off the OLI.",
                                                  "Mass, power and volume are simply wrong.", 
                                                  "numReadOutE is guessed."]
                                }
        landsat_oli_blue = PassiveOpticalScannerModel.from_dict(landsat_oli_blue_dict)
        # landsat 8 orbit at 10 Apr 2021 14:24:17.819 UTC  (NIGHT time)          
        sc_orbit_state = {'time [JDUT1]':2459315.100208333,  'x [km]': -7012.215259847972,    'y [km]': 981.6284579029395,    'z [km]': 16.62328546479549, 
                                                            'vx [km/s]': 0.1664588472531363, 'vy [km/s]': 1.055747095699285, 'vz [km/s]': 7.426472416008381 }
        target_coords = {'lat [deg]': 0.01942147899019397 , 'lon [deg]': 117.1899962481559} # nadir position of satellite
        obsv_metrics = landsat_oli_blue.calc_data_metrics(sc_orbit_state, target_coords)
        self.assertEqual(obsv_metrics['ground pixel along-track resolution [m]'], 29.96)
        self.assertEqual(obsv_metrics['ground pixel cross-track resolution [m]'], 30.0)
        self.assertEqual(obsv_metrics['SNR'], 0.0)
        self.assertEqual(obsv_metrics['dynamic range'], 0.0)
        assert math.isclose(obsv_metrics['noise-equivalent delta T [K]'], 1318917697165785.8, rel_tol=1e-9)

        # 10 Apr 2021 15:07:22.788   (Day time)                                    
        sc_orbit_state = {'time [JDUT1]':2459315.130127315,  'x [km]': 6512.435033854175,    'y [km]': -511.354859668807,    'z [km]': 2713.225164499847, 
                                                            'vx [km/s]': 2.748229230268995, 'vy [km/s]': -1.377488737466059, 'vz [km/s]': -6.850883979753837 }
        target_coords = {'lat [deg]': 22.7948561533915 , 'lon [deg]': -70.13495405345812} # nadir position of satellite
        obsv_metrics = landsat_oli_blue.calc_data_metrics(sc_orbit_state, target_coords)
        self.assertEqual(obsv_metrics, {'ground pixel along-track resolution [m]': 29.67, 'ground pixel cross-track resolution [m]': 29.73, 'SNR': 1065.8, 'dynamic range': 141998.68, 'noise-equivalent delta T [K]': np.inf})
        
        # disable LOWTRAN7 atmospheric loss model and evaluate results
        landsat_oli_blue_dict["atmosLossModel"] = None
        landsat_oli_blue = PassiveOpticalScannerModel.from_dict(landsat_oli_blue_dict)
        obsv_metrics = landsat_oli_blue.calc_data_metrics(sc_orbit_state, target_coords)
        self.assertEqual(obsv_metrics, {'ground pixel along-track resolution [m]': 29.67, 'ground pixel cross-track resolution [m]': 29.73, 'SNR': 1281.85, 'dynamic range': 205400.19, 'noise-equivalent delta T [K]': np.inf})

    def test_calc_data_metrics_MODISBand10_precomputed(self):
        """  Model instrument with MODIS Band 10 specs (Optical, whiskbroom) and test with the results as computed on 4 April 2021.

        Refer to the following article more context:
        V. Ravindra and S. Nag, "Instrument Data Metrics Evaluator for Tradespace Analysis of Earth Observing Constellations," 
        2020 IEEE Aerospace Conference, Big Sky, MT, USA, 2020.

        .. note:: The results using the LOWTRAN7 model are computed at resolution wav_step_percm = 5.

        """
        modis_band10_dict = {
                                "@type": "Passive Optical Scanner",
                                "name": "MODIS Band10",
                                "mass": 274,
                                "volume": 1.6, 
                                "power": 162.5, 
                                "fieldOfViewGeometry": {
                                    "shape": "RECTANGULAR",
                                    "angleHeight": 0.812366806011266,
                                    "angleWidth": 110
                                },
                                "scanTechnique": "WHISKBROOM",
                                "dataRate": 6.2,
                                "numberDetectorRows": 10,
                                "numberDetectorCols": 1,
                                "detectorWidth": 540e-6,
                                "focalLength": 380.859e-3,
                                "operatingWavelength": 490e-9,
                                "bandwidth": 10e-9,
                                "quantumEff": 0.33,
                                "targetBlackBodyTemp": 300,
                                "bitsPerPixel": 12,
                                "opticsSysEff": 0.8,
                                "numOfReadOutE": 25,
                                "apertureDia": 0.1778,
                                "Fnum": 2.1421,
                                "maxDetectorExposureTime": 323.333e-6,
                                "atmosLossModel": "LOWTRAN7",
                                "_comments": ["purpose is for observation of surface/ cloud temperature(note target temp)",
                                              "quantumEff, opticsSysEff, numofReadoutE are guessed."]
                            }
        modis_band10 = PassiveOpticalScannerModel.from_dict(modis_band10_dict)
        # Aqua orbit at 10 Apr 2021 15:07:56.800 UTC  (NIGHT time)                                                                          
        sc_orbit_state = {'time [JDUT1]':2459315.130520833,  'x [km]': -5054.315202286442,    'y [km]': -4878.491479401228,    'z [km]': 883.5310463297755, 
                                                            'vx [km/s]': -1.417318347731835, 'vy [km/s]': 0.1319708892386859, 'vz [km/s]': -7.367383505358474 }
        target_coords = {'lat [deg]': 7.127116160568699 , 'lon [deg]': 158.1924750010043} # nadir position of satellite
        obsv_metrics = modis_band10.calc_data_metrics(sc_orbit_state, target_coords)
        self.assertEqual(obsv_metrics['ground pixel along-track resolution [m]'], 996.11)
        self.assertEqual(obsv_metrics['ground pixel cross-track resolution [m]'], 997.19)
        self.assertEqual(obsv_metrics['SNR'], 0.0)
        self.assertEqual(obsv_metrics['dynamic range'], 0.0)
        assert math.isclose(obsv_metrics['noise-equivalent delta T [K]'], 224514118246476.62, rel_tol=1e-9)

        # 10 Apr 2021 15:55:53.269   (Day time)                                                                                                                   
        sc_orbit_state = {'time [JDUT1]':2459315.1638078704,  'x [km]': 4904.051098680667,    'y [km]': 4868.949787679997,    'z [km]': -1516.567875770611, 
                                                            'vx [km/s]': 1.903191094106026, 'vy [km/s]': 0.3436316797910688, 'vz [km/s]': 7.255195566766275 }
        target_coords = {'lat [deg]': -12.36694967995247 , 'lon [deg]': -33.02510031498068} # nadir position of satellite
        obsv_metrics = modis_band10.calc_data_metrics(sc_orbit_state, target_coords)
        self.assertEqual(obsv_metrics, {'ground pixel along-track resolution [m]': 988.98, 'ground pixel cross-track resolution [m]': 989.93, 'SNR': 3254.27, 'dynamic range': 423635.37, 'noise-equivalent delta T [K]': np.inf})
        
        
        # disable LOWTRAN7 atmospheric loss model and evaluate results
        modis_band10_dict["atmosLossModel"] = None
        modis_band10 = PassiveOpticalScannerModel.from_dict(modis_band10_dict)
        obsv_metrics = modis_band10.calc_data_metrics(sc_orbit_state, target_coords)
        self.assertEqual(obsv_metrics, {'ground pixel along-track resolution [m]': 988.98, 'ground pixel cross-track resolution [m]': 989.93, 'SNR': 3892.58, 'dynamic range': 606110.78, 'noise-equivalent delta T [K]': np.inf})

    def test_calc_data_metrics_MODISBand1_precomputed(self):
        """  Model instrument with MODIS Band 1 specs (Optical, whiskbroom) and test with the results as computed on 4 April 2021.

        Refer to the following article more context:
        V. Ravindra and S. Nag, "Instrument Data Metrics Evaluator for Tradespace Analysis of Earth Observing Constellations," 
        2020 IEEE Aerospace Conference, Big Sky, MT, USA, 2020.

        .. note:: The results using the LOWTRAN7 model are computed at resolution wav_step_percm = 5.

        """
        modis_band1_dict = {
                                "@type": "Passive Optical Scanner",
                                "name": "MODIS Band1",
                                "mass": 274,
                                "volume": 1.6, 
                                "power": 162.5, 
                                "fieldOfViewGeometry": {
                                    "shape": "RECTANGULAR",
                                    "angleHeight": 0.812366806011266,
                                    "angleWidth": 110
                                },
                                "scanTechnique": "WHISKBROOM",
                                "dataRate": 6.2,
                                "numberDetectorRows": 40,
                                "numberDetectorCols": 1,
                                "detectorWidth": 135e-6,
                                "focalLength": 380.859e-3,
                                "operatingWavelength": 645e-9,
                                "bandwidth": 50e-9,
                                "quantumEff": 0.33,
                                "targetBlackBodyTemp": 300,
                                "bitsPerPixel": 12,
                                "opticsSysEff": 0.8,
                                "numOfReadOutE": 25,
                                "apertureDia": 0.1778,
                                "Fnum": 2.1421,
                                "maxDetectorExposureTime": 73.3e-6,
                                "atmosLossModel": "LOWTRAN7",
                                "_comments": ["purpose for surface/ cloud temperature(note target temp)",
                                              "quantumEff, opticsSysEff, numofReadoutE are guessed."]
                            }
        modis_band1 = PassiveOpticalScannerModel.from_dict(modis_band1_dict)
        # Aqua orbit at 10 Apr 2021 15:07:56.800 UTC  (NIGHT time)                                                                          
        sc_orbit_state = {'time [JDUT1]':2459315.130520833,  'x [km]': -5054.315202286442,    'y [km]': -4878.491479401228,    'z [km]': 883.5310463297755, 
                                                            'vx [km/s]': -1.417318347731835, 'vy [km/s]': 0.1319708892386859, 'vz [km/s]': -7.367383505358474 }
        target_coords = {'lat [deg]': 7.127116160568699 , 'lon [deg]': 158.1924750010043} # nadir position of satellite
        obsv_metrics = modis_band1.calc_data_metrics(sc_orbit_state, target_coords)
        self.assertEqual(obsv_metrics['ground pixel along-track resolution [m]'], 249.03)
        self.assertEqual(obsv_metrics['ground pixel cross-track resolution [m]'], 249.3)
        self.assertEqual(obsv_metrics['SNR'], 0.0)
        self.assertEqual(obsv_metrics['dynamic range'], 0.0)
        assert math.isclose(obsv_metrics['noise-equivalent delta T [K]'], 10124816467.16259, rel_tol=1e-9)

        # 10 Apr 2021 15:55:53.269   (Day time)                                                                                                                   
        sc_orbit_state = {'time [JDUT1]':2459315.1638078704,  'x [km]': 4904.051098680667,    'y [km]': 4868.949787679997,    'z [km]': -1516.567875770611, 
                                                            'vx [km/s]': 1.903191094106026, 'vy [km/s]': 0.3436316797910688, 'vz [km/s]': 7.255195566766275 }
        target_coords = {'lat [deg]': -12.36694967995247 , 'lon [deg]': -33.02510031498068} # nadir position of satellite
        obsv_metrics = modis_band1.calc_data_metrics(sc_orbit_state, target_coords)
        self.assertEqual(obsv_metrics, {'ground pixel along-track resolution [m]': 247.24, 'ground pixel cross-track resolution [m]': 247.48, 'SNR': 986.24, 'dynamic range': 38932.1, 'noise-equivalent delta T [K]': np.inf})

        # disable LOWTRAN7 atmospheric loss model and evaluate results
        modis_band1_dict["atmosLossModel"] = None
        modis_band1 = PassiveOpticalScannerModel.from_dict(modis_band1_dict)
        obsv_metrics = modis_band1.calc_data_metrics(sc_orbit_state, target_coords)
        self.assertEqual(obsv_metrics, {'ground pixel along-track resolution [m]': 247.24, 'ground pixel cross-track resolution [m]': 247.48, 'SNR': 1085.11, 'dynamic range': 47123.85, 'noise-equivalent delta T [K]': np.inf})

    def test_calc_data_metrics_CCAMBlueBand_precomputed(self):
        """  Model instrument with CCAM Blue Band specs (Matrix, optical) and test with the results as computed on 4 April 2021.

        Refer to the following article more context:
        V. Ravindra and S. Nag, "Instrument Data Metrics Evaluator for Tradespace Analysis of Earth Observing Constellations," 
        2020 IEEE Aerospace Conference, Big Sky, MT, USA, 2020.

        .. note:: The results using the LOWTRAN7 model are computed at resolution wav_step_percm = 5.

        """
        ccam_blue_band_dict = {
                                "@type": "Passive Optical Scanner",
                                "name": "CCAM",
                                "fieldOfViewGeometry": {
                                    "shape": "RECTANGULAR",
                                    "angleHeight": 1.2,
                                    "angleWidth": 1.2
                                },
                                "scanTechnique": "MATRIX_IMAGER",
                                "numberDetectorRows": 2048,
                                "numberDetectorCols": 2048,
                                "detectorWidth": 5.5e-6,
                                "focalLength": 520e-3,
                                "operatingWavelength": 470e-9,
                                "bandwidth": 150e-9,
                                "quantumEff": 0.40,
                                "targetBlackBodyTemp": 290,
                                "opticsSysEff": 0.6,
                                "numOfReadOutE": 13,
                                "apertureDia": 94.6e-3,
                                "Fnum": 5.5,
                                "maxDetectorExposureTime": 678e-6,
                                "atmosLossModel": "LOWTRAN7"
                            }
        ccam_blue_band = PassiveOpticalScannerModel.from_dict(ccam_blue_band_dict)
        # Aqua orbit at 10 Apr 2021 15:07:56.800 UTC  (NIGHT time)                                                                          
        sc_orbit_state = {'time [JDUT1]':2459315.130520833,  'x [km]': -5054.315202286442,    'y [km]': -4878.491479401228,    'z [km]': 883.5310463297755, 
                                                            'vx [km/s]': -1.417318347731835, 'vy [km/s]': 0.1319708892386859, 'vz [km/s]': -7.367383505358474 }
        target_coords = {'lat [deg]': 7.127116160568699 , 'lon [deg]': 158.1924750010043} # nadir position of satellite
        obsv_metrics = ccam_blue_band.calc_data_metrics(sc_orbit_state, target_coords)
        self.assertEqual(obsv_metrics['ground pixel along-track resolution [m]'], 7.43)
        self.assertEqual(obsv_metrics['ground pixel cross-track resolution [m]'], 7.44)
        self.assertEqual(obsv_metrics['SNR'], 0.0)
        self.assertEqual(obsv_metrics['dynamic range'], 0.0)
        assert math.isclose(obsv_metrics['noise-equivalent delta T [K]'], 2302356852773662.0, rel_tol=1e-9)

        # 10 Apr 2021 15:55:53.269   (Day time)                                                                                                                   
        sc_orbit_state = {'time [JDUT1]':2459315.1638078704,  'x [km]': 4904.051098680667,    'y [km]': 4868.949787679997,    'z [km]': -1516.567875770611, 
                                                            'vx [km/s]': 1.903191094106026, 'vy [km/s]': 0.3436316797910688, 'vz [km/s]': 7.255195566766275 }
        target_coords = {'lat [deg]': -12.36694967995247 , 'lon [deg]': -33.02510031498068} # nadir position of satellite
        obsv_metrics = ccam_blue_band.calc_data_metrics(sc_orbit_state, target_coords)
        self.assertEqual(obsv_metrics, {'ground pixel along-track resolution [m]': 7.38, 'ground pixel cross-track resolution [m]': 7.38, 'SNR': 63.3, 'dynamic range': 320.71, 'noise-equivalent delta T [K]': np.inf})

        # disable LOWTRAN7 atmospheric loss model and evaluate results
        ccam_blue_band_dict["atmosLossModel"] = None
        ccam_blue_band = PassiveOpticalScannerModel.from_dict(ccam_blue_band_dict)
        obsv_metrics = ccam_blue_band.calc_data_metrics(sc_orbit_state, target_coords)
        self.assertEqual(obsv_metrics, {'ground pixel along-track resolution [m]': 7.38, 'ground pixel cross-track resolution [m]': 7.38, 'SNR': 79.01, 'dynamic range': 492.91, 'noise-equivalent delta T [K]': np.inf})
