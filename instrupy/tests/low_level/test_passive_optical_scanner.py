"""Unit tests for instrupy.passive_optical_sensor module.
"""

import unittest
import json
import numpy
import sys, os


from instrupy.passive_optical_scanner import *
from instrupy.util import Orientation, FieldOfView


class TestPassiveOpticalScanner(unittest.TestCase):



    firesat = PassiveOpticalScanner.from_json('{"@type": "Passive Optical Scanner",'
                                                '"name": "FireSat",'
                                                '"mass": 28,'
                                                '"volume": 0.12,' 
                                                '"power": 32,'
                                                '"fieldOfView": {'
                                                '   "sensorGeometry": "RECTanGULAR",'
                                                '   "alongTrackFieldOfView": 0.628,'
                                                '   "crossTrackFieldOfView": 115.8'
                                                ' },'
                                                '"scanTechnique": "WhiskBROOM",'
                                                '"orientation": {'
                                                '   "convention": "SIDE_loOK",'
                                                '   "sideLookAngle": 0'
                                                ' },'
                                                '"dataRate": 85,'
                                                '"numberOfDetectorsRowsAlongTrack": 256,'
                                                '"numberOfDetectorsColsCrossTrack": 1,'
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
                                                '"snrThreshold": 10,'
                                                '"considerAtmosLoss": true}')

    def assertNearlyZeroErrorFraction(self,a,b,fraction=0.01,msg=None):
        if abs(a-b) > abs(fraction*a):
            if msg is None:
                self.fail("The given numbers %s and %s are not near each other."%(a,b))
            else:
                self.fail(msg)

    def test_from_json_basic(self):
        
        # Test: Typical case      
        self.assertEqual(self.firesat._type, "Passive Optical Scanner")
        self.assertEqual(self.firesat.name, "FireSat")
        self.assertIsInstance(self.firesat.name, str)
        self.assertEqual(self.firesat.mass, 28)
        self.assertIsInstance(self.firesat.mass, float)
        self.assertEqual(self.firesat.volume, 0.12)
        self.assertIsInstance(self.firesat.volume, float)
        self.assertEqual(self.firesat.power, 32)
        self.assertIsInstance(self.firesat.power, float)
        self.assertEqual(self.firesat.dataRate, 85)
        self.assertIsInstance(self.firesat.dataRate, float)
        self.assertIsInstance(self.firesat.orientation, Orientation)
        self.assertEqual(self.firesat.orientation.x_rot_deg, 0)
        self.assertEqual(self.firesat.orientation.y_rot_deg, 0)
        self.assertEqual(self.firesat.orientation.z_rot_deg, 0)
        self.assertIsInstance(self.firesat.fieldOfView, FieldOfView)
        self.assertAlmostEqual(self.firesat.fieldOfView.get_rectangular_fov_specs()[0], 0.628)
        self.assertAlmostEqual(self.firesat.fieldOfView.get_rectangular_fov_specs()[1], 115.8)
        self.assertEqual(self.firesat.scanTechnique, "WHISKBROOM")
        self.assertIsInstance(self.firesat.scanTechnique, str)        
        self.assertEqual(self.firesat.detectorWidth, 30e-6)
        self.assertIsInstance(self.firesat.detectorWidth, float)
        self.assertEqual(self.firesat.focalLength, 0.7)
        self.assertIsInstance(self.firesat.focalLength, float)
        self.assertEqual(self.firesat.operatingWavelength, 4.2e-6)
        self.assertIsInstance(self.firesat.operatingWavelength, float)
        self.assertEqual(self.firesat.bandwidth, 1.9e-6)
        self.assertIsInstance(self.firesat.bandwidth, float)
        self.assertEqual(self.firesat.quantumEff, 0.5)
        self.assertIsInstance(self.firesat.quantumEff, float)
        self.assertEqual(self.firesat.targetBlackBodyTemp, 290)
        self.assertIsInstance(self.firesat.targetBlackBodyTemp, float)
        self.assertEqual(self.firesat.bitsPerPixel, 8)
        self.assertIsInstance(self.firesat.bitsPerPixel, int)
        self.assertEqual(self.firesat.opticsSysEff, 0.75)
        self.assertIsInstance(self.firesat.opticsSysEff, float)
        self.assertEqual(self.firesat.numOfReadOutE, 25)
        self.assertIsInstance(self.firesat.numOfReadOutE, float)
        self.assertEqual(self.firesat.apertureDia, 0.26)
        self.assertIsInstance(self.firesat.apertureDia, float)
        self.assertEqual(self.firesat.Fnum,  2.7)
        self.assertIsInstance(self.firesat.Fnum, float)
        self.assertEqual(self.firesat.snrThreshold, 10)
        self.assertIsInstance(self.firesat.snrThreshold, float)
        self.assertTrue(self.firesat.considerAtmosLoss)
        self.assertIsInstance(self.firesat.considerAtmosLoss, bool)
        self.assertIsInstance(self.firesat, PassiveOpticalScanner)
        self.assertIsNone(self.firesat._id)
        self.assertEqual(self.firesat._type, "Passive Optical Scanner")

        # Test of an improper orientation specification. Although physically specifing XYZ as (0,10,0) degrees is the same as specifying 
        # the side-look angle as 10 deg, the behavior of the code is to throw an error.
        with self.assertRaises(Exception):
            o = PassiveOpticalScanner.from_json('{"@type": "Passive Optical Scanner",'
                                                '"name": "FireSat",'
                                                '"mass": 28,'
                                                '"volume": 0.12,' 
                                                '"power": 32,'
                                                '"fieldOfView": {'
                                                '   "sensorGeometry": "RECTANGULAR",'
                                                '   "alongTrackFieldOfView": 0.628,'
                                                '   "crossTrackFieldOfView": 115.8'
                                                ' },'
                                                '"scanTechnique": "WHISKBROOM",'
                                                '"orientation": {'
                                                '   "convention": "XYZ",'
                                                '   "xRotation": 0,'
                                                '   "yRotation": 10,'
                                                '   "zRotation": 0'
                                                ' },'                                             
                                                '"dataRate": 85,'
                                                '"numberOfDetectorsRowsAlongTrack": 256,'
                                                '"numberOfDetectorsColsCrossTrack": 1,'
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
                                                '"snrThreshold": 10,'
                                                '"considerAtmosLoss": true}')

        
        # Test of an improper field-of-view specification. 
        with self.assertRaises(Exception):
            o = PassiveOpticalScanner.from_json('{"@type": "Passive Optical Scanner",'
                                                '"name": "FireSat",'
                                                '"mass": 28,'
                                                '"volume": 0.12,' 
                                                '"power": 32,'
                                                '"fieldOfView": {'
                                                '   "sensorGeometry": "CONICAL",'
                                                '   "fullConeAngle": 30'
                                                ' },'
                                                '"scanTechnique": "WHISKBROOM",'
                                                '"orientation": {'
                                                '   "convention": "SIDE_LOOK",'
                                                '   "sideLookAngle": 0'
                                                ' },'
                                                '"dataRate": 85,'
                                                '"numberOfDetectorsRowsAlongTrack": 256,'
                                                '"numberOfDetectorsColsCrossTrack": 1,'
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
                                                '"snrThreshold": 10,'
                                                '"considerAtmosLoss": true}')


        # Test of an improper scanning technique specification
        with self.assertRaises(Exception):
            o = PassiveOpticalScanner.from_json('{"@type": "Passive Optical Scanner",'
                                                '"fieldOfView": {'
                                                '   "sensorGeometry": "RECTANGULAR",'
                                                '   "alongTrackFieldOfView": 0.628,'
                                                '   "crossTrackFieldOfView": 115.8'
                                                ' },'
                                                '"scanTechnique": "Circular",'
                                                '"orientation": {'
                                                '   "convention": "SIDE_LOOK",'
                                                '   "sideLookAngle": 0'
                                                ' },'
                                                '"dataRate": 85,'
                                                '"numberOfDetectorsRowsAlongTrack": 256,'
                                                '"numberOfDetectorsColsCrossTrack": 1,'
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
                                                '"snrThreshold": 10,'
                                                '"considerAtmosLoss": true}')


        # Test of an PUSHBROOM scanning technique specification and more than one :code:`numberOfDetectorsRowsAlongTrack` specification.
        with self.assertRaises(Exception):
            o = PassiveOpticalScanner.from_json('{"@type": "Passive Optical Scanner",'
                                                '"name": "FireSat",'
                                                '"mass": 28,'
                                                '"volume": 0.12,' 
                                                '"power": 32,'
                                                '"fieldOfView": {'
                                                '   "sensorGeometry": "RECTANGULAR",'
                                                '   "alongTrackFieldOfView": 0.628,'
                                                '   "crossTrackFieldOfView": 115.8'
                                                ' },'
                                                '"scanTechnique": "PUSHBROOM",'
                                                '"orientation": {'
                                                '   "convention": "SIDE_LOOK",'
                                                '   "sideLookAngle": 0'
                                                ' },'
                                                '"dataRate": 85,'
                                                '"numberOfDetectorsRowsAlongTrack": 10,'
                                                '"numberOfDetectorsColsCrossTrack": 1000,'
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
                                                '"snrThreshold": 10,'
                                                '"considerAtmosLoss": true}')

        # Test of an WHISKBROOM scanning technique specification and more than one :code:`numberOfDetectorsColsCrossTrack` specification.
        with self.assertRaises(Exception):
            o = PassiveOpticalScanner.from_json('{"@type": "Passive Optical Scanner",'
                                                '"name": "FireSat",'
                                                '"mass": 28,'
                                                '"volume": 0.12,' 
                                                '"power": 32,'
                                                '"fieldOfView": {'
                                                '   "sensorGeometry": "RECTANGULAR",'
                                                '   "alongTrackFieldOfView": 0.628,'
                                                '   "crossTrackFieldOfView": 115.8'
                                                ' },'
                                                '"scanTechnique": "PUSHBROOM",'
                                                '"orientation": {'
                                                '   "convention": "SIDE_LOOK",'
                                                '   "sideLookAngle": 0'
                                                ' },'
                                                '"dataRate": 85,'
                                                '"numberOfDetectorsRowsAlongTrack": 10,'
                                                '"numberOfDetectorsColsCrossTrack": 100,'
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
                                                '"snrThreshold": 10,'
                                                '"considerAtmosLoss": true}')

    def test_planck_photon_integral(self):

        # Test trivial case with 0 wavelength
        self.assertAlmostEqual(PassiveOpticalScanner.planck_photon_integral(0,290), 0)

        """ Tests using online calculator from <https://www.opticsthewebsite.com/OpticsCalculators.aspx> as truth data.
            Note that the online calculator requires minimum wavelength to be set as 1e-9 um which is nearly 0 wavelength.

        """
        self.assertNearlyZeroErrorFraction(PassiveOpticalScanner.planck_photon_integral(12e-6,1500), 1.46801e+20 * 1e4)
        self.assertNearlyZeroErrorFraction(PassiveOpticalScanner.planck_photon_integral(2e-6,500), 3.36875e+14 * 1e4) 
        self.assertNearlyZeroErrorFraction(PassiveOpticalScanner.planck_photon_integral(500e-6,45), 4.10754e+15 * 1e4) 
        self.assertNearlyZeroErrorFraction(PassiveOpticalScanner.planck_photon_integral(1e9,45), 4.40891e+15 * 1e4)  # specifying 1e9 m as wavelength is to get approximately the radiance over entire spectrum
    
    def test_radianceWithEarthAsBlackBodyRadiator(self):


        """ Tests using online calculator from <https://www.opticsthewebsite.com/OpticsCalculators.aspx> as truth data for in-band radiance calculation.
            Note that the online calculator requires minimum wavelength to be set as 1e-9 um which is nearly 0 wavelength.

        """
        self.assertNearlyZeroErrorFraction(PassiveOpticalScanner.radianceWithEarthAsBlackBodyRadiator(16.25e-6, 12.5e-6, 290, 0, considerAtmosLoss = False),  5.08113e+17 * 1e4) # 10 um to 22.5 um at 290 K, 0 deg incidence angle
        self.assertNearlyZeroErrorFraction(PassiveOpticalScanner.radianceWithEarthAsBlackBodyRadiator(6e-6, 4e-6, 180.5, 0, considerAtmosLoss = False),  6.74299e+14 * 1e4) 
        self.assertNearlyZeroErrorFraction(PassiveOpticalScanner.radianceWithEarthAsBlackBodyRadiator(0.5e-6, 0.8e-6, 270, 0, considerAtmosLoss = False),  2.74990e-5 * 1e4) 
        self.assertNearlyZeroErrorFraction(PassiveOpticalScanner.radianceWithEarthAsBlackBodyRadiator(140e-6, 40e-6, 330, 0, considerAtmosLoss = False),  1.77269e+16 * 1e4)

        # Tests with 10 um to 22.5 um at 290 K at different incidence angles
        self.assertNearlyZeroErrorFraction(PassiveOpticalScanner.radianceWithEarthAsBlackBodyRadiator(16.25e-6, 12.5e-6, 290, 0.1, considerAtmosLoss = False),  5.08113e+17 * 1e4 * numpy.cos(0.1)) 
        self.assertNearlyZeroErrorFraction(PassiveOpticalScanner.radianceWithEarthAsBlackBodyRadiator(16.25e-6, 12.5e-6, 290, -0.1, considerAtmosLoss = False),  5.08113e+17 * 1e4 * numpy.cos(0.1))
        self.assertNearlyZeroErrorFraction(PassiveOpticalScanner.radianceWithEarthAsBlackBodyRadiator(16.25e-6, 12.5e-6, 290, 0.5 + 2*3.141, considerAtmosLoss = False),  5.08113e+17 * 1e4 * numpy.cos(0.5))
        self.assertNearlyZeroErrorFraction(PassiveOpticalScanner.radianceWithEarthAsBlackBodyRadiator(16.25e-6, 12.5e-6, 290, -0.5 - 8*3.141, considerAtmosLoss = False),  5.08113e+17 * 1e4 * numpy.cos(0.5))


        # Tests with unrealistic observation incidence angles
        with self.assertRaises(Exception):
            PassiveOpticalScanner.radianceWithEarthAsBlackBodyRadiator(16.25e-6, 12.5e-6, 290, numpy.deg2rad(91), considerAtmosLoss = False)
        with self.assertRaises(Exception):
            PassiveOpticalScanner.radianceWithEarthAsBlackBodyRadiator(16.25e-6, 12.5e-6, 290, numpy.deg2rad(-91), considerAtmosLoss = False)
        with self.assertRaises(Exception):
            PassiveOpticalScanner.radianceWithEarthAsBlackBodyRadiator(16.25e-6, 12.5e-6, 290, numpy.deg2rad(-91 + 360*8), considerAtmosLoss = False)

    
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
        self.assertAlmostEqual(PassiveOpticalScanner.radianceWithEarthAsReflector(0.5e-9, 0.2e-9, tObs_JDUT1, obs_pos_km, tar_pos_km, obs_area_m2, considerAtmosLoss = True), 0)
        self.assertAlmostEqual(PassiveOpticalScanner.radianceWithEarthAsReflector(1, 1e-2, tObs_JDUT1, obs_pos_km, tar_pos_km, obs_area_m2, considerAtmosLoss = True), 0, places = 3)

        """ Test: Reflected energy for visible wavelengths must be greater than that of other wavelengths, keeping bandwidth same.
        """        
        self.assertGreater(PassiveOpticalScanner.radianceWithEarthAsReflector(0.5e-6, 0.2e-6, tObs_JDUT1, obs_pos_km, tar_pos_km, obs_area_m2, considerAtmosLoss = True), PassiveOpticalScanner.radianceWithEarthAsReflector(6e-6, 0.2e-6, tObs_JDUT1, obs_pos_km, tar_pos_km, obs_area_m2, considerAtmosLoss = True)) # longer wavelengths than visible
        self.assertGreater(PassiveOpticalScanner.radianceWithEarthAsReflector(0.5e-6, 0.2e-6, tObs_JDUT1, obs_pos_km, tar_pos_km, obs_area_m2, considerAtmosLoss = True), PassiveOpticalScanner.radianceWithEarthAsReflector(0.25e-6, 0.2e-6, tObs_JDUT1, obs_pos_km, tar_pos_km, obs_area_m2, considerAtmosLoss = True)) # shorter wavelengths than visible

        """
            Test with another observer position, one where the observer sees the target pixel (of fixed area) at a larger angle.
        """
        opWav_m = 0.5e-6
        bw_m = 0.2e-6
        obs2_pos_km = [6893.654271085462, -9.186593534864809e-013, 1215.537243668513]
        self.assertGreater(PassiveOpticalScanner.radianceWithEarthAsReflector(opWav_m, bw_m, tObs_JDUT1, obs_pos_km, tar_pos_km, obs_area_m2, considerAtmosLoss = True), PassiveOpticalScanner.radianceWithEarthAsReflector(opWav_m, bw_m, tObs_JDUT1, obs2_pos_km, tar_pos_km, obs_area_m2, considerAtmosLoss = True)) # longer wavelengths than visible


    
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
        # The InstruPy computed value must be greater than truth value since SMAD does not cosdier the energy reflected off Sun, and the date corresponds to satellite over target at day-time.
        self.assertGreater(PassiveOpticalScanner.calculate_number_of_signal_electrons(opWav_m, bw_m, bbT_K, apDia_m, opTrns, QE, tObs_JDUT1, obs_pos_km, tar_pos_km, pixelArea_m2, Ti_s, considerAtmosLoss = True), 8286.104444633884)        

    def test_calc_typ_data_metrics(self):
        
        # Test: Truth data from SMAD 3rd edition, Table 9-15. Note that the  inputs are made so that they are consistent with the example.
        # Further note that SMAD does not consider the energy from Sun reflected off Sun in their calculations.
        epoch_JDUT1 = 2451623.999630
        SpacecraftOrbitState = {'Time[JDUT1]':epoch_JDUT1, 'x[km]': 7078.137, 'y[km]': 0, 'z[km]': 0, 'vx[km/s]': 0, 'vy[km/s]': 7.5, 'vz[km/s]': 0} # equatorial orbit, altitude about 700 km
        TargetCoords = {'Lat [deg]': 0, 'Lon [deg]': 0} # lat = 0, lon = 0 corresponds to [6378, 0, 0] km in ECI for observer position, check using Matlab function: eci2lla([6378, 0, 0] ,[2000 3 20 11 59 28.000])
        obsv_metrics = self.firesat.calc_typ_data_metrics(SpacecraftOrbitState, TargetCoords)
        self.assertAlmostEqual(obsv_metrics["Ground Pixel Along-Track Resolution [m]"], 30, delta = 3)
        self.assertAlmostEqual(obsv_metrics["Ground Pixel Cross-Track Resolution [m]"], 30, delta = 3)
        # A (positive) deviation is expected since SMAD does not consider the energy from Sun reflected off Sun
        self.assertGreater(obsv_metrics["SNR"], 88)
        self.assertGreater(obsv_metrics["DR"], 332.9)
        self.assertAlmostEqual(obsv_metrics["Noise-Equivalent Delta T [K]"], 0.3 , places = 2)
        self.assertTrue(obsv_metrics["Coverage [T/F]"])

        ## To do, make test with satellite in night region, thus no reflected enregy of Sun. Result should match SMAD.

    
    def test_calculate_integration_time(self):

        # Test: PUSHBROOM scanning
        self.assertAlmostEqual(PassiveOpticalScanner.calculate_integration_time("PUSHBROOM",1, 1, 12.5, 1), 12.5)
        self.assertAlmostEqual(PassiveOpticalScanner.calculate_integration_time("PushbroOM",1, 150, 12.5, 0.1), 12.5)
        self.assertAlmostEqual(PassiveOpticalScanner.calculate_integration_time("PUSHBROOM",1, 150, 12.5, 0.1, crossTrack_fov_deg = 30), 12.5)
        # check max exposure time functionality
        self.assertAlmostEqual(PassiveOpticalScanner.calculate_integration_time("PUSHBROOM",1, 150, 12.5, 0.1, maxDetectorExposureTime = 5), 5)
        
        with self.assertRaises(Exception): # Exception expected if number of detector rows in along-track direction is not 1.
            PassiveOpticalScanner.calculate_integration_time("PUSHBROOM",10, 150, 12.5, 1)

        # Test: Whiskroom scanning
        self.assertAlmostEqual(PassiveOpticalScanner.calculate_integration_time("WHISKBROOM",1, 1, 12.5, 1, maxDetectorExposureTime = None, crossTrack_fov_deg = 30), 12.5/30)
        self.assertAlmostEqual(PassiveOpticalScanner.calculate_integration_time("WhiskbrOOM",1, 1, 12.5, 0.1, maxDetectorExposureTime = None, crossTrack_fov_deg = 30), 12.5/300)
        self.assertAlmostEqual(PassiveOpticalScanner.calculate_integration_time("WHISKBRooM",20, 1, 12.5, 0.1, maxDetectorExposureTime = None, crossTrack_fov_deg = 30), 12.5/300)
        # check max exposure time functionality
        self.assertAlmostEqual(PassiveOpticalScanner.calculate_integration_time("WHISKBRooM",20, 1, 12.5, 0.1, maxDetectorExposureTime = (12.5/3000), crossTrack_fov_deg = 30), 12.5/3000)

        with self.assertRaises(Exception): # Exception expected if number of detector columns in cross-track direction is not 1.
            PassiveOpticalScanner.calculate_integration_time("WHISKBROOM",10, 150, 12.5, 1, maxDetectorExposureTime = None, crossTrack_fov_deg = 30)

        with self.assertRaises(Exception): # Exception expected if cross-track-fov is not specified
            PassiveOpticalScanner.calculate_integration_time("WHISKBROOM",20, 1, 12.5, 0.1, maxDetectorExposureTime = None, crossTrack_fov_deg = None)

        # Test: MATRIX_IMAGER scanning
        self.assertAlmostEqual(PassiveOpticalScanner.calculate_integration_time("MATRIX_IMAGER",1, 1, 12.5, 1), 12.5)
        self.assertAlmostEqual(PassiveOpticalScanner.calculate_integration_time("MATRIX_IMAGER",1, 1, 12.5, 0.1), 12.5)
        self.assertAlmostEqual(PassiveOpticalScanner.calculate_integration_time("Matrix_Imager",20, 1, 12.5, 0.1), 12.5)
        self.assertAlmostEqual(PassiveOpticalScanner.calculate_integration_time("matrix_imager",1, 20, 12.5, 0.1), 12.5)
        self.assertAlmostEqual(PassiveOpticalScanner.calculate_integration_time("mATrIX_ImAGEr",20, 20, 12.5, 0.1, 30), 12.5)
        # check max exposure time functionality
        self.assertAlmostEqual(PassiveOpticalScanner.calculate_integration_time("matrix_imager",1, 20, 12.5, 0.1, maxDetectorExposureTime = 1.2), 1.2)
        
        