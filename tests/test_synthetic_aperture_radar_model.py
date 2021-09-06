"""Unit tests for instrupy.synthetic_aperture_radar_model.

Tests:

* test_from_json_basic_1: Test typical initialization of the synthetic aperture radar

* test_prf_constraint_eval: Validate (single/compact pol) PRF selection with external data. Truth data from the "Spaceborne SAR Study: LDRD 92 Final Report" SANDIA Report.

* test_calc_data_metrics_seasat: Validate with Seasat data in the "Spaceborne SAR Study: LDRD 92 Final Report" SANDIA Report. Large difference observed in case of azimuthal resolution.

* test_calc_data_metrics_microxsar: Validate with Mixroxsar data in IEICE transaction article.

* test_from_json_calc_data_metrics_test_sar1: Test with precomputed datametrics. Aspects: 

    * Altitude = 500km, 
    * Different incidence angles: 60 deg, 35 deg, 45 deg, 55 deg.
    * L-band.
    * Flag ``instru_look_angle_from_target_inc_angle`` is set to ``True``.
    * Full-swath.
    * Compact-pol (same result as Single Pol).

* test_from_json_calc_data_metrics_test_sar2: Test with precomputed datametrics. Aspects:

    * ALtitude = 500km.
    * Different incidence angles: 60 deg, 35 deg, 45 deg, 55 deg.
    * L-band.
    * Flag ``instru_look_angle_from_target_inc_angle`` is set to ``True``.
    * Full-swath.
    * Dual-pol, AIRSAR.

* test_from_json_calc_data_metrics_test_sar3: Test with precomputed datametrics. Aspects:

    * ALtitude = 500km.
    * Different incidence angles: 60 deg, 35 deg, 45 deg, 55 deg.
    * L-band.
    * Flag ``instru_look_angle_from_target_inc_angle`` is set to ``True``.
    * Full-swath.
    * Dual-pol, SMAP.

* test_from_json_calc_data_metrics_test_sar4: Test with precomputed datametrics. Aspects:

    * ALtitude = 500km.
    * Different incidence angles: 60 deg, 35 deg, 45 deg, 55 deg.
    * P-band.
    * Flag ``instru_look_angle_from_target_inc_angle`` is set to ``True``.
    * Fixed-swath 50km.
    * Dual-pol, SMAP.

* test_from_json_calc_data_metrics_test_sar5: Test with precomputed datametrics. Aspects:

    * ALtitude = 600km.
    * Different incidence angles: 60 deg, 35 deg, 45 deg, 55 deg.
    * L-band.
    * Flag ``instru_look_angle_from_target_inc_angle`` is set to ``True``.
    * Fixed-swath 25km.
    * Dual-pol, AIRSAR.

* test_from_dict_calc_data_metrics_scansar: Test ScanSAR operation. Aspects:
            
    * ScanSAR, numSubSwaths = 1, 2, 5. When  numSubSwaths=1, the result is the same as from Stripmap.
    * ALtitude = 600km.
    * Different incidence angles: 15 deg and 30 deg.
    * Flag ``instru_look_angle_from_target_inc_angle`` is set to ``True``.
    * Fixed-swath 25km.            
    * Dual-pol, SMAP, custom pulse-seperation
    * Custom atmospheric loss

"""

import unittest
import json
import numpy as np
import sys, os

from instrupy.synthetic_aperture_radar_model import SyntheticApertureRadarModel, ScanTech, PolTypeSAR, DualPolPulseConfig, SwathTypeSAR
from instrupy.util import Orientation, SphericalGeometry, ViewGeometry, FileUtilityFunctions, Maneuver

RE = 6378.137
def orbital_speed(alt_km):
    return np.sqrt(398600.5/(RE + alt_km))

############# Some SAR JSON definitions to be used in the tests. #############
# MICROXSAR
# H. Saito et al., “Compact x-band synthetic aperture radar for 100kg class satellite,” IEICE Transactions on
# Communications, vol. E100.B, no. 9, pp. 1653–1660, 2017
microxsar_json_str =   '{"@type": "Synthetic Aperture Radar",' \
                        '"name": "MiroXSAR",'  \
                        '"mass": 130,' \
                        '"volume": 0.343,' \
                        '"power": 1100,' \
                        '"orientation": {' \
                        '   "referenceFrame": "SC_BODY_FIXED",' \
                        '    "convention": "SIDE_LOOK",' \
                        '    "sideLookAngle": 30' \
                        '},' \
                        '"dataRate": 2000,' \
                        '"bitsPerPixel": 16,' \
                        '"pulseWidth": 31e-6,' \
                        '"antenna":{"shape": "RECTANGULAR", "height": 4.9, "width": 0.7, "apertureEfficiency": 0.5, "apertureExcitationProfile": "UNIFORM"},' \
                        '"operatingFrequency": 9.65e9,' \
                        '"peakTransmitPower": 1e3,' \
                        '"chirpBandwidth": 75e6,' \
                        '"minimumPRF": 3000,' \
                        '"maximumPRF": 8000,' \
                        '"radarLoss": 3.5,' \
                        '"sceneNoiseTemp": 290,' \
                        '"systemNoiseFigure": 4.3' \
                        '}'

# SEASAT
seasat_json_str =  '{"@type": "Synthetic Aperture Radar",' \
                    '"orientation": {' \
                    '   "referenceFrame": "SC_BODY_FIXED",' \
                    '    "convention": "SIDE_LOOK",' \
                    '    "sideLookAngle": 20.5' \
                    '},' \
                    '"pulseWidth": 33.4e-6,' \
                    '"antenna":{"shape": "RECTANGULAR", "height": 10.7, "width": 2.16, "apertureEfficiency": 0.6, "apertureExcitationProfile": "UNIFORM"},' \
                    '"operatingFrequency": 1.2757e9,' \
                    '"peakTransmitPower": 1000,' \
                    '"chirpBandwidth": 19e6,' \
                    '"minimumPRF": 1463,' \
                    '"maximumPRF": 1686,' \
                    '"radarLoss": 3.5,' \
                    '"systemNoiseFigure": 5.11' \
                    '}'

# ERS1
ers1_json_str = '{"@type": "Synthetic Aperture Radar",' \
                '   "orientation": {' \
                '   "referenceFrame": "SC_BODY_FIXED",' \
                '      "convention": "SIDE_LOOK",' \
                '      "sideLookAngle": 20' \
                '  },' \
                '  "pulseWidth": 37.1e-6,' \
                '  "antenna":{"shape": "RECTANGULAR", "height": 10, "width": 1, "apertureEfficiency": 0.26, "apertureExcitationProfile": "UNIFORM"},' \
                '  "operatingFrequency": 5.25e9,' \
                '  "peakTransmitPower": 4800,' \
                '  "chirpBandwidth": 15.5e6,' \
                '  "minimumPRF": 1680,' \
                '  "maximumPRF": 1700,' \
                '  "radarLoss": 3.5,' \
                '  "systemNoiseFigure": 3.4' \
                '}'

class TestSyntheticApertureRadarModel(unittest.TestCase):   

    def __init__(self, *args, **kwargs):
        # 
        super(TestSyntheticApertureRadarModel, self).__init__(*args, **kwargs)

    def test_from_json_basic_1(self):
        """ Test typical initialization of the synthetic aperture radar. 
        """   
        microxsar= SyntheticApertureRadarModel.from_json(microxsar_json_str) 
        # SAR without specification of maneuverability, numStripsInScene, altitude fields
        self.assertIsNotNone(microxsar._id) # random id is assigned when "@id" field is missing in the specs.
        self.assertEqual(microxsar._type, "Synthetic Aperture Radar")
        self.assertIsInstance(microxsar.name, str)
        self.assertEqual(microxsar.name, "MiroXSAR")
        self.assertIsInstance(microxsar.mass, float)
        self.assertEqual(microxsar.mass, 130)
        self.assertIsInstance(microxsar.volume, float)
        self.assertEqual(microxsar.volume, 0.343)
        self.assertIsInstance(microxsar.power, float)
        self.assertEqual(microxsar.power, 1100)
        self.assertEqual(microxsar.orientation, Orientation.from_json({"referenceFrame": "SC_BODY_FIXED", "convention":"Euler", "eulerSeq1":1, "eulerSeq2":2, "eulerSeq3":3, "eulerAngle1":0, "eulerAngle2":30, "eulerAngle3":0}))
        self.assertIsInstance(microxsar.orientation, Orientation)
        self.assertIsInstance(microxsar.fieldOfView, ViewGeometry)         
        self.assertAlmostEqual(microxsar.fieldOfView.sph_geom.get_fov_height_and_width()[0], 0.3197, places = 2)
        self.assertAlmostEqual(microxsar.fieldOfView.sph_geom.get_fov_height_and_width()[1], 2.2377, places = 2)
        self.assertEqual(microxsar.sceneFieldOfView, microxsar.fieldOfView)
        self.assertIsNone(microxsar.maneuver)
        self.assertIsNone(microxsar.fieldOfRegard)
        self.assertIsNone(microxsar.pointingOption)
        self.assertIsInstance(microxsar.dataRate, float)
        self.assertEqual(microxsar.dataRate, 2000)
        self.assertIsInstance(microxsar.bitsPerPixel, int)
        self.assertEqual(microxsar.bitsPerPixel, 16)        
        self.assertIsInstance(microxsar.pulseWidth, float)
        self.assertEqual(microxsar.pulseWidth, 31e-6)
        self.assertIsInstance(microxsar.antenna.height, float)
        self.assertEqual(microxsar.antenna.height, 4.9)
        self.assertIsInstance(microxsar.antenna.width, float)
        self.assertEqual(microxsar.antenna.width, 0.7)
        self.assertIsInstance(microxsar.antenna.apertureEfficiency, float)
        self.assertEqual(microxsar.antenna.apertureEfficiency, 0.5)
        self.assertEqual(microxsar.antenna.apertureExcitationProfile.value, "UNIFORM")
        self.assertIsInstance(microxsar.operatingFrequency, float)
        self.assertEqual(microxsar.operatingFrequency, 9.65e9)
        self.assertIsInstance(microxsar.peakTransmitPower, float)
        self.assertEqual(microxsar.peakTransmitPower, 1e3)
        self.assertIsInstance(microxsar.chirpBandwidth, float)
        self.assertEqual(microxsar.chirpBandwidth, 75e6)
        self.assertIsInstance(microxsar.minimumPRF, float)
        self.assertEqual(microxsar.minimumPRF, 3000)
        self.assertIsInstance(microxsar.maximumPRF, float)
        self.assertEqual(microxsar.maximumPRF, 8000)
        self.assertIsInstance(microxsar.sceneNoiseTemp, float)
        self.assertEqual(microxsar.sceneNoiseTemp, 290)
        self.assertIsInstance(microxsar.radarLoss, float)
        self.assertEqual(microxsar.radarLoss, 3.5)
        self.assertIsInstance(microxsar.atmosLoss, float)
        self.assertEqual(microxsar.atmosLoss, 2) # default value        
        self.assertIsInstance(microxsar.systemNoiseFigure, float)
        self.assertEqual(microxsar.systemNoiseFigure, 4.3)
        self.assertIsInstance(microxsar.polType, PolTypeSAR) 
        self.assertEqual(microxsar.polType, PolTypeSAR.SINGLE)
        self.assertIsNone(microxsar.dualPolPulseConfig)
        self.assertIsNone(microxsar.dualPolPulseSep)        
        self.assertIsInstance(microxsar.scanTechnique, ScanTech) 
        self.assertEqual(microxsar.scanTechnique, ScanTech.STRIPMAP) 
        self.assertIsInstance(microxsar.swathType, SwathTypeSAR) 
        self.assertEqual(microxsar.swathType, SwathTypeSAR.FULL) # default value
        self.assertIsNone(microxsar.fixedSwathSize) 
        self.assertIsInstance(microxsar.numSubSwaths, int) 
        self.assertEqual(microxsar.numSubSwaths, 1) # default value 

        # Test of an improper orientation specification. Although physically specifing XYZ as (0,10,0) degrees is the same as specifying 
        # the side-look angle as 10 deg, the behavior of the code is to throw an error.
        with self.assertRaises(Exception):
            o = SyntheticApertureRadarModel.from_json( '{"@type": "Synthetic Aperture Radar",'
                                                  '"name": "MiroXSAR",'  
                                                  '"mass": 130,' 
                                                  '"volume": 0.343,' 
                                                  '"power": 1100,' 
                                                  '"orientation": {'
                                                  '   "referenceFrame": "SC_BODY_FIXED",'
                                                  '   "convention": "Euler",'
                                                  '   "eulerAngle1": 0,'
                                                  '   "eulerAngle2": 10,'
                                                  '   "eulerAngle3": 0'
                                                  ' },'
                                                  '"dataRate": 2000,'
                                                  '"bitsPerPixel": 16,'
                                                  '"pulseWidth": 31e-6,'
                                                  '"antenna":{"shape": "RECTANGULAR", "height": 4.9, "width": 0.7, "apertureEfficiency": 0.5, "apertureExcitationProfile": "UNIFORM"},'
                                                  '"operatingFrequency": 9.65e9,' 
                                                  '"peakTransmitPower": 1e3,' 
                                                  '"chirpBandwidth": 75e6,'      
                                                  '"minimumPRF": 3000,' 
                                                  '"maximumPRF": 8000,' 
                                                  '"radarLosses": 3.5,' 
                                                  '"sceneNoiseTemp": 290,' 
                                                  '"systemNoiseFigure": 4.3,'
                                                  '}')

        # Check for invalid PRF min, max specification
        with self.assertRaises(Exception):
            o = SyntheticApertureRadarModel.from_json( '{"@type": "Synthetic Aperture Radar",'
                                                  '"name": "MiroXSAR",'  
                                                  '"mass": 130,' 
                                                  '"volume": 0.343,' 
                                                  '"power": 1100,' 
                                                  '"orientation": {'
                                                  '   "referenceFrame": "SC_BODY_FIXED",'
                                                  '   "convention": "SIDE_LOOK",'
                                                  '   "sideLookAngle": 30'
                                                  '},'
                                                  '"dataRate": 2000,'
                                                  '"bitsPerPixel": 16,'
                                                  '"pulseWidth": 31e-6,'
                                                  '"antenna":{"shape": "RECTANGULAR", "height": 4.9, "width": 0.7, "apertureEfficiency": 0.5},'
                                                  '"operatingFrequency": 9.65e9,' 
                                                  '"peakTransmitPower": 1e3,' 
                                                  '"chirpBandwidth": 75e6,'      
                                                  '"minimumPRF": 8000,' 
                                                  '"maximumPRF": 3000,' 
                                                  '"radarLosses": 3.5,' 
                                                  '"sceneNoiseTemp": 290,' 
                                                  '"systemNoiseFigure": 4.3,'
                                                  '}')

    def test_prf_constraint_eval(self):
        """ Validate (single/compact pol) PRF selection with external data.
            Truth data from Spaceborne SAR Study: LDRD 92 Final Report SANDIA Report March 1993, Page 35
            Using look angle = 18.5 deg corresponding to incidence angle of 20 deg.
            Choose Delev (antennaWidth) and fc (center frequency) to get a swath of 10 km.      

        """
        # PRF min is 2538 Hz
        f_P_min = 1
        f_P_max = 2538
        v_sc = 7.613
        v_x = 7.0596
        alt_km = 500
        instru_look_angle_rad = np.deg2rad(18.5) 
        tau_p = 30e-6
        D_az = 6
        D_elv = 1.75
        fc = 9.65e9
        self.assertIsNone(SyntheticApertureRadarModel.prf_constraint_eval(f_P_min, f_P_max, v_sc, v_x, alt_km, instru_look_angle_rad, tau_p, D_az, D_elv, fc, pol_type=PolTypeSAR.SINGLE, swath_type=SwathTypeSAR.FULL)[0])

        # PRF max is 9331 Hz
        f_P_min = 9331
        f_P_max = 20000
        v_sc = 7.613
        v_x = 7.0596
        alt_km = 500
        instru_look_angle_rad = np.deg2rad(18.5) 
        tau_p = 30e-6
        D_az = 6
        D_elv = 1.75
        fc = 9.65e9
        self.assertIsNone(SyntheticApertureRadarModel.prf_constraint_eval(f_P_min, f_P_max, v_sc, v_x, alt_km, instru_look_angle_rad, tau_p, D_az, D_elv, fc, pol_type=PolTypeSAR.COMPACT, swath_type=SwathTypeSAR.FULL)[0])

        f_P_min = 2500
        f_P_max = 3000
        v_sc = 7.613
        v_x = 7.0596
        alt_km = 500
        instru_look_angle_rad = np.deg2rad(18.5)
        tau_p = 30e-6
        D_az = 6
        D_elv = 1.75 
        fc = 9.65e9
        self.assertEqual(SyntheticApertureRadarModel.prf_constraint_eval(f_P_min, f_P_max, v_sc, v_x, alt_km, instru_look_angle_rad, tau_p, D_az, D_elv, fc, pol_type=PolTypeSAR.SINGLE)[0], 3000)

        f_P_min = 2500
        f_P_max = 4000
        v_sc = 7.613
        v_x = 7.0596        
        alt_km = 500
        instru_look_angle_rad = np.deg2rad(18.5)
        tau_p = 30e-6
        D_az = 6
        D_elv = 1.75 
        fc = 9.65e9
        self.assertEqual(SyntheticApertureRadarModel.prf_constraint_eval(f_P_min, f_P_max, v_sc, v_x, alt_km, instru_look_angle_rad, tau_p, D_az, D_elv, fc, swath_type=SwathTypeSAR.FULL)[0], 3916)

        f_P_min = 2500
        f_P_max = 4500
        v_sc = 7.613
        v_x = 7.0596        
        alt_km = 500
        instru_look_angle_rad = np.deg2rad(18.5)
        tau_p = 30e-6
        D_az = 6
        D_elv = 1.75 
        fc = 9.65e9
        self.assertEqual(SyntheticApertureRadarModel.prf_constraint_eval(f_P_min, f_P_max, v_sc, v_x, alt_km, instru_look_angle_rad, tau_p, D_az, D_elv, fc)[0], 4184)

    def test_calc_data_metrics_seasat(self):
        """ SeaSatA typical data-metrics. Truth values taken from the reference (below) are approximately equal to those computed by this
            test. There is ambiguity in the actual operating point used in computation of the truth value (eg: exact PRF within
            the given range, pulse widths.
           
            Ref: D. Bickel, B. Brock, and C. Allen, Spaceborne SAR study: LDRD 92 final report, Mar 1993.
        """
        # 
        seasat = SyntheticApertureRadarModel.from_json(seasat_json_str)
        epoch_JDUT1 = 2451623.999630 # 2000 3 20 11 59 28.000
        # lat = 0, lon = 0 corresponds to [6378, 0, 0] km in ECI for observer position, check using Matlab function: eci2lla([6378, 0, 0] ,[2000 3 20 11 59 28.000])
        sc_orbit_state = {'time [JDUT1]':epoch_JDUT1, 'x [km]': 6378.137 + 775, 'y [km]': 0, 'z [km]': 0, 'vx [km/s]': 0, 'vy [km/s]': orbital_speed(775), 'vz [km/s]': 0} # equatorial orbit, altitude 600 km
        target_coords = {'lat [deg]': 1.89, 'lon [deg]': 0} # incidence angle to the target is 23 deg
        obsv_metrics = seasat.calc_data_metrics(sc_orbit_state, target_coords, instru_look_angle_from_target_inc_angle=False)
        self.assertAlmostEqual(obsv_metrics["ground pixel along-track resolution [m]"]*4, 25, delta=8) # the 4-look resolution is specified in the LDRD92 report
        self.assertAlmostEqual(obsv_metrics["ground pixel cross-track resolution [m]"], 25, delta=1)
        self.assertAlmostEqual(obsv_metrics["NESZ [dB]"], -24, delta=2)
        self.assertAlmostEqual(obsv_metrics["incidence angle [deg]"], 23, delta=0.1)
        self.assertAlmostEqual(obsv_metrics["swath-width [km]"], 100, delta=5)
        self.assertAlmostEqual(obsv_metrics["PRF [Hz]"], 1647, delta=50) # upper limit in the reference is 1647 Hz

    def test_calc_data_metrics_microxsar(self):
        """ MicroXSAR typical data-metrics. Truth values taken from the reference (below) are approximately equal to those computed by this
            test. There is ambiguity in the actual operating point used in computation of the truth value (eg: exact PRF within
            the given range, pulse width.)s.
           
            Ref: H. Saito et al., “Compact x-band synthetic aperture radar for 100kg class satellite,” IEICE Transactions on
            Communications, vol. E100.B, no. 9, pp. 1653–1660, 2017

        """
        microxsar = SyntheticApertureRadarModel.from_json(microxsar_json_str)
        epoch_JDUT1 = 2451623.999630# 2000 3 20 11 59 28.000
        # lat = 0, lon = 0 corresponds to [6378, 0, 0] km in ECI for observer position, check using Matlab function: eci2lla([6378, 0, 0] ,[2000 3 20 11 59 28.000])
        sc_orbit_state = {'time [JDUT1]':epoch_JDUT1, 'x [km]': 6378.137 + 600, 'y [km]': 0, 'z [km]': 0, 'vx [km/s]': 0, 'vy [km/s]': orbital_speed(600), 'vz [km/s]': 0} # equatorial orbit, altitude 600 km
        target_coords = {'lat [deg]': 2, 'lon [deg]': 0} # incidence angle to the target is 28.9855 deg
        obsv_metrics = microxsar.calc_data_metrics(sc_orbit_state, target_coords)
        self.assertAlmostEqual(obsv_metrics["ground pixel along-track resolution [m]"], 2.09, delta = 0.1)
        self.assertAlmostEqual(obsv_metrics["ground pixel cross-track resolution [m]"], 4.95, delta = 0.1)
        self.assertAlmostEqual(obsv_metrics["NESZ [dB]"], -13, delta = 0.5)
        self.assertAlmostEqual(obsv_metrics["incidence angle [deg]"], 28.9855, delta = 0.1)
        self.assertAlmostEqual(obsv_metrics["swath-width [km]"], 37, delta = 1)
        self.assertAlmostEqual(obsv_metrics["PRF [Hz]"], 3751)

        # test with the instru_look_angle_from_target_inc_angle flag set to True
        # the instrument look angle is considered from the incidence angle at the target ground-point (28.99 deg).
        # Note that the calculated swath-width is smaller. A different PRF is evaluated (higher PRF) which results in a different (improved) NESZ. 
        obsv_metrics = microxsar.calc_data_metrics(sc_orbit_state, target_coords, instru_look_angle_from_target_inc_angle=True)
        self.assertAlmostEqual(obsv_metrics["ground pixel along-track resolution [m]"], 2.09, delta = 0.1)
        self.assertAlmostEqual(obsv_metrics["ground pixel cross-track resolution [m]"], 4.95, delta = 0.1)
        self.assertAlmostEqual(obsv_metrics["NESZ [dB]"], -14.12, delta = 0.5)
        self.assertAlmostEqual(obsv_metrics["incidence angle [deg]"], 28.9855, delta = 0.1)
        self.assertAlmostEqual(obsv_metrics["swath-width [km]"], 34.4, delta = 1)
        self.assertAlmostEqual(obsv_metrics["PRF [Hz]"], 4976)

    def test_from_json_calc_data_metrics_test_sar1(self):
        """ Test of SAR with specs as given by ``test_sar1_json_str``. Following aspects of the SAR are in focus:
            
            * ALtitude = 500km.
            * Different incidence angles: 60 deg, 35 deg, 45 deg, 55 deg.
            * L-band.
            * Flag ``instru_look_angle_from_target_inc_angle`` is set to ``True``.
            * Full-swath.
            * Single/compact -pol (same result).

            Tests:

            * Initialization (partly).
            * Resulting datametrics are compared with precomputed values (generated on 2 April 2021). 

        """
        
        test_sar1_json_str = '{"@type": "Synthetic Aperture Radar",' \
                    '"orientation": {' \
                    '   "referenceFrame": "SC_BODY_FIXED",' \
                    '    "convention": "SIDE_LOOK",' \
                    '    "sideLookAngle": 30' \
                    '},' \
                    '"pulseWidth": 6.61e-6,' \
                    '"antenna":{"shape": "RECTANGULAR", "height": 7.01, "width": 6.58, "apertureEfficiency": 0.6},' \
                    '"operatingFrequency": 1275.7e6,' \
                    '"peakTransmitPower": 1000,' \
                    '"chirpBandwidth": 3.23e6,' \
                    '"minimumPRF": 1,' \
                    '"maximumPRF": 20000,' \
                    '"radarLoss": 2,' \
                    '"systemNoiseFigure": 2,' \
                    '"swathConfig": {' \
                    '   "@type": "full"' \
                    '},' \
                    '"polarization": {' \
                    '   "@type": "compact"' \
                    '}' \
                    '}'
        
        test_sar1 = SyntheticApertureRadarModel.from_json(test_sar1_json_str)
        
        self.assertIsInstance(test_sar1.swathType, SwathTypeSAR)
        self.assertEqual(test_sar1.swathType, SwathTypeSAR.FULL)
        self.assertIsInstance(test_sar1.polType, PolTypeSAR)
        self.assertEqual(test_sar1.polType, PolTypeSAR.COMPACT)
        
        h = 500 # [km]
        orb_speed = orbital_speed(h)

        inc_deg = 60
        obsv_metrics = test_sar1.calc_data_metrics(alt_km=h, sc_speed_kmps=orb_speed, sc_gnd_speed_kmps=orb_speed*(RE/(RE+h)), inc_angle_deg=inc_deg, 
                                                    instru_look_angle_from_target_inc_angle=True)
        self.assertEqual(obsv_metrics, {'ground pixel along-track resolution [m]': 3.25, 'ground pixel cross-track resolution [m]': 64.3, 'NESZ [dB]': -30.64, 'incidence angle [deg]': 60.0, 'swath-width [km]': 65.1, 'PRF [Hz]': 2385})

        inc_deg = 35
        obsv_metrics = test_sar1.calc_data_metrics(alt_km=h, sc_speed_kmps=orb_speed, sc_gnd_speed_kmps=orb_speed*(RE/(RE+h)), inc_angle_deg=inc_deg, 
                                                    instru_look_angle_from_target_inc_angle=True)
        self.assertEqual(obsv_metrics, {'ground pixel along-track resolution [m]': 3.25, 'ground pixel cross-track resolution [m]': 97.09, 'NESZ [dB]': -42.46, 'incidence angle [deg]': 35.0, 'swath-width [km]': 26.2, 'PRF [Hz]': 6896})

        inc_deg = 45
        obsv_metrics = test_sar1.calc_data_metrics(alt_km=h, sc_speed_kmps=orb_speed, sc_gnd_speed_kmps=orb_speed*(RE/(RE+h)), inc_angle_deg=inc_deg, 
                                                    instru_look_angle_from_target_inc_angle=True)
        self.assertEqual(obsv_metrics,{'ground pixel along-track resolution [m]': 3.25, 'ground pixel cross-track resolution [m]': 78.76, 'NESZ [dB]': -38.03, 'incidence angle [deg]': 45.0, 'swath-width [km]': 34.5, 'PRF [Hz]': 4518})

        inc_deg = 55
        obsv_metrics = test_sar1.calc_data_metrics(alt_km=h, sc_speed_kmps=orb_speed, sc_gnd_speed_kmps=orb_speed*(RE/(RE+h)), inc_angle_deg=inc_deg, 
                                                    instru_look_angle_from_target_inc_angle=True)
        self.assertEqual(obsv_metrics,{'ground pixel along-track resolution [m]': 3.25, 'ground pixel cross-track resolution [m]': 67.98, 'NESZ [dB]': -32.79, 'incidence angle [deg]': 55.0, 'swath-width [km]': 50.8, 'PRF [Hz]': 2663})

        
    def test_from_json_calc_data_metrics_test_sar2(self):
        """ Test of SAR with specs as given by ``test_sar2_json_str``. Following aspects of the SAR are in focus:
            
            * ALtitude = 500km.
            * Different incidence angles: 60 deg, 35 deg, 45 deg, 55 deg.
            * L-band.
            * Flag ``instru_look_angle_from_target_inc_angle`` is set to ``True``.
            * Full-swath.
            * Dual-pol, AIRSAR.

            Tests:

            * Initialization (partly).
            * Resulting datametrics are compared with precomputed values (generated on 2 April 2021). 

        """
        
        test_sar2_json_str = '{"@type": "Synthetic Aperture Radar",' \
                    '"pulseWidth": 1e-6,' \
                    '"antenna":{"shape": "RECTANGULAR", "height": 13.80, "width": 6.30, "apertureEfficiency": 0.6, "apertureExcitationProfile": "UNIFORM"},' \
                    '"operatingFrequency": 1275.7e6,' \
                    '"peakTransmitPower": 1000,' \
                    '"chirpBandwidth": 1.80e6,' \
                    '"minimumPRF": 1,' \
                    '"maximumPRF": 20000,' \
                    '"radarLoss": 2,' \
                    '"systemNoiseFigure": 2,' \
                    '"swathConfig": {' \
                    '   "@type": "full"' \
                    '},' \
                    '"polarization": {' \
                    '   "@type": "dual",' \
                    '   "pulseConfig": {' \
                    '        "@type": "AIRSAR"' \
                    '}' \
                    '}' \
                    '}'
        
        test_sar2 = SyntheticApertureRadarModel.from_json(test_sar2_json_str)
        
        self.assertIsInstance(test_sar2.swathType, SwathTypeSAR)
        self.assertEqual(test_sar2.swathType, SwathTypeSAR.FULL)
        self.assertIsInstance(test_sar2.polType, PolTypeSAR)
        self.assertEqual(test_sar2.polType, PolTypeSAR.DUAL)
        self.assertIsInstance(test_sar2.dualPolPulseConfig, DualPolPulseConfig)
        self.assertEqual(test_sar2.dualPolPulseConfig, DualPolPulseConfig.AIRSAR)
        self.assertIsInstance(test_sar2.orientation, Orientation)
        self.assertEqual(test_sar2.orientation, Orientation.from_dict({"referenceFrame": "SC_BODY_FIXED", "convention": "SIDE_LOOK", "sideLookAngle":25})) # default orientation
        
        h = 500 # [km]
        orb_speed = orbital_speed(h)

        inc_deg = 60
        obsv_metrics = test_sar2.calc_data_metrics(alt_km=h, sc_speed_kmps=orb_speed, sc_gnd_speed_kmps=orb_speed*(RE/(RE+h)), inc_angle_deg=inc_deg, 
                                                    instru_look_angle_from_target_inc_angle=True)
        self.assertEqual(obsv_metrics, {'ground pixel along-track resolution [m]': 6.4, 'ground pixel cross-track resolution [m]': 115.39, 'NESZ [dB]': -27.47, 'incidence angle [deg]': 60.0, 'swath-width [km]': 68.0, 'PRF [Hz]': 2382})

        inc_deg = 35
        obsv_metrics = test_sar2.calc_data_metrics(alt_km=h, sc_speed_kmps=orb_speed, sc_gnd_speed_kmps=orb_speed*(RE/(RE+h)), inc_angle_deg=inc_deg, 
                                                    instru_look_angle_from_target_inc_angle=True)
        self.assertEqual(obsv_metrics, {'ground pixel along-track resolution [m]': 6.4, 'ground pixel cross-track resolution [m]': 174.22, 'NESZ [dB]': -39.3, 'incidence angle [deg]': 35.0, 'swath-width [km]': 27.3, 'PRF [Hz]': 6902})

        inc_deg = 45
        obsv_metrics = test_sar2.calc_data_metrics(alt_km=h, sc_speed_kmps=orb_speed, sc_gnd_speed_kmps=orb_speed*(RE/(RE+h)), inc_angle_deg=inc_deg, 
                                                    instru_look_angle_from_target_inc_angle=True)
        self.assertEqual(obsv_metrics,{'ground pixel along-track resolution [m]': 6.4, 'ground pixel cross-track resolution [m]': 141.32, 'NESZ [dB]': -34.86, 'incidence angle [deg]': 45.0, 'swath-width [km]': 36.1, 'PRF [Hz]': 4520})

        inc_deg = 55
        obsv_metrics = test_sar2.calc_data_metrics(alt_km=h, sc_speed_kmps=orb_speed, sc_gnd_speed_kmps=orb_speed*(RE/(RE+h)), inc_angle_deg=inc_deg, 
                                                    instru_look_angle_from_target_inc_angle=True)
        self.assertEqual(obsv_metrics,{'ground pixel along-track resolution [m]': 6.4, 'ground pixel cross-track resolution [m]': 121.99, 'NESZ [dB]': -29.62, 'incidence angle [deg]': 55.0, 'swath-width [km]': 53.1, 'PRF [Hz]': 2662})
        
    def test_from_json_calc_data_metrics_test_sar3(self):
        """ Test of SAR with specs as given by ``test_sar3_json_str``. Following aspects of the SAR are in focus:
            
            * ALtitude = 500km.
            * Different incidence angles: 60 deg, 35 deg, 45 deg, 55 deg.
            * L-band.
            * Flag ``instru_look_angle_from_target_inc_angle`` is set to ``True``.
            * Full-swath.
            * Dual-pol, SMAP.

            Tests:

            * Initialization (partly).
            * Resulting datametrics are compared with precomputed values (generated on 2 April 2021). 

        """
        
        test_sar3_json_str = '{"@type": "Synthetic Aperture Radar",' \
                    '"pulseWidth": 1.64e-6,' \
                    '"antenna":{"shape": "RECTANGULAR", "height": 3.98, "width": 11.03, "apertureEfficiency": 0.6, "apertureExcitationProfile": "UNIFORM"},' \
                    '"operatingFrequency": 1275.7e6,' \
                    '"peakTransmitPower": 1000,' \
                    '"chirpBandwidth": 1.290e6,' \
                    '"minimumPRF": 1,' \
                    '"maximumPRF": 20000,' \
                    '"radarLoss": 2,' \
                    '"systemNoiseFigure": 2,' \
                    '"swathConfig": {' \
                    '   "@type": "full"' \
                    '},' \
                    '"polarization": {' \
                    '   "@type": "dual",' \
                    '   "pulseConfig": {' \
                    '        "@type": "SMAP"' \
                    '}' \
                    '}' \
                    '}'
        
        test_sar3 = SyntheticApertureRadarModel.from_json(test_sar3_json_str)
        
        self.assertIsInstance(test_sar3.swathType, SwathTypeSAR)
        self.assertEqual(test_sar3.swathType, SwathTypeSAR.FULL)
        self.assertIsInstance(test_sar3.polType, PolTypeSAR)
        self.assertEqual(test_sar3.polType, PolTypeSAR.DUAL)
        self.assertIsInstance(test_sar3.dualPolPulseConfig, DualPolPulseConfig)
        self.assertEqual(test_sar3.dualPolPulseConfig, DualPolPulseConfig.SMAP)
        self.assertIsInstance(test_sar3.dualPolPulseSep, float)
        self.assertEqual(test_sar3.dualPolPulseSep, 0.5*test_sar3.pulseWidth) # default value
        self.assertIsInstance(test_sar3.orientation, Orientation)
        self.assertEqual(test_sar3.orientation, Orientation.from_dict({"referenceFrame": "SC_BODY_FIXED", "convention": "SIDE_LOOK", "sideLookAngle":25})) # default orientation
        
        h = 500 # [km]
        orb_speed = orbital_speed(h)

        inc_deg = 60
        obsv_metrics = test_sar3.calc_data_metrics(alt_km=h, sc_speed_kmps=orb_speed, sc_gnd_speed_kmps=orb_speed*(RE/(RE+h)), inc_angle_deg=inc_deg, 
                                                    instru_look_angle_from_target_inc_angle=True)
        self.assertEqual(obsv_metrics, {'ground pixel along-track resolution [m]': 1.85, 'ground pixel cross-track resolution [m]': 161.01, 'NESZ [dB]': -30.6, 'incidence angle [deg]': 60.0, 'swath-width [km]': 38.8, 'PRF [Hz]': 4202})

        inc_deg = 35
        obsv_metrics = test_sar3.calc_data_metrics(alt_km=h, sc_speed_kmps=orb_speed, sc_gnd_speed_kmps=orb_speed*(RE/(RE+h)), inc_angle_deg=inc_deg, 
                                                    instru_look_angle_from_target_inc_angle=True)
        self.assertEqual(obsv_metrics, {'ground pixel along-track resolution [m]': 1.85, 'ground pixel cross-track resolution [m]': 243.1, 'NESZ [dB]': -42.15, 'incidence angle [deg]': 35.0, 'swath-width [km]': 15.6, 'PRF [Hz]': 11396})

        inc_deg = 45
        obsv_metrics = test_sar3.calc_data_metrics(alt_km=h, sc_speed_kmps=orb_speed, sc_gnd_speed_kmps=orb_speed*(RE/(RE+h)), inc_angle_deg=inc_deg, 
                                                    instru_look_angle_from_target_inc_angle=True)
        self.assertEqual(obsv_metrics,{'ground pixel along-track resolution [m]': 1.85, 'ground pixel cross-track resolution [m]': 197.2, 'NESZ [dB]': -37.91, 'incidence angle [deg]': 45.0, 'swath-width [km]': 20.6, 'PRF [Hz]': 7808})

        inc_deg = 55
        obsv_metrics = test_sar3.calc_data_metrics(alt_km=h, sc_speed_kmps=orb_speed, sc_gnd_speed_kmps=orb_speed*(RE/(RE+h)), inc_angle_deg=inc_deg, 
                                                    instru_look_angle_from_target_inc_angle=True)
        self.assertEqual(obsv_metrics,{'ground pixel along-track resolution [m]': 1.85, 'ground pixel cross-track resolution [m]': 170.22, 'NESZ [dB]': -32.59, 'incidence angle [deg]': 55.0, 'swath-width [km]': 30.3, 'PRF [Hz]': 4523})

    def test_from_json_calc_data_metrics_test_sar4(self):
        """ Test of SAR with specs as given by ``test_sar4_json_str``. Following aspects of the SAR are in focus:
            
            * ALtitude = 500km.
            * Different incidence angles: 60 deg, 35 deg, 45 deg, 55 deg.
            * P-band.
            * Flag ``instru_look_angle_from_target_inc_angle`` is set to ``True``.
            * Fixed-swath 50km.
            * Dual-pol, SMAP.

            Tests:

            * Initialization (partly).
            * Resulting datametrics are compared with precomputed values (generated on 2 April 2021). 

        """        
        test_sar4_json_str = '{"@type": "Synthetic Aperture Radar",' \
                    '"pulseWidth": 8.81e-6,' \
                    '"antenna":{"shape": "RECTANGULAR", "height": 19.07, "width": 2.86, "apertureEfficiency": 0.6, "apertureExcitationProfile": "UNIFORM"},' \
                    '"operatingFrequency": 435e6,' \
                    '"peakTransmitPower": 1000,' \
                    '"chirpBandwidth": 0.616e6,' \
                    '"minimumPRF": 1,' \
                    '"maximumPRF": 20000,' \
                    '"radarLoss": 2,' \
                    '"systemNoiseFigure": 2,' \
                    '"swathConfig": {' \
                    '   "@type": "fixed",' \
                    '   "fixedSwathSize": 50' \
                    '},' \
                    '"polarization": {' \
                    '   "@type": "dual",' \
                    '   "pulseConfig": {' \
                    '        "@type": "SMAP"' \
                    '}' \
                    '},' \
                    '"scanTechnique": "STripmap"' \
                    '}'
        
        test_sar4 = SyntheticApertureRadarModel.from_json(test_sar4_json_str)
        
        self.assertIsInstance(test_sar4.swathType, SwathTypeSAR)
        self.assertEqual(test_sar4.swathType, SwathTypeSAR.FIXED)
        self.assertIsInstance(test_sar4.fixedSwathSize, float)
        self.assertEqual(test_sar4.fixedSwathSize, 50)
        self.assertIsInstance(test_sar4.polType, PolTypeSAR)
        self.assertEqual(test_sar4.polType, PolTypeSAR.DUAL)
        self.assertIsInstance(test_sar4.dualPolPulseConfig, DualPolPulseConfig)
        self.assertEqual(test_sar4.dualPolPulseConfig, DualPolPulseConfig.SMAP)
        self.assertIsInstance(test_sar4.dualPolPulseSep, float)
        self.assertEqual(test_sar4.dualPolPulseSep, 0.5*test_sar4.pulseWidth) # default value
        self.assertIsInstance(test_sar4.orientation, Orientation)
        self.assertEqual(test_sar4.orientation, Orientation.from_dict({"referenceFrame": "SC_BODY_FIXED", "convention": "SIDE_LOOK", "sideLookAngle":25})) # default orientation
        self.assertIsInstance(test_sar4.scanTechnique, ScanTech)
        self.assertEqual(test_sar4.scanTechnique, ScanTech.STRIPMAP)

        h = 500 # [km]
        orb_speed = orbital_speed(h)

        inc_deg = 60
        obsv_metrics = test_sar4.calc_data_metrics(alt_km=h, sc_speed_kmps=orb_speed, sc_gnd_speed_kmps=orb_speed*(RE/(RE+h)), inc_angle_deg=inc_deg, 
                                                    instru_look_angle_from_target_inc_angle=True)
        self.assertEqual(obsv_metrics, {'ground pixel along-track resolution [m]': 8.84, 'ground pixel cross-track resolution [m]': 337.18, 'NESZ [dB]': -31.47, 'incidence angle [deg]': 60.0, 'swath-width [km]': 50.0, 'PRF [Hz]': 866})
        
        inc_deg = 35
        obsv_metrics = test_sar4.calc_data_metrics(alt_km=h, sc_speed_kmps=orb_speed, sc_gnd_speed_kmps=orb_speed*(RE/(RE+h)), inc_angle_deg=inc_deg, 
                                                    instru_look_angle_from_target_inc_angle=True)
        self.assertEqual(obsv_metrics, {'ground pixel along-track resolution [m]': 8.84, 'ground pixel cross-track resolution [m]': 509.1, 'NESZ [dB]': -43.15, 'incidence angle [deg]': 35.0, 'swath-width [km]': 50.0, 'PRF [Hz]': 2425})
        
        inc_deg = 45
        obsv_metrics = test_sar4.calc_data_metrics(alt_km=h, sc_speed_kmps=orb_speed, sc_gnd_speed_kmps=orb_speed*(RE/(RE+h)), inc_angle_deg=inc_deg, 
                                                    instru_look_angle_from_target_inc_angle=True)
        self.assertEqual(obsv_metrics,{'ground pixel along-track resolution [m]': 8.84, 'ground pixel cross-track resolution [m]': 412.96, 'NESZ [dB]': -38.37, 'incidence angle [deg]': 45.0, 'swath-width [km]': 50.0, 'PRF [Hz]': 1467})
        
        inc_deg = 55
        obsv_metrics = test_sar4.calc_data_metrics(alt_km=h, sc_speed_kmps=orb_speed, sc_gnd_speed_kmps=orb_speed*(RE/(RE+h)), inc_angle_deg=inc_deg, 
                                                    instru_look_angle_from_target_inc_angle=True)
        self.assertEqual(obsv_metrics,{'ground pixel along-track resolution [m]': 8.84, 'ground pixel cross-track resolution [m]': 356.47, 'NESZ [dB]': -34.06, 'incidence angle [deg]': 55.0, 'swath-width [km]': 50.0, 'PRF [Hz]': 1071})

    def test_from_json_calc_data_metrics_test_sar5(self):
        """ Test of SAR with specs as given by ``test_sar5_json_str``. Following aspects of the SAR are in focus:
            
            * ALtitude = 600km.
            * Different incidence angles: 60 deg, 35 deg, 45 deg, 55 deg.
            * L-band.
            * Flag ``instru_look_angle_from_target_inc_angle`` is set to ``True``.
            * Fixed-swath 25km.
            * Dual-pol, AIRSAR.

            Tests:

            * Initialization (partly).
            * Resulting datametrics are compared with precomputed values (generated on 2 April 2021). 

        """        
        test_sar5_json_str = '{"@type": "Synthetic Aperture Radar",' \
                    '"pulseWidth": 3.73e-6,' \
                    '"antenna":{"shape": "RECTANGULAR", "height": 13.57, "width": 4.15, "apertureEfficiency": 0.6, "apertureExcitationProfile": "UNIFORM"},' \
                    '"operatingFrequency": 1257.7e6,' \
                    '"peakTransmitPower": 1000,' \
                    '"chirpBandwidth": 0.75e6,' \
                    '"minimumPRF": 1,' \
                    '"maximumPRF": 20000,' \
                    '"radarLoss": 2,' \
                    '"systemNoiseFigure": 2,' \
                    '"swathConfig": {' \
                    '   "@type": "fixed",' \
                    '   "fixedSwathSize": 25' \
                    '},' \
                    '"polarization": {' \
                    '   "@type": "dual",' \
                    '   "pulseConfig": {' \
                    '        "@type": "AIRSAR"' \
                    '}' \
                    '},' \
                    '"scanTechnique": "STripmap"' \
                    '}'
        
        test_sar5 = SyntheticApertureRadarModel.from_json(test_sar5_json_str)
        
        self.assertIsInstance(test_sar5.swathType, SwathTypeSAR)
        self.assertEqual(test_sar5.swathType, SwathTypeSAR.FIXED)
        self.assertIsInstance(test_sar5.fixedSwathSize, float)
        self.assertEqual(test_sar5.fixedSwathSize, 25)
        self.assertIsInstance(test_sar5.polType, PolTypeSAR)
        self.assertEqual(test_sar5.polType, PolTypeSAR.DUAL)
        self.assertIsInstance(test_sar5.dualPolPulseConfig, DualPolPulseConfig)
        self.assertEqual(test_sar5.dualPolPulseConfig, DualPolPulseConfig.AIRSAR)
        self.assertIsNone(test_sar5.dualPolPulseSep)
        self.assertIsInstance(test_sar5.orientation, Orientation)
        self.assertEqual(test_sar5.orientation, Orientation.from_dict({"referenceFrame": "SC_BODY_FIXED", "convention": "SIDE_LOOK", "sideLookAngle":25})) # default orientation
        self.assertIsInstance(test_sar5.scanTechnique, ScanTech)
        self.assertEqual(test_sar5.scanTechnique, ScanTech.STRIPMAP)

        h = 600 # [km]
        orb_speed = orbital_speed(h)

        inc_deg = 60
        obsv_metrics = test_sar5.calc_data_metrics(alt_km=h, sc_speed_kmps=orb_speed, sc_gnd_speed_kmps=orb_speed*(RE/(RE+h)), inc_angle_deg=inc_deg, 
                                                    instru_look_angle_from_target_inc_angle=True)
        self.assertEqual(obsv_metrics, {'ground pixel along-track resolution [m]': 6.2, 'ground pixel cross-track resolution [m]': 276.94, 'NESZ [dB]': -31.11, 'incidence angle [deg]': 60.0, 'swath-width [km]': 25.0, 'PRF [Hz]': 2439})

        inc_deg = 35
        obsv_metrics = test_sar5.calc_data_metrics(alt_km=h, sc_speed_kmps=orb_speed, sc_gnd_speed_kmps=orb_speed*(RE/(RE+h)), inc_angle_deg=inc_deg, 
                                                    instru_look_angle_from_target_inc_angle=True)
        self.assertEqual(obsv_metrics, {'ground pixel along-track resolution [m]': 6.2, 'ground pixel cross-track resolution [m]': 418.14, 'NESZ [dB]': -42.63, 'incidence angle [deg]': 35.0, 'swath-width [km]': 25.0, 'PRF [Hz]': 6818})
        
        inc_deg = 45
        obsv_metrics = test_sar5.calc_data_metrics(alt_km=h, sc_speed_kmps=orb_speed, sc_gnd_speed_kmps=orb_speed*(RE/(RE+h)), inc_angle_deg=inc_deg, 
                                                    instru_look_angle_from_target_inc_angle=True)
        self.assertEqual(obsv_metrics,{'ground pixel along-track resolution [m]': 6.2, 'ground pixel cross-track resolution [m]': 339.18, 'NESZ [dB]': -38.43, 'incidence angle [deg]': 45.0, 'swath-width [km]': 25.0, 'PRF [Hz]': 4678})

        inc_deg = 55
        obsv_metrics = test_sar5.calc_data_metrics(alt_km=h, sc_speed_kmps=orb_speed, sc_gnd_speed_kmps=orb_speed*(RE/(RE+h)), inc_angle_deg=inc_deg, 
                                                    instru_look_angle_from_target_inc_angle=True)
        self.assertEqual(obsv_metrics,{'ground pixel along-track resolution [m]': 6.2, 'ground pixel cross-track resolution [m]': 292.78, 'NESZ [dB]': -33.72, 'incidence angle [deg]': 55.0, 'swath-width [km]': 25.0, 'PRF [Hz]': 3064})

    def test_from_dict_calc_data_metrics_scansar(self):
        """ Test ScanSAR operation. Following aspects of the SAR are in focus:
            
            * ScanSAR, numSubSwaths = 1, 2, 5. When  numSubSwaths=1, the result is the same as from Stripmap.
            * ALtitude = 600km.
            * Different incidence angles: 15 deg and 30 deg.
            * Flag ``instru_look_angle_from_target_inc_angle`` is set to ``True``.
            * Fixed-swath 25km.            
            * Dual-pol, SMAP, custom pulse-seperation
            * Custom atmospheric loss

            Tests:

            * Initialization (partly).
            * Resulting datametrics are compared with precomputed values (generated on 2 April 2021). 

        """                
        h = 600 # [km]
        orb_speed = orbital_speed(h)

        test_sar_dict = {"@type": "Synthetic Aperture Radar",
                    "pulseWidth": 12e-6,
                    "antenna":{"shape": "RECTANGULAR", "height": 4.9, "width": 1, "apertureEfficiency": 0.6, "apertureExcitationProfile": "UNIFORM"},
                    "operatingFrequency": 9.65e9,
                    "peakTransmitPower": 1000,
                    "chirpBandwidth": 75e6,
                    "minimumPRF": 1,
                    "maximumPRF": 20000,
                    "radarLoss": 3.5,
                    "systemNoiseFigure": 4.5,
                    "polarization": {
                       "@type": "dual",
                       "pulseConfig": {
                            "@type": "SMAP",
                            "pulseSeparation": 2e-6
                    }
                    },
                    "scanTechnique": "Scansar",
                    "numSubSwaths": 1,
                    "atmosLoss": 1
                    }

        ######## numSubSwaths = 1 ########
        test_sar = SyntheticApertureRadarModel.from_dict(test_sar_dict)
        
        self.assertIsInstance(test_sar.scanTechnique, ScanTech)
        self.assertEqual(test_sar.scanTechnique, ScanTech.SCANSAR)
        self.assertIsInstance(test_sar.numSubSwaths, int)
        self.assertEqual(test_sar.numSubSwaths, 1)
        self.assertIsInstance(test_sar.swathType, SwathTypeSAR)
        self.assertEqual(test_sar.swathType, SwathTypeSAR.FULL)
        self.assertIsInstance(test_sar.polType, PolTypeSAR)
        self.assertEqual(test_sar.polType, PolTypeSAR.DUAL)
        self.assertIsInstance(test_sar.dualPolPulseConfig, DualPolPulseConfig)
        self.assertEqual(test_sar.dualPolPulseConfig, DualPolPulseConfig.SMAP)
        self.assertIsInstance(test_sar.dualPolPulseSep, float)
        self.assertEqual(test_sar.dualPolPulseSep, 2e-6)
        self.assertIsInstance(test_sar.orientation, Orientation)
        self.assertEqual(test_sar.orientation, Orientation.from_dict({"referenceFrame": "SC_BODY_FIXED", "convention": "SIDE_LOOK", "sideLookAngle":25})) # default orientation    
        self.assertIsInstance(test_sar.atmosLoss, float)
        self.assertEqual(test_sar.atmosLoss, 1)     

        inc_deg = 15
        obsv_metrics = test_sar.calc_data_metrics(alt_km=h, sc_speed_kmps=orb_speed, sc_gnd_speed_kmps=orb_speed*(RE/(RE+h)), inc_angle_deg=inc_deg, 
                                                    instru_look_angle_from_target_inc_angle=True)
        self.assertEqual(obsv_metrics, {'ground pixel along-track resolution [m]': 2.24, 'ground pixel cross-track resolution [m]': 9.27, 'NESZ [dB]': -19.99, 'incidence angle [deg]': 15.0, 'swath-width [km]': 19.9, 'PRF [Hz]': 5748})

        inc_deg = 30
        obsv_metrics = test_sar.calc_data_metrics(alt_km=h, sc_speed_kmps=orb_speed, sc_gnd_speed_kmps=orb_speed*(RE/(RE+h)), inc_angle_deg=inc_deg, 
                                                    instru_look_angle_from_target_inc_angle=True)
        self.assertEqual(obsv_metrics,{'ground pixel along-track resolution [m]': 2.24, 'ground pixel cross-track resolution [m]': 4.8, 'NESZ [dB]': -16.23, 'incidence angle [deg]': 30.0, 'swath-width [km]': 24.5, 'PRF [Hz]': 6269})


        ######## numSubSwaths = 2 ########
        test_sar_dict['numSubSwaths'] = 2
        test_sar = SyntheticApertureRadarModel.from_dict(test_sar_dict)
        
        self.assertIsInstance(test_sar.numSubSwaths, int)
        self.assertEqual(test_sar.numSubSwaths, 2)

        inc_deg = 15
        obsv_metrics = test_sar.calc_data_metrics(alt_km=h, sc_speed_kmps=orb_speed, sc_gnd_speed_kmps=orb_speed*(RE/(RE+h)), inc_angle_deg=inc_deg, 
                                                    instru_look_angle_from_target_inc_angle=True)
        self.assertEqual(obsv_metrics, {'ground pixel along-track resolution [m]': 4.48, 'ground pixel cross-track resolution [m]': 9.27, 'NESZ [dB]': -22.5, 'incidence angle [deg]': 15.0, 'swath-width [km]': 39.9, 'PRF [Hz]': 10252})

        inc_deg = 30
        obsv_metrics = test_sar.calc_data_metrics(alt_km=h, sc_speed_kmps=orb_speed, sc_gnd_speed_kmps=orb_speed*(RE/(RE+h)), inc_angle_deg=inc_deg, 
                                                    instru_look_angle_from_target_inc_angle=True)
        self.assertEqual(obsv_metrics,{'ground pixel along-track resolution [m]': 4.48, 'ground pixel cross-track resolution [m]': 4.8, 'NESZ [dB]': -16.04, 'incidence angle [deg]': 30.0, 'swath-width [km]': 49.1, 'PRF [Hz]': 5996})

        ######## numSubSwaths = 5 ########
        test_sar_dict['numSubSwaths'] = 5
        test_sar = SyntheticApertureRadarModel.from_dict(test_sar_dict)
        
        self.assertIsInstance(test_sar.numSubSwaths, int)
        self.assertEqual(test_sar.numSubSwaths, 5)

        inc_deg = 15
        obsv_metrics = test_sar.calc_data_metrics(alt_km=h, sc_speed_kmps=orb_speed, sc_gnd_speed_kmps=orb_speed*(RE/(RE+h)), inc_angle_deg=inc_deg, 
                                                    instru_look_angle_from_target_inc_angle=True)
        self.assertEqual(obsv_metrics, {'ground pixel along-track resolution [m]': 11.2, 'ground pixel cross-track resolution [m]': 9.27, 'NESZ [dB]': -21.29, 'incidence angle [deg]': 15.0, 'swath-width [km]': 99.9, 'PRF [Hz]': 7747})

        inc_deg = 30
        obsv_metrics = test_sar.calc_data_metrics(alt_km=h, sc_speed_kmps=orb_speed, sc_gnd_speed_kmps=orb_speed*(RE/(RE+h)), inc_angle_deg=inc_deg, 
                                                    instru_look_angle_from_target_inc_angle=True)
        self.assertEqual(obsv_metrics,{'ground pixel along-track resolution [m]': 11.2, 'ground pixel cross-track resolution [m]': 4.8, 'NESZ [dB]': -15.05, 'incidence angle [deg]': 30.0, 'swath-width [km]': 123.1, 'PRF [Hz]': 4776})

    def test_get_id(self): #@TODO
        pass

    def test_get_field_of_view(self): #@TODO
        pass

    def test_get_scene_field_of_view(self): #@TODO
        pass

    def test_get_field_of_regard(self): #@TODO
        pass

    def test_get_orientation(self): #@TODO
        pass

    def test_get_pointing_option(self): #TODO
        pass
    
    def test_to_dict(self): #@TODO
        pass

