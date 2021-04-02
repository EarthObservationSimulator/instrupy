"""Unit tests for instrupy.synthetic_aperture_radar.synthetic_aperture_radar_model.
"""

import unittest
import json
import numpy as np
import sys, os


from instrupy.synthetic_aperture_radar_model import SyntheticApertureRadarModel, ScanTechSAR, PolTypeSAR, DualPolPulseConfig, SwathTypeSAR
from instrupy.util import Orientation, SphericalGeometry, ViewGeometry, FileUtilityFunctions, Maneuver

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
                        '    "convention": "SIDE_LOOK",' \
                        '    "sideLookAngle": 30' \
                        '},' \
                        '"dataRate": 2000,' \
                        '"bitsPerPixel": 16,' \
                        '"pulseWidth": 31e-6,' \
                        '"antennaHeight": 4.9,' \
                        '"antennaWidth": 0.7,' \
                        '"antennaApertureEfficiency": 0.5,' \
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
                    '    "convention": "SIDE_LOOK",' \
                    '    "sideLookAngle": 20.5' \
                    '},' \
                    '"pulseWidth": 33.4e-6,' \
                    '"antennaHeight": 10.7,' \
                    '"antennaWidth": 2.16,' \
                    '"antennaApertureEfficiency": 0.6,' \
                    '"operatingFrequency": 1.2757e9,' \
                    '"peakTransmitPower": 1000,' \
                    '"chirpBandwidth": 19e6,' \
                    '"minimumPRF": 1463,' \
                    '"maximumPRF": 1640,' \
                    '"radarLoss": 3.5,' \
                    '"systemNoiseFigure": 5.11' \
                    '}'

# ERS1
ers1_json_str = '{"@type": "Synthetic Aperture Radar",' \
                '   "orientation": {' \
                '      "convention": "SIDE_LOOK",' \
                '      "sideLookAngle": 20' \
                '  },' \
                '  "pulseWidth": 37.1e-6,' \
                '  "antennaHeight": 10,' \
                '  "antennaWidth": 1,' \
                '  "antennaApertureEfficiency": 0.26,' \
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
        """ Test initialization of the synthetic aperture radar in the many different ways allowed. 
        """   
        microxsar= SyntheticApertureRadarModel.from_json(microxsar_json_str) 
        # SAR without specification of maneuverability, numStripsInScene, altitude fields
        self.assertEqual(microxsar._type, "Synthetic Aperture Radar")
        self.assertIsInstance(microxsar.name, str)
        self.assertEqual(microxsar.name, "MiroXSAR")
        self.assertIsInstance(microxsar.mass, float)
        self.assertEqual(microxsar.mass, 130)
        self.assertIsInstance(microxsar.volume, float)
        self.assertEqual(microxsar.volume, 0.343)
        self.assertIsInstance(microxsar.power, float)
        self.assertEqual(microxsar.power, 1100)
        self.assertEqual(microxsar.orientation, Orientation.from_json({"convention":"Euler", "eulerSeq1":1, "eulerSeq2":2, "eulerSeq3":3, "eulerAngle1":0, "eulerAngle2":30, "eulerAngle3":0}))
        self.assertIsInstance(microxsar.orientation, Orientation)
        self.assertIsInstance(microxsar.fieldOfView, ViewGeometry)         
        self.assertAlmostEqual(microxsar.fieldOfView.sph_geom.get_fov_height_and_width()[0], 0.3635, places = 2)
        self.assertAlmostEqual(microxsar.fieldOfView.sph_geom.get_fov_height_and_width()[1], 2.5446, places = 2)
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
        self.assertIsInstance(microxsar.antennaHeight, float)
        self.assertEqual(microxsar.antennaHeight, 4.9)
        self.assertIsInstance(microxsar.antennaWidth, float)
        self.assertEqual(microxsar.antennaWidth, 0.7)
        self.assertIsInstance(microxsar.antennaApertureEfficiency, float)
        self.assertEqual(microxsar.antennaApertureEfficiency, 0.5)
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
        self.assertIsInstance(microxsar.scanTechnique, ScanTechSAR) 
        self.assertEqual(microxsar.scanTechnique, ScanTechSAR.STRIPMAP) 
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
                                                  '   "convention": "Euler",'
                                                  '   "eulerAngle1": 0,'
                                                  '   "eulerAngle2": 10,'
                                                  '   "eulerAngle3": 0'
                                                  ' },'
                                                  '"dataRate": 2000,'
                                                  '"bitsPerPixel": 16,'
                                                  '"pulseWidth": 31e-6,'
                                                  '"antennaHeight": 4.9,'
                                                  '"antennaWidth": 0.7,' 
                                                  '"antennaApertureEfficiency": 0.5,' 
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
                                                  '    "convention": "SIDE_LOOK",'
                                                  '    "sideLookAngle": 30'
                                                  '},'
                                                  '"dataRate": 2000,'
                                                  '"bitsPerPixel": 16,'
                                                  '"pulseWidth": 31e-6,'
                                                  '"antennaHeight": 4.9,'
                                                  '"antennaWidth": 0.7,' 
                                                  '"antennaApertureEfficiency": 0.5,' 
                                                  '"operatingFrequency": 9.65e9,' 
                                                  '"peakTransmitPower": 1e3,' 
                                                  '"chirpBandwidth": 75e6,'      
                                                  '"minimumPRF": 8000,' 
                                                  '"maximumPRF": 3000,' 
                                                  '"radarLosses": 3.5,' 
                                                  '"sceneNoiseTemp": 290,' 
                                                  '"systemNoiseFigure": 4.3,'
                                                  '}')


    '''   
    def test_from_json_basic_2(self):
        """ Test initialization of the synthetic aperture radar in the many different ways allowed.
        """
        # SAR with maneuverability, numStripsInScene, altitude fields and without other optional fields
        self.assertEqual(self.seasat._type, "Synthetic Aperture Radar")
        self.assertIsNone(self.seasat.name)
        self.assertIsNone(self.seasat.acronym)
        self.assertIsNone(self.seasat.mass)
        self.assertIsNone(self.seasat.volume)
        self.assertIsNone(self.seasat.power)
        self.assertIsNone(self.seasat.dataRate)
        self.assertIsNone(self.seasat.bitsPerPixel)
        self.assertIsInstance(self.seasat.orientation, Orientation)
        self.assertEqual(self.seasat.orientation.euler_seq1, 1)
        self.assertEqual(self.seasat.orientation.euler_seq2, 2)
        self.assertEqual(self.seasat.orientation.euler_seq3, 3)
        self.assertEqual(self.seasat.orientation.euler_angle1, 0)
        self.assertEqual(self.seasat.orientation.euler_angle2, 20.5)
        self.assertEqual(self.seasat.orientation.euler_angle3, 0)
        self.assertEqual(self.seasat.pulseWidth, 33.4e-6)
        self.assertIsInstance(self.seasat.pulseWidth, float)
        self.assertEqual(self.seasat.antennaHeight, 10.7)
        self.assertIsInstance(self.seasat.antennaHeight, float)
        self.assertEqual(self.seasat.antennaWidth, 2.16)
        self.assertIsInstance(self.seasat.antennaWidth, float)
        self.assertEqual(self.seasat.antennaApertureEfficiency, 0.6)
        self.assertIsInstance(self.seasat.antennaApertureEfficiency, float)
        self.assertEqual(self.seasat.operatingFrequency, 1.2757e9)
        self.assertIsInstance(self.seasat.operatingFrequency, float)
        self.assertEqual(self.seasat.peakTransmitPower, 1e3)
        self.assertIsInstance(self.seasat.peakTransmitPower, float)
        self.assertEqual(self.seasat.chirpBandwidth, 19e6)
        self.assertIsInstance(self.seasat.chirpBandwidth, float)
        self.assertEqual(self.seasat.minimumPRF, 1463)
        self.assertIsInstance(self.seasat.minimumPRF, float)
        self.assertEqual(self.seasat.maximumPRF, 1640)
        self.assertIsInstance(self.seasat.maximumPRF, float)
        self.assertEqual(self.seasat.radarLosses, 3.5)
        self.assertIsInstance(self.seasat.radarLosses, float)
        self.assertEqual(self.seasat.sceneNoiseTemp, 290) # default value
        self.assertIsInstance(self.seasat.sceneNoiseTemp, float)
        self.assertEqual(self.seasat.systemNoiseFigure, 5.11)
        self.assertIsInstance(self.seasat.systemNoiseFigure, float)
        # along-track fov and cross-track fov are calculated from antenna size and initialized in the object
        self.assertAlmostEqual(self.seasat.fieldOfView.get_rectangular_fov_specs_from_custom_fov_specs()[0], 1.2584, places = 2)
        self.assertAlmostEqual(self.seasat.fieldOfView.get_rectangular_fov_specs_from_custom_fov_specs()[1], 6.2336, places = 2)
        # sceneFOV
        self.assertAlmostEqual(self.seasat.sceneFieldOfView._AT_fov_deg, 5.5145, places = 2)
        self.assertAlmostEqual(self.seasat.sceneFieldOfView._CT_fov_deg, 6.2336, places = 2)
        # FOR is built from sceneFOV
        self.assertAlmostEqual(self.seasat.fieldOfRegard.get_ATCT_fov()[0], 5.5145, places = 2)
        self.assertAlmostEqual(self.seasat.fieldOfRegard.get_ATCT_fov()[1], 20+6.2336, places = 2)
        self.assertTrue(self.seasat.fieldOfRegard._yaw180_flag)

    def test_from_json_basic_3(self):
        """ Test initialization of the synthetic aperture radar in the many different ways allowed.
        """
         # SAR with maneuverability field and without other optional fields
          # SAR with maneuverability, numStripsInScene, altitude fields and without other optional fields
        self.assertEqual(self.ers1._type, "Synthetic Aperture Radar")
        self.assertIsNone(self.ers1.name)
        self.assertIsNone(self.ers1.acronym)
        self.assertIsNone(self.ers1.mass)
        self.assertIsNone(self.ers1.volume)
        self.assertIsNone(self.ers1.power)
        self.assertIsNone(self.ers1.dataRate)
        self.assertIsNone(self.ers1.bitsPerPixel)
        self.assertIsInstance(self.ers1.orientation, Orientation)
        self.assertEqual(self.ers1.orientation.euler_seq1, 1)
        self.assertEqual(self.ers1.orientation.euler_seq2, 2)
        self.assertEqual(self.ers1.orientation.euler_seq3, 3)
        self.assertEqual(self.ers1.orientation.euler_angle1, 0)
        self.assertEqual(self.ers1.orientation.euler_angle2, 20)
        self.assertEqual(self.ers1.orientation.euler_angle3, 0)
        self.assertEqual(self.ers1.pulseWidth, 0.0000371)
        self.assertIsInstance(self.ers1.pulseWidth, float)
        self.assertEqual(self.ers1.antennaHeight, 10)
        self.assertIsInstance(self.ers1.antennaHeight, float)
        self.assertEqual(self.ers1.antennaWidth, 1)
        self.assertIsInstance(self.ers1.antennaWidth, float)
        self.assertEqual(self.ers1.antennaApertureEfficiency, 0.26)
        self.assertIsInstance(self.ers1.antennaApertureEfficiency, float)
        self.assertEqual(self.ers1.operatingFrequency, 5.25e9)
        self.assertIsInstance(self.ers1.operatingFrequency, float)
        self.assertEqual(self.ers1.peakTransmitPower, 4.8e3)
        self.assertIsInstance(self.ers1.peakTransmitPower, float)
        self.assertEqual(self.ers1.chirpBandwidth, 15500000)
        self.assertIsInstance(self.ers1.chirpBandwidth, float)
        self.assertEqual(self.ers1.minimumPRF, 1680)
        self.assertIsInstance(self.ers1.minimumPRF, float)
        self.assertEqual(self.ers1.maximumPRF, 1700)
        self.assertIsInstance(self.ers1.maximumPRF, float)
        self.assertEqual(self.ers1.radarLosses, 3.5)
        self.assertIsInstance(self.ers1.radarLosses, float)
        self.assertEqual(self.ers1.systemNoiseFigure, 3.4)
        self.assertIsInstance(self.ers1.systemNoiseFigure, float)
        # along-track fov and cross-track fov are calculated from antenna size and initialized in the object
        self.assertAlmostEqual(self.ers1.fieldOfView.get_ATCT_fov()[0], 0.3272, places = 2)
        self.assertAlmostEqual(self.ers1.fieldOfView.get_ATCT_fov()[1], 3.2718, places = 2)
        # sceneFOV
        self.assertIsNone(self.ers1.sceneFieldOfView)
        # FOR is built from FOV
        self.assertAlmostEqual(self.ers1.fieldOfRegard.get_rectangular_fov_specs_from_custom_fov_specs()[0], 0.3272, places = 2)
        self.assertAlmostEqual(self.ers1.fieldOfRegard.get_rectangular_fov_specs_from_custom_fov_specs()[1], 30+3.2718, places = 2)
        self.assertFalse(self.ers1.fieldOfRegard._yaw180_flag)

    def test_find_valid_highest_possible_PRF(self):
        """ Truth test data from Spaceborne SAR Study: LDRD 92 Final Report SANDIA Report March 1993, Page 35
            Using look angle = 18.5 deg corresponding to incidence angle of 20 deg.
            Choose Delev and fc to get a swath of 10 km.            
        """
        # PRF min is 2538 Hz
        f_Pmin = 1
        f_Pmax = 2538
        v_sc = 7.613
        v_x = 7.0596
        alt_km = 500
        instru_look_angle_rad = np.deg2rad(18.5) 
        tau_p = 30e-6
        D_az = 6
        D_elv = 1.75
        fc = 9.65e9
        self.assertIsNone(self.microxsar.find_valid_highest_possible_PRF(f_Pmin, f_Pmax, v_sc, v_x, alt_km, instru_look_angle_rad, tau_p, D_az, D_elv, fc)[0])

        # PRF max is 9331 Hz
        f_Pmin = 9331
        f_Pmax = 20000
        v_sc = 7.613
        v_x = 7.0596
        alt_km = 500
        instru_look_angle_rad = np.deg2rad(18.5) 
        tau_p = 30e-6
        D_az = 6
        D_elv = 1.75
        fc = 9.65e9
        self.assertIsNone(self.microxsar.find_valid_highest_possible_PRF(f_Pmin, f_Pmax, v_sc, v_x, alt_km, instru_look_angle_rad, tau_p, D_az, D_elv, fc)[0])

        f_Pmin = 2500
        f_Pmax = 3000
        v_sc = 7.613
        v_x = 7.0596
        alt_km = 500
        instru_look_angle_rad = np.deg2rad(18.5)
        tau_p = 30e-6
        D_az = 6
        D_elv = 1.75 
        fc = 9.65e9
        self.assertEqual(self.microxsar.find_valid_highest_possible_PRF(f_Pmin, f_Pmax, v_sc, v_x, alt_km, instru_look_angle_rad, tau_p, D_az, D_elv, fc)[0], 3000)

        f_Pmin = 2500
        f_Pmax = 4000
        v_sc = 7.613
        v_x = 7.0596        
        alt_km = 500
        instru_look_angle_rad = np.deg2rad(18.5)
        tau_p = 30e-6
        D_az = 6
        D_elv = 1.75 
        fc = 9.65e9
        self.assertEqual(self.microxsar.find_valid_highest_possible_PRF(f_Pmin, f_Pmax, v_sc, v_x, alt_km, instru_look_angle_rad, tau_p, D_az, D_elv, fc)[0], 3916)

        f_Pmin = 2500
        f_Pmax = 4500
        v_sc = 7.613
        v_x = 7.0596        
        alt_km = 500
        instru_look_angle_rad = np.deg2rad(18.5)
        tau_p = 30e-6
        D_az = 6
        D_elv = 1.75 
        fc = 9.65e9
        self.assertEqual(self.microxsar.find_valid_highest_possible_PRF(f_Pmin, f_Pmax, v_sc, v_x, alt_km, instru_look_angle_rad, tau_p, D_az, D_elv, fc)[0], 4184)  

    def test_calc_typ_data_metrics_20200915data(self):
        """ Validate with previous run of the synthetic_aperture_radar_model module. 
            The test data is stored in a json file as a list of instruments and other params (input), 
            with the respective resultant metrics (truth data). 
        """
        Re = 6378.137e3 # [m]
        c = 299792458 # m/s
        dir_path = os.path.dirname(os.path.realpath(__file__))

        with open(dir_path+'/20200915_synthetic_aperture_radar_model_refdata.json', 'r') as f:
            data = FileUtilityFunctions.from_json(f)["foos"]
        
        for _d in data:
            test_sar = SyntheticApertureRadarModel.from_json(_d)
            params  = _d["otherSimParams"]
            metrics = _d["resultantMetrics"]
            if("incidenceAngle" in params):
                h = params["altitude"]*1e3
                orb_speed = np.sqrt(3.986004418e14/(Re + h)) # [m/s]
                inc_deg = params["incidenceAngle"]
                obsv_metrics = test_sar.calc_typ_data_metrics(alt_km = h*1e-3, v_sc_kmps = orb_speed*1e-3, v_g_kmps = orb_speed*1e-3*(Re/(Re+h)), incidence_angle_deg = inc_deg, 
                                              instru_look_angle_from_GP_inc_angle = True)
                self.assertAlmostEqual(obsv_metrics["NESZ [dB]"], metrics["NESZ"], delta = 0.01)
            else:
                pass


    # TODO
    def xtest_calc_typ_data_metrics_1(self):
        """ MicroXSAR typical data-metrics. Truth values taken from the reference (below) are approximately equal to those computed by this
            test. There is ambiguity in the actual operating point used in computation of the truth value (eg: exact PRF within
            the given range, pulse width.)s.
           
            Ref: H. Saito et al., “Compact x-band synthetic aperture radar for 100kg class satellite,” IEICE Transactions on
            Communications, vol. E100.B, no. 9, pp. 1653–1660, 2017
        """
        # 
        """ # 2019 Feb 28 13:27:40 is time at which the ECEF and ECI frames approximately align, hence ECEF to ECI rotation is identity. See <https://www.celnav.de/longterm.htm> online calculator of GMST.
        epoch_JDUT1 = 2458543.06088
        SpacecraftOrbitState = {'Time[JDUT1]':epoch_JDUT1, 'x[km]': (Re + h)*1e-3, 'y[km]': 0, 'z[km]': 0, 'vx[km/s]': 0, 'vy[km/s]': orb_speed*1e-3, 'vz[km/s]': 0} # equatorial orbit, altitude 600 km
        TargetCoords = {'Lat [deg]': 0.798, 'Lon [deg]': 0} # g
        """
        epoch_JDUT1 = 2451623.999630
        # lat = 0, lon = 0 corresponds to [6378, 0, 0] km in ECI for observer position, check using Matlab function: eci2lla([6378, 0, 0] ,[2000 3 20 11 59 28.000])
        SpacecraftOrbitState = {'Time[JDUT1]':epoch_JDUT1, 'x[km]': 6378.137 + 600, 'y[km]': 0, 'z[km]': 0, 'vx[km/s]': 0, 'vy[km/s]': 7.6126, 'vz[km/s]': 0} # equatorial orbit, altitude 600 km
        TargetCoords = {'Lat [deg]': 2, 'Lon [deg]': 0} # incidence angle to the target is 28.9855 deg
        obsv_metrics = self.microxsar.calc_typ_data_metrics(SpacecraftOrbitState, TargetCoords)
        self.assertAlmostEqual(obsv_metrics["Ground Pixel Along-Track Resolution [m]"], 2.118, delta = 0.1)
        self.assertAlmostEqual(obsv_metrics["Ground Pixel Cross-Track Resolution [m]"], 4.9492, delta = 0.1)
        self.assertAlmostEqual(obsv_metrics["NESZ [dB]"], -13, delta = 0.5)
        self.assertAlmostEqual(obsv_metrics["Incidence angle [deg]"], 28.9855, delta = 0.1)
        self.assertAlmostEqual(obsv_metrics["(Nominal) Swath-width [km]"], 37, delta = 1)
        self.assertFalse(obsv_metrics["Coverage [T/F]"])
    
    '''
    '''
    def test_calc_typ_data_metrics_2(self):
        """ SeaSatA typical data-metrics. Truth values taken from the reference (below) are approximately equal to those computed by this
            test. There is ambiguity in the actual operating point used in computation of the truth value (eg: exact PRF within
            the given range, pulse widths.
           
            Ref: D. Bickel, B. Brock, and C. Allen, Spaceborne SAR study: LDRD 92 final report, Mar 1993.
        """
        # 

        epoch_JDUT1 = 2451623.999630
        # lat = 0, lon = 0 corresponds to [6378, 0, 0] km in ECI for observer position, check using Matlab function: eci2lla([6378, 0, 0] ,[2000 3 20 11 59 28.000])
        SpacecraftOrbitState = {'Time[JDUT1]':epoch_JDUT1, 'x[km]': 6378.137 + 800, 'y[km]': 0, 'z[km]': 0, 'vx[km/s]': 0, 'vy[km/s]': 7.4518, 'vz[km/s]': 0} # equatorial orbit, altitude 600 km
        TargetCoords = {'Lat [deg]': 1.99, 'Lon [deg]': 0} # incidence angle to the target is 23 deg
        obsv_metrics = self._seasat.calc_typ_data_metrics(SpacecraftOrbitState, TargetCoords)
        print(obsv_metrics)
        self.assertAlmostEqual(obsv_metrics["Ground Pixel Along-Track Resolution [m]"], 25, delta = 0.1)
        self.assertAlmostEqual(obsv_metrics["Ground Pixel Cross-Track Resolution [m]"], 25, delta = 0.1)
        self.assertAlmostEqual(obsv_metrics["NESZ [dB]"], -25, delta = 0.5)
        self.assertAlmostEqual(obsv_metrics["Incidence angle [deg]"], 23, delta = 0.1)
        self.assertAlmostEqual(obsv_metrics["(Nominal) Swath-width [km]"], 37, delta = 1)
        self.assertFalse(obsv_metrics["Coverage [T/F]"])
    '''

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

