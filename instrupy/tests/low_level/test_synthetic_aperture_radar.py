"""Unit tests for instrupy.synthetic_aperture_radar module.
"""

import unittest
import json
import numpy
import sys, os


from instrupy.synthetic_aperture_radar import *
from instrupy.util import Orientation, FieldOfView


class TestSyntheticApertureRadar(unittest.TestCase):

    microxsar = SyntheticApertureRadar.from_json( '{"@type": "Synthetic Aperture Radar",'
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
                                                  '"pulseWidth": 1e-6,'
                                                  '"antennaAlongTrackDim": 4.9,'
                                                  '"antennaCrossTrackDim": 0.7,' 
                                                  '"antennaApertureEfficiency": 0.7,' 
                                                  '"operatingFrequency": 9.65e9,' 
                                                  '"peakTransmitPower": 1e3,' 
                                                  '"chirpBandwidth": 300e6,'      
                                                  '"minimumPRF": 3000,' 
                                                  '"maximumPRF": 8000,' 
                                                  '"radarLosses": 3.5,' 
                                                  '"sceneNoiseTemp": 290,' 
                                                  '"systemNoiseFigure": 4.3,'
                                                  '"sigmaNEZ0threshold": -10}')


    def test_from_json_basic(self):
         # Test: Typical case      
        self.assertEqual(self.microxsar._type, "Synthetic Aperture Radar")
        self.assertEqual(self.microxsar.name, "MiroXSAR")
        self.assertIsInstance(self.microxsar.name, str)
        self.assertEqual(self.microxsar.mass, 130)
        self.assertIsInstance(self.microxsar.mass, float)
        self.assertEqual(self.microxsar.volume, 0.343)
        self.assertIsInstance(self.microxsar.volume, float)
        self.assertEqual(self.microxsar.power, 1100)
        self.assertIsInstance(self.microxsar.power, float)
        self.assertEqual(self.microxsar.dataRate, 2000)
        self.assertIsInstance(self.microxsar.dataRate, float)
        self.assertEqual(self.microxsar.bitsPerPixel, 16)
        self.assertIsInstance(self.microxsar.bitsPerPixel, int)
        self.assertIsInstance(self.microxsar.orientation, Orientation)
        self.assertEqual(self.microxsar.orientation.euler_angle1, 0)
        self.assertEqual(self.microxsar.orientation.euler_angle2, 30)
        self.assertEqual(self.microxsar.orientation.euler_angle3, 0)
        self.assertEqual(self.microxsar.pulseWidth, 1e-6)
        self.assertIsInstance(self.microxsar.pulseWidth, float)
        self.assertEqual(self.microxsar.antennaAlongTrackDim, 4.9)
        self.assertIsInstance(self.microxsar.antennaAlongTrackDim, float)
        self.assertEqual(self.microxsar.antennaCrossTrackDim, 0.7)
        self.assertIsInstance(self.microxsar.antennaCrossTrackDim, float)
        self.assertEqual(self.microxsar.antennaApertureEfficiency, 0.7)
        self.assertIsInstance(self.microxsar.antennaApertureEfficiency, float)
        self.assertEqual(self.microxsar.operatingFrequency, 9.65e9)
        self.assertIsInstance(self.microxsar.operatingFrequency, float)
        self.assertEqual(self.microxsar.peakTransmitPower, 1e3)
        self.assertIsInstance(self.microxsar.peakTransmitPower, float)
        self.assertEqual(self.microxsar.chirpBandwidth, 300e6)
        self.assertIsInstance(self.microxsar.chirpBandwidth, float)
        self.assertEqual(self.microxsar.minimumPRF, 3000)
        self.assertIsInstance(self.microxsar.minimumPRF, float)
        self.assertEqual(self.microxsar.maximumPRF, 8000)
        self.assertIsInstance(self.microxsar.maximumPRF, float)
        self.assertEqual(self.microxsar.radarLosses, 3.5)
        self.assertIsInstance(self.microxsar.radarLosses, float)
        self.assertEqual(self.microxsar.sceneNoiseTemp, 290)
        self.assertIsInstance(self.microxsar.sceneNoiseTemp, float)
        self.assertEqual(self.microxsar.systemNoiseFigure, 4.3)
        self.assertIsInstance(self.microxsar.systemNoiseFigure, float)
        self.assertEqual(self.microxsar.sigmaNEZ0threshold, -10)
        self.assertIsInstance(self.microxsar.sigmaNEZ0threshold, float)
        # along-track fov and cross-track fov are calculated from antenna size and initialized in the object
        self.assertAlmostEqual(self.microxsar.fieldOfView.get_rectangular_fov_specs_from_custom_fov_specs()[0], 0.3635, places = 2)
        self.assertAlmostEqual(self.microxsar.fieldOfView.get_rectangular_fov_specs_from_custom_fov_specs()[1], 2.5446, places = 2)


        
        # Test of an improper orientation specification. Although physically specifing XYZ as (0,10,0) degrees is the same as specifying 
        # the side-look angle as 10 deg, the behavior of the code is to throw an error.
        with self.assertRaises(Exception):
            o = SyntheticApertureRadar.from_json( '{"@type": "Synthetic Aperture Radar",'
                                                  '"name": "MiroXSAR",'  
                                                  '"mass": 130,' 
                                                  '"volume": 0.343,' 
                                                  '"power": 1100,' 
                                                  '"orientation": {'
                                                  '   "convention": "XYZ",'
                                                  '   "xRotation": 0,'
                                                  '   "yRotation": 10,'
                                                  '   "zRotation": 0'
                                                  ' },'
                                                  '"dataRate": 2000,'
                                                  '"bitsPerPixel": 16,'
                                                  '"pulseWidth": 1e-6,'
                                                  '"antennaAlongTrackDim": 4.9,'
                                                  '"antennaCrossTrackDim": 0.7,' 
                                                  '"antennaApertureEfficiency": 0.7,' 
                                                  '"operatingFrequency": 9.65e9,' 
                                                  '"peakTransmitPower": 1e3,' 
                                                  '"chirpBandwidth": 300e6,'      
                                                  '"minimumPRF": 3000,' 
                                                  '"maximumPRF": 8000,' 
                                                  '"radarLosses": 3.5,' 
                                                  '"sceneNoiseTemp": 290,' 
                                                  '"systemNoiseFigure": 4.3,'
                                                  '"sigmaNEZ0threshold": -10}')

        # Check for invalid PRF min, max specification
        with self.assertRaises(Exception):
            o = SyntheticApertureRadar.from_json( '{"@type": "Synthetic Aperture Radar",'
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
                                                  '"pulseWidth": 1e-6,'
                                                  '"antennaAlongTrackDim": 4.9,'
                                                  '"antennaCrossTrackDim": 0.7,' 
                                                  '"antennaApertureEfficiency": 0.7,' 
                                                  '"operatingFrequency": 9.65e9,' 
                                                  '"peakTransmitPower": 1e3,' 
                                                  '"chirpBandwidth": 300e6,'      
                                                  '"minimumPRF": 8000,' 
                                                  '"maximumPRF": 3000,' 
                                                  '"radarLosses": 3.5,' 
                                                  '"sceneNoiseTemp": 290,' 
                                                  '"systemNoiseFigure": 4.3,'
                                                  '"sigmaNEZ0threshold": -10}')


        
    def test_find_valid_highest_possible_PRF(self):
        """ Truth test data from Spaceborne SAR Study: LDRD 92 Final Report SANDIA Report March 1993, Page 35
            Using look angle corresponding to incidence angle of 20 deg.
            Choose Delev and fc to get a swath of 10 km.
            
        """
        # PRF min is 2538 Hz
        f_Pmin = 1
        f_Pmax = 2538
        v_sc = 7.613
        v_x = 7.0596
        alt_km = 500
        instru_look_angle_rad = numpy.deg2rad(18.5) 
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
        instru_look_angle_rad = numpy.deg2rad(18.5) 
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
        instru_look_angle_rad = numpy.deg2rad(18.5)
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
        instru_look_angle_rad = numpy.deg2rad(18.5)
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
        instru_look_angle_rad = numpy.deg2rad(18.5)
        tau_p = 30e-6
        D_az = 6
        D_elv = 1.75 
        fc = 9.65e9
        self.assertEqual(self.microxsar.find_valid_highest_possible_PRF(f_Pmin, f_Pmax, v_sc, v_x, alt_km, instru_look_angle_rad, tau_p, D_az, D_elv, fc)[0], 4184)
        
    

    def test_calc_typ_data_metrics(self):
        epoch_JDUT1 = 2451623.999630
        SpacecraftOrbitState = {'Time[JDUT1]':epoch_JDUT1, 'x[km]': 6878.137, 'y[km]': 0, 'z[km]': 0, 'vx[km/s]': 0, 'vy[km/s]': 7.6126, 'vz[km/s]': 0} # equatorial orbit, altitude about 500 km
        TargetCoords = {'Lat [deg]': 0, 'Lon [deg]': 0} # lat = 0, lon = 0 corresponds to [6378, 0, 0] km in ECI for observer position, check using Matlab function: eci2lla([6378, 0, 0] ,[2000 3 20 11 59 28.000])
        obsv_metrics = self.microxsar.calc_typ_data_metrics(SpacecraftOrbitState, TargetCoords)
        self.assertAlmostEqual(obsv_metrics["Ground Pixel Along-Track Resolution [m]"], 2.118, delta = 0.1)
        


