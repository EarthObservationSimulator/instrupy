"""Unit tests for instrupy.basic_sensor.basic_sensor_model.
"""

import unittest
import json
import numpy as np
import sys, os
import random

from instrupy.basic_sensor.basic_sensor_model import BasicSensorModel
from instrupy.util import Orientation, ViewGeometry, SphericalGeometry, ReferenceFrame, SyntheticDataConfiguration

RE = 6378.137 # [km] radius of Earth

class TestBasicSensorModel(unittest.TestCase):

    def __init__(self, *args, **kwargs):
        # 
        super(TestBasicSensorModel, self).__init__(*args, **kwargs)

    def test_from_json_basic(self):
        """ Test initialization of basic sensor in the many different ways allowed.
        """
        # Typical case
        o = BasicSensorModel.from_json('{"name": "Atom", "mass":10, "volume":12.45, "dataRate": 40, "bitsPerPixel": 8, "power": 12, \
                                  "orientation": {"referenceFrame": "SC_BODY_FIXED", "convention": "REF_FRAME_ALIGNED"}, \
                                  "fieldOfViewGeometry": {"shape": "CIRCULAR", "diameter":5 }, \
                                  "maneuver":{"@type": "CIRCULAR", "diameter":10}, \
                                  "numberDetectorRows":5, "numberDetectorCols":10, "@id": "bs1", \
                                  "syntheticDataConfig": {"sourceFilePaths":   ["C:/workspace/gfs_forecast_data/gfs.t12z.pgrb2.0p25.f000.nc", \
                                                                                "C:/workspace/gfs_forecast_data/gfs.t12z.pgrb2.0p25.f001.nc", \
                                                                                "C:/workspace/gfs_forecast_data/gfs.t12z.pgrb2.0p25.f002.nc", \
                                                                                "C:/workspace/gfs_forecast_data/gfs.t12z.pgrb2.0p25.f003.nc", \
                                                                                "C:/workspace/gfs_forecast_data/gfs.t12z.pgrb2.0p25.f004.nc"], \
                                                            "geophysicalVar": "TMP_P0_L1_GLL0", "interpolMethod":"SCIPY_LINEAR"}}')

        self.assertEqual(o._type, "Basic Sensor")

        self.assertEqual(o.name, "Atom")
        self.assertIsInstance(o.name, str)

        self.assertIsInstance(o._id, str)
        self.assertEqual(o._id, "bs1")    

        self.assertEqual(o.mass, 10)
        self.assertIsInstance(o.mass, float)

        self.assertEqual(o.volume, 12.45)
        self.assertIsInstance(o.volume, float)

        self.assertEqual(o.power, 12)
        self.assertIsInstance(o.power, float)

        self.assertEqual(o.dataRate, 40)
        self.assertIsInstance(o.dataRate, float)

        self.assertEqual(o.bitsPerPixel, 8)
        self.assertIsInstance(o.bitsPerPixel, int)

        self.assertIsInstance(o.orientation, Orientation)
        self.assertEqual(o.orientation.ref_frame, ReferenceFrame.SC_BODY_FIXED)
        self.assertEqual(o.orientation.euler_angle1, 0)
        self.assertEqual(o.orientation.euler_angle2, 0)
        self.assertEqual(o.orientation.euler_angle3, 0)
        self.assertEqual(o.orientation.euler_seq1, 1)
        self.assertEqual(o.orientation.euler_seq2, 2)
        self.assertEqual(o.orientation.euler_seq3, 3)

        self.assertIsInstance(o.fieldOfViewGeometry, SphericalGeometry)
        
        self.assertIsInstance(o.fieldOfView, ViewGeometry)
        self.assertEqual(o.fieldOfView, ViewGeometry(orien=Orientation.from_dict({"referenceFrame":"SC_BODY_FIXED", "convention": "REF_FRAME_ALIGNED"}), sph_geom=SphericalGeometry.from_dict({"shape":"Circular", "diameter":5})))

        self.assertIsInstance(o.fieldOfRegard[0], ViewGeometry)
        self.assertEqual(o.fieldOfRegard[0], ViewGeometry(orien=Orientation.from_dict({"referenceFrame":"NADIR_POINTING", "convention": "REF_FRAME_ALIGNED"}), sph_geom=SphericalGeometry.from_dict({"shape":"Circular", "diameter":15})))

        self.assertEqual(o.numberDetectorRows, 5)
        self.assertIsInstance(o.numberDetectorRows, int)

        self.assertEqual(o.numberDetectorCols, 10)
        self.assertIsInstance(o.numberDetectorCols, int)

        self.assertIsInstance(o.syntheticDataConfig, SyntheticDataConfiguration)
        self.assertEqual(o.syntheticDataConfig, SyntheticDataConfiguration.from_dict({"sourceFilePaths":   ["C:/workspace/gfs_forecast_data/gfs.t12z.pgrb2.0p25.f000.nc", \
                                                                                                            "C:/workspace/gfs_forecast_data/gfs.t12z.pgrb2.0p25.f001.nc", \
                                                                                                            "C:/workspace/gfs_forecast_data/gfs.t12z.pgrb2.0p25.f002.nc", \
                                                                                                            "C:/workspace/gfs_forecast_data/gfs.t12z.pgrb2.0p25.f003.nc", \
                                                                                                            "C:/workspace/gfs_forecast_data/gfs.t12z.pgrb2.0p25.f004.nc"], \
                                                                                        "geophysicalVar": "TMP_P0_L1_GLL0", "interpolMethod":"SCIPY_LINEAR"}))
                                    
        self.assertIsInstance(o, BasicSensorModel)
'''
        # Test: Test default initialization of orientation, FOV and FOR fields.        
        o = BasicSensorModel.from_json('{}')
        self.assertIsNone(o.name)
        self.assertIsNone(o.acronym)
        self.assertIsNone(o.mass)
        self.assertIsNone(o.volume)
        self.assertIsNone(o.power)
        self.assertIsNone(o.dataRate)

        self.assertIsInstance(o.fieldOfView, FieldOfView)
        self.assertEqual(o.fieldOfView._coneAngleVec_deg, [25/2])
        self.assertIsNone(o.fieldOfView._clockAngleVec_deg)
        self.assertEqual(o.fieldOfView._AT_fov_deg, 25)
        self.assertEqual(o.fieldOfView._CT_fov_deg, 25)

        self.assertEqual(o.fieldOfRegard._geometry, SensorGeometry.CONICAL) # missing manuever field => FOR = FOV
        self.assertEqual(o.fieldOfRegard._geometry, SensorGeometry.CONICAL)
        self.assertEqual(o.fieldOfRegard._coneAngleVec_deg, [25/2])
        self.assertIsNone(o.fieldOfRegard._clockAngleVec_deg)
        self.assertEqual(o.fieldOfRegard._AT_fov_deg, 25)
        self.assertEqual(o.fieldOfRegard._CT_fov_deg, 25)

        self.assertIsInstance(o.orientation, Orientation)
        self.assertEqual(o.orientation.euler_angle1, 0)
        self.assertEqual(o.orientation.euler_angle2, 0)
        self.assertEqual(o.orientation.euler_angle3, 0)
        self.assertEqual(o.orientation.euler_seq1, 1)
        self.assertEqual(o.orientation.euler_seq2, 2)
        self.assertEqual(o.orientation.euler_seq3, 3)

        self.assertIsNone(o._id)

        # test acronym is initialized to name when no specified
        o = BasicSensorModel.from_json('{"name": "Atom"}' )
        self.assertEqual(o._type, "Basic Sensor")
        self.assertEqual(o.name, "Atom")
        self.assertIsInstance(o.name, str)
        
        # Test: Incomplete field-of-view specification, test that Exception is raised
        with self.assertRaises(Exception):
            BasicSensorModel.from_json('{"name": "Atom","mass":10,"volume":12.45, "fieldOfView": {"sensorGeometry": "RECTANGULAR", "alongTrackFieldOfView": 10 }}')


    def test_calc_typ_data_metrics_1(self):
        """ Simple test involving satellite above POI at (lat = 0,lon = 0). Date chosen so that ECEF and ECI frames are aligned.
            Sensor specs do not influence the below calcs. They do however shall influence the coverage calcs (which is not covered by this test).
            Velocity vector do not influence the calcs.
        """
        o = BasicSensorModel.from_json('{}')  
        epoch_JDUT1 =  2458543.06088 # 2019 Feb 28 13:27:40 is time at which the ECEF and ECI frames approximately align, hence ECEF to ECI rotation is identity. See <https://www.celnav.de/longterm.htm> online calculator of GMST.
        
        SpacecraftOrbitState = {'Time[JDUT1]':epoch_JDUT1, 'x[km]': 6878.137, 'y[km]': 0, 'z[km]': 0, 'vx[km/s]': None, 'vy[km/s]': None, 'vz[km/s]': None} # altitude 500 km
        TargetCoords = {'Lat [deg]': 0, 'Lon [deg]': 0}
        obsv_metrics = o.calc_typ_data_metrics(SpacecraftOrbitState, TargetCoords)
        self.assertTrue(obsv_metrics["Coverage [T/F]"])
        self.assertAlmostEqual(obsv_metrics["Observation Range [km]"], 500, delta = 1)
        self.assertAlmostEqual(obsv_metrics["Incidence angle [deg]"], 0, delta = 0.1)
        self.assertAlmostEqual(obsv_metrics["Look angle [deg]"], 0, delta = 0.1)
        self.assertAlmostEqual(obsv_metrics["Solar Zenith [deg]"], 20.335, delta = 0.1) # precomputed value at the epoch and (lat=0, lon=0) position

    def test_calc_typ_data_metrics_2(self):
        """ Test with satellite above POI at (lat = 0,lon = 0), and making observation to the East.
            Date chosen so that ECEF and ECI frames are aligned.
            Sensor specs do not influence the below calcs. They do however shall influence the coverage calcs (which is not covered by this test).
            Velocity vector do not influence the calcs.
        """
        o = BasicSensorModel.from_json('{}')
        epoch_JDUT1 =  2458543.06088 # 2019 Feb 28 13:27:40 is time at which the ECEF and ECI frames approximately align, hence ECEF to ECI rotation is identity. See <https://www.celnav.de/longterm.htm> online calculator of GMST.
        
        # test with reference model which is good for small angles
        poi_lon_deg = random.uniform(0.01, 0.1)
        nadir_angle_deg = abs(np.rad2deg(RE*np.deg2rad(poi_lon_deg)/ 500)) # approximate model, good for small nadir angles
        range_km = 500/np.cos(np.deg2rad(nadir_angle_deg))
        SpacecraftOrbitState = {'Time[JDUT1]':epoch_JDUT1, 'x[km]': 6878.137, 'y[km]': 0, 'z[km]': 0, 'vx[km/s]': None, 'vy[km/s]': None, 'vz[km/s]': None} # altitude 500 km
        TargetCoords = {'Lat [deg]': 0, 'Lon [deg]': poi_lon_deg}
        obsv_metrics = o.calc_typ_data_metrics(SpacecraftOrbitState, TargetCoords)
        self.assertTrue(obsv_metrics["Coverage [T/F]"])
        self.assertAlmostEqual(obsv_metrics["Observation Range [km]"], range_km, delta = 1)
        self.assertAlmostEqual(obsv_metrics["Incidence angle [deg]"], nadir_angle_deg, delta = 0.15)
        self.assertAlmostEqual(obsv_metrics["Look angle [deg]"], nadir_angle_deg, delta = 0.15)
        self.assertAlmostEqual(obsv_metrics["Solar Zenith [deg]"], 20.335, delta = 0.1) # precomputed value at the epoch and near the (lat=0, lon=0) position

        # test with alternate reference model which works for larger angles for a special case with range = RE
        poi_lon_deg = random.uniform(10, 45) 
        range_km = RE # fix range to RE, an isosceles triangle forms
        alt = np.sqrt(RE*RE*(2-2*np.cos(np.deg2rad(180-2*poi_lon_deg)))) - RE
        SpacecraftOrbitState = {'Time[JDUT1]':epoch_JDUT1, 'x[km]': 6378.137+alt, 'y[km]': 0, 'z[km]': 0, 'vx[km/s]': None, 'vy[km/s]': None, 'vz[km/s]': None} # altitude 500 km
        TargetCoords = {'Lat [deg]': 0, 'Lon [deg]': poi_lon_deg}
        obsv_metrics = o.calc_typ_data_metrics(SpacecraftOrbitState, TargetCoords)
        self.assertTrue(obsv_metrics["Coverage [T/F]"])
        self.assertAlmostEqual(obsv_metrics["Observation Range [km]"], range_km, delta = 1)
        self.assertAlmostEqual(obsv_metrics["Incidence angle [deg]"], 2*poi_lon_deg, delta = 0.15)
        self.assertAlmostEqual(obsv_metrics["Look angle [deg]"], poi_lon_deg, delta = 0.15)


    def test_calc_typ_data_metrics_3(self):
        """ Test with satellite above POI at (lat = 0,lon = 0), and making observation to the West.
            Date chosen so that ECEF and ECI frames are aligned.
            Sensor specs do not influence the below calcs. They do however shall influence the coverage calcs (which is not covered by this test).
            Velocity vector do not influence the calcs.
        """
        o = BasicSensorModel.from_json('{}')
        epoch_JDUT1 =  2458543.06088 # 2019 Feb 28 13:27:40 is time at which the ECEF and ECI frames approximately align, hence ECEF to ECI rotation is identity. See <https://www.celnav.de/longterm.htm> online calculator of GMST.
        
        # test with reference model which is good for small angles
        poi_lon_deg = -1 * random.uniform(0.01, 0.1) 
        nadir_angle_deg = abs(np.rad2deg(RE*np.deg2rad(poi_lon_deg)/ 500)) # approximate model, good for small nadir angles
        range_km = 500/np.cos(np.deg2rad(nadir_angle_deg))
        SpacecraftOrbitState = {'Time[JDUT1]':epoch_JDUT1, 'x[km]': 6878.137, 'y[km]': 0, 'z[km]': 0, 'vx[km/s]': None, 'vy[km/s]': None, 'vz[km/s]': None} # altitude 500 km
        TargetCoords = {'Lat [deg]': 0, 'Lon [deg]': poi_lon_deg}
        obsv_metrics = o.calc_typ_data_metrics(SpacecraftOrbitState, TargetCoords)
        self.assertTrue(obsv_metrics["Coverage [T/F]"])
        self.assertAlmostEqual(obsv_metrics["Observation Range [km]"], range_km, delta = 1)
        self.assertAlmostEqual(obsv_metrics["Incidence angle [deg]"], nadir_angle_deg, delta = 0.15) 
        self.assertAlmostEqual(obsv_metrics["Look angle [deg]"], nadir_angle_deg, delta = 0.15)
        self.assertAlmostEqual(obsv_metrics["Solar Zenith [deg]"], 20.335, delta = 0.1) # precomputed value at the epoch and near the (lat=0, lon=0) position

        # test with alternate reference model which works for larger angles for a special case with range = RE
        poi_lon_deg = -1 * random.uniform(10, 45)
        range_km = RE # fix range to RE, an isosceles triangle forms
        alt = np.sqrt(RE*RE*(2-2*np.cos(np.deg2rad(180-2*abs(poi_lon_deg))))) - RE
        SpacecraftOrbitState = {'Time[JDUT1]':epoch_JDUT1, 'x[km]': 6378.137+alt, 'y[km]': 0, 'z[km]': 0, 'vx[km/s]': None, 'vy[km/s]': None, 'vz[km/s]': None} # altitude 500 km
        TargetCoords = {'Lat [deg]': 0, 'Lon [deg]': poi_lon_deg}
        obsv_metrics = o.calc_typ_data_metrics(SpacecraftOrbitState, TargetCoords)
        self.assertTrue(obsv_metrics["Coverage [T/F]"])
        self.assertAlmostEqual(obsv_metrics["Observation Range [km]"], range_km, delta = 1)
        self.assertAlmostEqual(obsv_metrics["Incidence angle [deg]"], 2*abs(poi_lon_deg), delta = 0.15)
        self.assertAlmostEqual(obsv_metrics["Look angle [deg]"], abs(poi_lon_deg), delta = 0.15)


    def test_calc_typ_data_metrics_4(self):
        """ Test with satellite above POI at (lat = 0,lon = 0), and making observation to South.
            Date chosen so that ECEF and ECI frames are aligned.
            Sensor specs do not influence the below calcs. They do however shall influence the coverage calcs (which is not covered by this test).
            Velocity vector do not influence the calcs.
        """
        o = BasicSensorModel.from_json('{}')
        epoch_JDUT1 =  2458543.06088 # 2019 Feb 28 13:27:40 is time at which the ECEF and ECI frames approximately align, hence ECEF to ECI rotation is identity. See <https://www.celnav.de/longterm.htm> online calculator of GMST.
        
        # test with reference model which is good for small angles
        poi_lat_deg = -1 * random.uniform(0.01, 0.1)
        nadir_angle_deg = abs(np.rad2deg(RE*np.deg2rad(poi_lat_deg)/ 500))
        range_km = 500/np.cos(np.deg2rad(nadir_angle_deg))
        SpacecraftOrbitState = {'Time[JDUT1]':epoch_JDUT1, 'x[km]': 6878.137, 'y[km]': 0, 'z[km]': 0, 'vx[km/s]': None, 'vy[km/s]': None, 'vz[km/s]': None} # altitude 500 km
        TargetCoords = {'Lat [deg]': poi_lat_deg, 'Lon [deg]': 0}
        obsv_metrics = o.calc_typ_data_metrics(SpacecraftOrbitState, TargetCoords)
        self.assertTrue(obsv_metrics["Coverage [T/F]"])
        self.assertAlmostEqual(obsv_metrics["Observation Range [km]"], range_km, delta = 1)
        self.assertAlmostEqual(obsv_metrics["Incidence angle [deg]"], nadir_angle_deg, delta = 0.15)
        self.assertAlmostEqual(obsv_metrics["Look angle [deg]"], nadir_angle_deg, delta = 0.15)
        self.assertAlmostEqual(obsv_metrics["Solar Zenith [deg]"], 20.335, delta = 0.1) # precomputed value at the epoch and near the (lat=0, lon=0) position

        # test with alternate reference model which works for larger angles for a special case with range = RE
        poi_lat_deg = random.uniform(10, 45) 
        range_km = RE # fix range to RE, an isosceles triangle forms
        alt = np.sqrt(RE*RE*(2-2*np.cos(np.deg2rad(180-2*poi_lat_deg)))) - RE
        SpacecraftOrbitState = {'Time[JDUT1]':epoch_JDUT1, 'x[km]': 6378.137+alt, 'y[km]': 0, 'z[km]': 0, 'vx[km/s]': None, 'vy[km/s]': None, 'vz[km/s]': None} # altitude 500 km
        TargetCoords = {'Lat [deg]': poi_lat_deg, 'Lon [deg]': 0}
        obsv_metrics = o.calc_typ_data_metrics(SpacecraftOrbitState, TargetCoords)
        self.assertTrue(obsv_metrics["Coverage [T/F]"])
        self.assertAlmostEqual(obsv_metrics["Observation Range [km]"], range_km, delta = 1)
        self.assertAlmostEqual(obsv_metrics["Incidence angle [deg]"], 2*poi_lat_deg, delta = 0.15)
        self.assertAlmostEqual(obsv_metrics["Look angle [deg]"], poi_lat_deg, delta = 0.15)
        

    def test_calc_typ_data_metrics_5(self):
        """ Test with satellite above POI at (lat = 0,lon = 0), and making observation to North.
            Date chosen so that ECEF and ECI frames are aligned.
            Sensor specs do not influence the below calcs. They do however shall influence the coverage calcs (which is not covered by this test).
            Velocity vector do not influence the calcs.
        """
        o = BasicSensorModel.from_json('{}')
        epoch_JDUT1 =  2458543.06088 # 2019 Feb 28 13:27:40 is time at which the ECEF and ECI frames approximately align, hence ECEF to ECI rotation is identity. See <https://www.celnav.de/longterm.htm> online calculator of GMST.
        
        # test with reference model which is good for small angles
        poi_lat_deg = random.uniform(0.01, 0.1)
        nadir_angle_deg = abs(np.rad2deg(RE*np.deg2rad(poi_lat_deg)/ 500))
        range_km = 500/np.cos(np.deg2rad(nadir_angle_deg))
        SpacecraftOrbitState = {'Time[JDUT1]':epoch_JDUT1, 'x[km]': 6878.137, 'y[km]': 0, 'z[km]': 0, 'vx[km/s]': None, 'vy[km/s]': None, 'vz[km/s]': None} # altitude 500 km
        TargetCoords = {'Lat [deg]': poi_lat_deg, 'Lon [deg]': 0}
        obsv_metrics = o.calc_typ_data_metrics(SpacecraftOrbitState, TargetCoords)
        self.assertTrue(obsv_metrics["Coverage [T/F]"])
        self.assertAlmostEqual(obsv_metrics["Observation Range [km]"], range_km, delta = 1)
        self.assertAlmostEqual(obsv_metrics["Incidence angle [deg]"], nadir_angle_deg, delta =  0.15)
        self.assertAlmostEqual(obsv_metrics["Look angle [deg]"], nadir_angle_deg, delta =  0.15)
        self.assertAlmostEqual(obsv_metrics["Solar Zenith [deg]"], 20.335, delta = 0.1) # precomputed value at the epoch and near the (lat=0, lon=0) position

        # test with alternate reference model which works for larger angles for a special case with range = RE
        poi_lat_deg = -1 * random.uniform(10, 45)
        range_km = RE # fix range to RE, an isosceles triangle forms
        alt = np.sqrt(RE*RE*(2-2*np.cos(np.deg2rad(180-2*abs(poi_lat_deg))))) - RE
        SpacecraftOrbitState = {'Time[JDUT1]':epoch_JDUT1, 'x[km]': 6378.137+alt, 'y[km]': 0, 'z[km]': 0, 'vx[km/s]': None, 'vy[km/s]': None, 'vz[km/s]': None} # altitude 500 km
        TargetCoords = {'Lat [deg]': poi_lat_deg, 'Lon [deg]': 0}
        obsv_metrics = o.calc_typ_data_metrics(SpacecraftOrbitState, TargetCoords)
        self.assertTrue(obsv_metrics["Coverage [T/F]"])
        self.assertAlmostEqual(obsv_metrics["Observation Range [km]"], range_km, delta = 1)
        self.assertAlmostEqual(obsv_metrics["Incidence angle [deg]"], 2*abs(poi_lat_deg), delta = 0.15)
        self.assertAlmostEqual(obsv_metrics["Look angle [deg]"], abs(poi_lat_deg), delta = 0.15)

'''