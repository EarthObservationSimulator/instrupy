"""Unit tests for instrupy.basic_sensor_model.
"""
import unittest
import numpy as np
import random
from deepdiff import DeepDiff

from instrupy.basic_sensor_model import BasicSensorModel
from instrupy.util import Orientation, ViewGeometry, SphericalGeometry, ReferenceFrame, SyntheticDataConfiguration, Maneuver

RE = 6378.137 # [km] radius of Earth

def orbital_speed(alt_km):
    return np.sqrt(398600.5/(RE + alt_km))

class TestBasicSensorModel(unittest.TestCase):

    def test_from_json(self):
        """ Test initialization of basic sensor in the many different ways allowed.
        """
        # Typical case
        o = BasicSensorModel.from_json('{"name": "Atom", "mass":10, "volume":12.45, "dataRate": 40, "bitsPerPixel": 8, "power": 12, \
                                  "orientation": {"referenceFrame": "SC_BODY_FIXED", "convention": "REF_FRAME_ALIGNED"}, \
                                  "fieldOfViewGeometry": {"shape": "CIRCULAR", "diameter":5 }, \
                                  "maneuver":{"maneuverType": "CIRCULAR", "diameter":10}, \
                                  "numberDetectorRows":5, "numberDetectorCols":10, "@id": "bs1", \
                                  "pointingOption": [{"referenceFrame": "NADIR_POINTING", "convention": "REF_FRAME_ALIGNED"}, \
                                                     {"referenceFrame": "NADIR_POINTING", "convention": "SIDE_LOOK","sideLookAngle":10}, \
                                                     {"referenceFrame": "NADIR_POINTING", "convention": "SIDE_LOOK","sideLookAngle":-10} \
                                                    ], \
                                  "syntheticDataConfig": {"sourceFilePaths":   ["C:/workspace/gfs_forecast_data/gfs.t12z.pgrb2.0p25.f000.nc", \
                                                                                "C:/workspace/gfs_forecast_data/gfs.t12z.pgrb2.0p25.f001.nc", \
                                                                                "C:/workspace/gfs_forecast_data/gfs.t12z.pgrb2.0p25.f002.nc", \
                                                                                "C:/workspace/gfs_forecast_data/gfs.t12z.pgrb2.0p25.f003.nc", \
                                                                                "C:/workspace/gfs_forecast_data/gfs.t12z.pgrb2.0p25.f004.nc"], \
                                                            "geophysicalVar": "TMP_P0_L1_GLL0", "interpolMethod":"SCIPY_LINEAR"}}')

        self.assertIsInstance(o, BasicSensorModel)
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
        self.assertEqual(o.orientation, Orientation.from_dict({"referenceFrame": "SC_BODY_FIXED", "convention": "REF_FRAME_ALIGNED"}))
        
        self.assertIsInstance(o.fieldOfView, ViewGeometry)
        self.assertEqual(o.fieldOfView, ViewGeometry(orien=Orientation.from_dict({"referenceFrame":"SC_BODY_FIXED", "convention": "REF_FRAME_ALIGNED"}), sph_geom=SphericalGeometry.from_dict({"shape":"Circular", "diameter":5})))
       
        self.assertIsInstance(o.sceneFieldOfView, ViewGeometry)
        self.assertEqual(o.sceneFieldOfView, o.fieldOfView) # since sceneFieldOfView is not initialized, it has is equal to the instrument fieldOfView

        self.assertIsInstance(o.maneuver, Maneuver)
        self.assertEqual(o.maneuver, Maneuver.from_dict({"maneuverType": "CIRCULAR", "diameter":10}))

        self.assertIsInstance(o.pointingOption, list)
        self.assertIsInstance(o.pointingOption[0], Orientation)
        self.assertIsInstance(o.pointingOption[1], Orientation)
        self.assertIsInstance(o.pointingOption[2], Orientation)
        self.assertEqual(o.pointingOption[0], Orientation.from_dict({"referenceFrame": "NADIR_POINTING", "convention": "REF_FRAME_ALIGNED"}))
        self.assertEqual(o.pointingOption[1], Orientation.from_dict({"referenceFrame": "NADIR_POINTING", "convention": "SIDE_LOOK","sideLookAngle":10}))
        self.assertEqual(o.pointingOption[2], Orientation.from_dict({"referenceFrame": "NADIR_POINTING", "convention": "SIDE_LOOK","sideLookAngle":-10}))

        self.assertIsInstance(o.fieldOfRegard[0], ViewGeometry)
        self.assertEqual(o.fieldOfRegard[0], ViewGeometry(orien=Orientation.from_dict({"referenceFrame":"NADIR_POINTING", "convention": "REF_FRAME_ALIGNED"}), sph_geom=SphericalGeometry.from_dict({"shape":"Circular", "diameter":15})))
        self.assertEqual(len(o.fieldOfRegard), 1)

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
                                    
        # Test: Test default initialization of orientation, fieldOfViewGeometry, sceneFieldOfViewGeometry, maneuver, numberDetectorRows, numberDetectorCols        
        o = BasicSensorModel.from_json('{}')

        self.assertIsInstance(o, BasicSensorModel)
        self.assertEqual(o._type, "Basic Sensor")
        self.assertIsNone(o.name)
        self.assertIsNotNone(o._id) # random id is assigned
        self.assertIsNone(o.mass)
        self.assertIsNone(o.volume)
        self.assertIsNone(o.power)
        self.assertIsNone(o.dataRate)
        self.assertIsNone(o.bitsPerPixel)

        self.assertIsInstance(o.orientation, Orientation)
        self.assertEqual(o.orientation, Orientation.from_dict({"referenceFrame": "SC_BODY_FIXED", "convention": "REF_FRAME_ALIGNED"}))

        self.assertEqual(o.fieldOfView, ViewGeometry(orien=Orientation.from_dict({"referenceFrame":"SC_BODY_FIXED", "convention": "REF_FRAME_ALIGNED"}), sph_geom=SphericalGeometry.from_dict({"shape":"Circular", "diameter": 25})))

        self.assertIsInstance(o.sceneFieldOfView, ViewGeometry)
        self.assertEqual(o.sceneFieldOfView, o.fieldOfView) # since sceneFieldOfView is not initialized, it has is equal to the instrument fieldOfView

        self.assertIsNone(o.maneuver, Maneuver)
        
        self.assertIsNone(o.fieldOfRegard)

        self.assertEqual(o.numberDetectorRows, 4)
        self.assertIsInstance(o.numberDetectorRows, int)

        self.assertEqual(o.numberDetectorCols, 4)
        self.assertIsInstance(o.numberDetectorCols, int)

        self.assertIsNone(o.syntheticDataConfig)

        # Test sceneFieldOfFViewGeometry and corresponding fieldOfRegard initialization (fieldOfRegard must be built considering sceneFieldOfViewGeometry and not fieldOfViewGeometry)
        o = BasicSensorModel.from_json('{"orientation": {"referenceFrame": "SC_BODY_FIXED", "convention": "REF_FRAME_ALIGNED"}, \
                                         "fieldOfViewGeometry": {"shape": "RECTANGULAR", "angleHeight":0.1, "angleWidth":60 }, \
                                         "sceneFieldOfViewGeometry": {"shape": "RECTANGULAR", "angleHeight":5, "angleWidth":60}, \
                                         "maneuver":{"maneuverType": "SINGLE_ROLL_ONLY", "A_rollMin":0, "A_rollMax": 30} \
                                         }')
        self.assertIsInstance(o, BasicSensorModel)
        self.assertEqual(o._type, "Basic Sensor")
        self.assertEqual(o.fieldOfView, ViewGeometry.from_dict({"orientation": {"referenceFrame":"SC_BODY_FIXED", "convention": "REF_FRAME_ALIGNED"}, 
                                                                "sphericalGeometry": {"shape": "RECTANGULAR", "angleHeight":0.1, "angleWidth":60}}))
        self.assertEqual(o.sceneFieldOfView, ViewGeometry.from_dict({"orientation": {"referenceFrame":"SC_BODY_FIXED", "convention": "REF_FRAME_ALIGNED"}, 
                                                                "sphericalGeometry": {"shape": "RECTANGULAR", "angleHeight":5, "angleWidth":60}}))
        
        ddiff = DeepDiff(o.fieldOfRegard, [ViewGeometry.from_dict({"orientation": {"referenceFrame":"NADIR_POINTING", "convention": "SIDE_LOOK", "sideLookAngle":15}, 
                                                                  "sphericalGeometry": {"shape": "RECTANGULAR", "angleHeight":5, "angleWidth":90}})], ignore_numeric_type_changes=True)
        self.assertEqual(ddiff, {}, msg=ddiff)
        
        # Test: Incomplete field-of-view geometry specification, test that Exception is raised
        with self.assertRaises(Exception):
            BasicSensorModel.from_json('{"name": "Atom","mass":10,"volume":12.45, "fieldOfViewGeometry": {"shape": "RECTANGULAR", "angleHeight": 10 }}')


    def test_calc_data_metrics_1(self):
        """ Simple test involving satellite above POI at (lat = 0,lon = 0). Date chosen so that ECEF and ECI frames are aligned.
            Sensor specs do not influence the below calcs. They do however shall influence the coverage calcs (which is not covered by this test).
        """
        o = BasicSensorModel.from_json('{}')  
        epoch_JDUT1 =  2458543.06088 # 2019 Feb 28 13:27:40 is time at which the ECEF and ECI frames approximately align, hence ECEF to ECI rotation is identity. See <https://www.celnav.de/longterm.htm> online calculator of GMST.
        
        SpacecraftOrbitState = {'time [JDUT1]':epoch_JDUT1, 'x [km]': 6878.137, 'y [km]': 0, 'z [km]': 0, 'vx [km/s]': 0, 'vy [km/s]': 7.6126, 'vz [km/s]': 0} # altitude 500 km
        TargetCoords = {'lat [deg]': 0, 'lon [deg]': 0}
        obsv_metrics = o.calc_data_metrics(SpacecraftOrbitState, TargetCoords)
        self.assertAlmostEqual(obsv_metrics["observation range [km]"], 500, delta = 1)
        self.assertAlmostEqual(obsv_metrics["incidence angle [deg]"], 0, delta = 0.1)
        self.assertAlmostEqual(obsv_metrics["look angle [deg]"], 0, delta = 0.1)
        self.assertAlmostEqual(obsv_metrics["solar zenith [deg]"], 20.335, delta = 0.1) # precomputed value at the epoch and (lat=0, lon=0) position

    def test_calc_data_metrics_2_1(self):
        """ Test with satellite above POI at (lat = 0,lon = 0), and making observation to the East.
            Date chosen so that ECEF and ECI frames are aligned.
            Sensor specs do not influence the below calcs. They do however shall influence the coverage calcs (which is not covered by this test).
            Velocity vector may influence the sign of look angle.
            Test with reference model which is good for small angles.
        """
        o = BasicSensorModel.from_json('{}')
        epoch_JDUT1 =  2458543.06088 # 2019 Feb 28 13:27:40 is time at which the ECEF and ECI frames approximately align, hence ECEF to ECI rotation is identity. See <https://www.celnav.de/longterm.htm> online calculator of GMST.
        
        # 0 deg orbit inclination
        poi_lon_deg = random.uniform(0.01, 0.1)
        nadir_angle_deg = abs(np.rad2deg(RE*np.deg2rad(poi_lon_deg)/ 500)) # approximate model, good for small nadir angles
        range_km = 500/np.cos(np.deg2rad(nadir_angle_deg))
        SpacecraftOrbitState = {'time [JDUT1]':epoch_JDUT1, 'x [km]': 6878.137, 'y [km]': 0, 'z [km]': 0, 'vx [km/s]': 0, 'vy [km/s]': 7.6126, 'vz [km/s]': 0} # altitude 500 km
        TargetCoords = {'lat [deg]': 0, 'lon [deg]': poi_lon_deg}
        obsv_metrics = o.calc_data_metrics(SpacecraftOrbitState, TargetCoords)
        self.assertAlmostEqual(obsv_metrics["observation range [km]"], range_km, delta = 1)
        self.assertAlmostEqual(obsv_metrics["incidence angle [deg]"], nadir_angle_deg, delta = 0.15)
        self.assertAlmostEqual(obsv_metrics["look angle [deg]"], nadir_angle_deg, delta = 0.15)
        self.assertAlmostEqual(obsv_metrics["solar zenith [deg]"], 20.335, delta = 0.1) # precomputed value at the epoch and near the (lat=0, lon=0) position

        # 90 deg orbit inclination
        poi_lon_deg = random.uniform(0.01, 0.1)
        nadir_angle_deg = abs(np.rad2deg(RE*np.deg2rad(poi_lon_deg)/ 500)) # approximate model, good for small nadir angles
        range_km = 500/np.cos(np.deg2rad(nadir_angle_deg))
        SpacecraftOrbitState = {'time [JDUT1]':epoch_JDUT1, 'x [km]': 6878.137, 'y [km]': 0, 'z [km]': 0, 'vx [km/s]': 0, 'vy [km/s]': 0, 'vz [km/s]': -7.6126} # altitude 500 km
        TargetCoords = {'lat [deg]': 0, 'lon [deg]': poi_lon_deg}
        obsv_metrics = o.calc_data_metrics(SpacecraftOrbitState, TargetCoords)
        self.assertAlmostEqual(obsv_metrics["observation range [km]"], range_km, delta = 1)
        self.assertAlmostEqual(obsv_metrics["incidence angle [deg]"], nadir_angle_deg, delta = 0.15)
        self.assertAlmostEqual(obsv_metrics["look angle [deg]"], nadir_angle_deg, delta = 0.15)
        self.assertAlmostEqual(obsv_metrics["solar zenith [deg]"], 20.335, delta = 0.1) # precomputed value at the epoch and near the (lat=0, lon=0) position

        # 90 deg orbit inclination
        poi_lon_deg = random.uniform(0.01, 0.1)
        nadir_angle_deg = abs(np.rad2deg(RE*np.deg2rad(poi_lon_deg)/ 500)) # approximate model, good for small nadir angles
        range_km = 500/np.cos(np.deg2rad(nadir_angle_deg))
        SpacecraftOrbitState = {'time [JDUT1]':epoch_JDUT1, 'x [km]': 6878.137, 'y [km]': 0, 'z [km]': 0, 'vx [km/s]': 0, 'vy [km/s]': 0, 'vz [km/s]': 7.6126} # altitude 500 km
        TargetCoords = {'lat [deg]': 0, 'lon [deg]': poi_lon_deg}
        obsv_metrics = o.calc_data_metrics(SpacecraftOrbitState, TargetCoords)
        self.assertAlmostEqual(obsv_metrics["observation range [km]"], range_km, delta = 1)
        self.assertAlmostEqual(obsv_metrics["incidence angle [deg]"], nadir_angle_deg, delta = 0.15)
        self.assertAlmostEqual(obsv_metrics["look angle [deg]"], -1*nadir_angle_deg, delta = 0.15)
        self.assertAlmostEqual(obsv_metrics["solar zenith [deg]"], 20.335, delta = 0.1) # precomputed value at the epoch and near the (lat=0, lon=0) position

    def test_calc_data_metrics_2_2(self):
        """ Test with satellite above POI at (lat = 0,lon = 0), and making observation to the East.
            Date chosen so that ECEF and ECI frames are aligned.
            Sensor specs do not influence the below calcs. They do however shall influence the coverage calcs (which is not covered by this test).
            Velocity vector may influence the sign of look angle.
            Test with reference model which works for larger angles for a special case with range = RE
        """
        o = BasicSensorModel.from_json('{}')
        epoch_JDUT1 =  2458543.06088 # 2019 Feb 28 13:27:40 is time at which the ECEF and ECI frames approximately align, hence ECEF to ECI rotation is identity. See <https://www.celnav.de/longterm.htm> online calculator of GMST.
        
        # test 0 deg inclination
        poi_lon_deg = random.uniform(10, 45) 
        range_km = RE # fix range to RE, an isosceles triangle forms
        alt = np.sqrt(RE*RE*(2-2*np.cos(np.deg2rad(180-2*poi_lon_deg)))) - RE
        SpacecraftOrbitState = {'time [JDUT1]':epoch_JDUT1, 'x [km]': 6378.137+alt, 'y [km]': 0, 'z [km]': 0, 'vx [km/s]': 0, 'vy [km/s]': orbital_speed(alt*1e-3), 'vz [km/s]': 0}
        TargetCoords = {'lat [deg]': 0, 'lon [deg]': poi_lon_deg}
        obsv_metrics = o.calc_data_metrics(SpacecraftOrbitState, TargetCoords)
        self.assertAlmostEqual(obsv_metrics["observation range [km]"], range_km, delta = 1)
        self.assertAlmostEqual(obsv_metrics["incidence angle [deg]"], 2*poi_lon_deg, delta = 0.15)
        self.assertAlmostEqual(obsv_metrics["look angle [deg]"], poi_lon_deg, delta = 0.15)

        # test positive look angle, 90 deg orbit inclination
        poi_lon_deg = random.uniform(10, 45) 
        range_km = RE # fix range to RE, an isosceles triangle forms
        alt = np.sqrt(RE*RE*(2-2*np.cos(np.deg2rad(180-2*poi_lon_deg)))) - RE
        SpacecraftOrbitState = {'time [JDUT1]':epoch_JDUT1, 'x [km]': 6378.137+alt, 'y [km]': 0, 'z [km]': 0, 'vx [km/s]': 0, 'vy [km/s]': 0, 'vz [km/s]': -1*orbital_speed(alt*1e-3)}
        TargetCoords = {'lat [deg]': 0, 'lon [deg]': poi_lon_deg}
        obsv_metrics = o.calc_data_metrics(SpacecraftOrbitState, TargetCoords)
        self.assertAlmostEqual(obsv_metrics["observation range [km]"], range_km, delta = 1)
        self.assertAlmostEqual(obsv_metrics["incidence angle [deg]"], 2*poi_lon_deg, delta = 0.15)
        self.assertAlmostEqual(obsv_metrics["look angle [deg]"], poi_lon_deg, delta = 0.15)
        # test negative look angle, 90 deg orbit inclination
        poi_lon_deg = random.uniform(10, 45) 
        range_km = RE # fix range to RE, an isosceles triangle forms
        alt = np.sqrt(RE*RE*(2-2*np.cos(np.deg2rad(180-2*poi_lon_deg)))) - RE
        SpacecraftOrbitState = {'time [JDUT1]':epoch_JDUT1, 'x [km]': 6378.137+alt, 'y [km]': 0, 'z [km]': 0, 'vx [km/s]': 0, 'vy [km/s]': 0, 'vz [km/s]':  orbital_speed(alt*1e-3)}
        TargetCoords = {'lat [deg]': 0, 'lon [deg]': poi_lon_deg}
        obsv_metrics = o.calc_data_metrics(SpacecraftOrbitState, TargetCoords)
        self.assertAlmostEqual(obsv_metrics["observation range [km]"], range_km, delta = 1)
        self.assertAlmostEqual(obsv_metrics["incidence angle [deg]"], 2*poi_lon_deg, delta = 0.15)
        self.assertAlmostEqual(obsv_metrics["look angle [deg]"], -1*poi_lon_deg, delta = 0.15)

    def test_calc_data_metrics_3_1(self):
        """ Test with satellite above POI at (lat = 0,lon = 0), and making observation to the West.
            Date chosen so that ECEF and ECI frames are aligned.
            Sensor specs do not influence the below calcs. They do however shall influence the coverage calcs (which is not covered by this test).
            Velocity vector may influence the sign of look angle.
            Test with reference model which is good for small angles.
        """
        o = BasicSensorModel.from_json('{}')
        epoch_JDUT1 =  2458543.06088 # 2019 Feb 28 13:27:40 is time at which the ECEF and ECI frames approximately align, hence ECEF to ECI rotation is identity. See <https://www.celnav.de/longterm.htm> online calculator of GMST.

        # test 0 deg inclination     
        poi_lon_deg = -1 * random.uniform(0.01, 0.1) 
        nadir_angle_deg = abs(np.rad2deg(RE*np.deg2rad(poi_lon_deg)/ 500)) # approximate model, good for small nadir angles
        range_km = 500/np.cos(np.deg2rad(nadir_angle_deg))
        SpacecraftOrbitState = {'time [JDUT1]':epoch_JDUT1, 'x [km]': 6878.137, 'y [km]': 0, 'z [km]': 0, 'vx [km/s]': 0, 'vy [km/s]': 7.6126, 'vz [km/s]': 0} # altitude 500 km
        TargetCoords = {'lat [deg]': 0, 'lon [deg]': poi_lon_deg}
        obsv_metrics = o.calc_data_metrics(SpacecraftOrbitState, TargetCoords)
        self.assertAlmostEqual(obsv_metrics["observation range [km]"], range_km, delta = 1)
        self.assertAlmostEqual(obsv_metrics["incidence angle [deg]"], nadir_angle_deg, delta = 0.15) 
        self.assertAlmostEqual(obsv_metrics["look angle [deg]"], nadir_angle_deg, delta = 0.15)
        self.assertAlmostEqual(obsv_metrics["solar zenith [deg]"], 20.335, delta = 0.1) # precomputed value at the epoch and near the (lat=0, lon=0) position

        # test positive look angle, 90 deg orbit inclination     
        poi_lon_deg = -1 * random.uniform(0.01, 0.1) 
        nadir_angle_deg = abs(np.rad2deg(RE*np.deg2rad(poi_lon_deg)/ 500)) # approximate model, good for small nadir angles
        range_km = 500/np.cos(np.deg2rad(nadir_angle_deg))
        SpacecraftOrbitState = {'time [JDUT1]':epoch_JDUT1, 'x [km]': 6878.137, 'y [km]': 0, 'z [km]': 0, 'vx [km/s]': 0, 'vy [km/s]': 0, 'vz [km/s]': 7.6126} # altitude 500 km
        TargetCoords = {'lat [deg]': 0, 'lon [deg]': poi_lon_deg}
        obsv_metrics = o.calc_data_metrics(SpacecraftOrbitState, TargetCoords)
        self.assertAlmostEqual(obsv_metrics["observation range [km]"], range_km, delta = 1)
        self.assertAlmostEqual(obsv_metrics["incidence angle [deg]"], nadir_angle_deg, delta = 0.15) 
        self.assertAlmostEqual(obsv_metrics["look angle [deg]"], nadir_angle_deg, delta = 0.15)
        self.assertAlmostEqual(obsv_metrics["solar zenith [deg]"], 20.335, delta = 0.1) # precomputed value at the epoch and near the (lat=0, lon=0) position

        # test negative look angle, 90 deg orbit inclination     
        poi_lon_deg = -1 * random.uniform(0.01, 0.1) 
        nadir_angle_deg = abs(np.rad2deg(RE*np.deg2rad(poi_lon_deg)/ 500)) # approximate model, good for small nadir angles
        range_km = 500/np.cos(np.deg2rad(nadir_angle_deg))
        SpacecraftOrbitState = {'time [JDUT1]':epoch_JDUT1, 'x [km]': 6878.137, 'y [km]': 0, 'z [km]': 0, 'vx [km/s]': 0, 'vy [km/s]': 0, 'vz [km/s]': -7.6126} # altitude 500 km
        TargetCoords = {'lat [deg]': 0, 'lon [deg]': poi_lon_deg}
        obsv_metrics = o.calc_data_metrics(SpacecraftOrbitState, TargetCoords)
        self.assertAlmostEqual(obsv_metrics["observation range [km]"], range_km, delta = 1)
        self.assertAlmostEqual(obsv_metrics["incidence angle [deg]"], nadir_angle_deg, delta = 0.15) 
        self.assertAlmostEqual(obsv_metrics["look angle [deg]"], -1*nadir_angle_deg, delta = 0.15)
        self.assertAlmostEqual(obsv_metrics["solar zenith [deg]"], 20.335, delta = 0.1) # precomputed value at the epoch and near the (lat=0, lon=0) position

    def test_calc_data_metrics_3_2(self):
        """ Test with satellite above POI at (lat = 0,lon = 0), and making observation to the West.
            Date chosen so that ECEF and ECI frames are aligned.
            Sensor specs do not influence the below calcs. They do however shall influence the coverage calcs (which is not covered by this test).
            Velocity vector may influence the sign of look angle.
            Test with alternate reference model which works for larger angles for a special case with range = RE
        """
        o = BasicSensorModel.from_json('{}')
        epoch_JDUT1 =  2458543.06088 # 2019 Feb 28 13:27:40 is time at which the ECEF and ECI frames approximately align, hence ECEF to ECI rotation is identity. See <https://www.celnav.de/longterm.htm> online calculator of GMST.
            
        # test 0 deg inclination  
        poi_lon_deg = -1 * random.uniform(10, 45)
        range_km = RE # fix range to RE, an isosceles triangle forms
        alt = np.sqrt(RE*RE*(2-2*np.cos(np.deg2rad(180-2*abs(poi_lon_deg))))) - RE
        SpacecraftOrbitState = {'time [JDUT1]':epoch_JDUT1, 'x [km]': 6378.137+alt, 'y [km]': 0, 'z [km]': 0, 'vx [km/s]': 0, 'vy [km/s]': orbital_speed(alt*1e-3), 'vz [km/s]': 0} 
        TargetCoords = {'lat [deg]': 0, 'lon [deg]': poi_lon_deg}
        obsv_metrics = o.calc_data_metrics(SpacecraftOrbitState, TargetCoords)
        self.assertAlmostEqual(obsv_metrics["observation range [km]"], range_km, delta = 1)
        self.assertAlmostEqual(obsv_metrics["incidence angle [deg]"], 2*abs(poi_lon_deg), delta = 0.15)
        self.assertAlmostEqual(obsv_metrics["look angle [deg]"], abs(poi_lon_deg), delta = 0.15)

        # test positive look angle, 90 deg orbit inclination  
        poi_lon_deg = -1 * random.uniform(10, 45)
        range_km = RE # fix range to RE, an isosceles triangle forms
        alt = np.sqrt(RE*RE*(2-2*np.cos(np.deg2rad(180-2*abs(poi_lon_deg))))) - RE
        SpacecraftOrbitState = {'time [JDUT1]':epoch_JDUT1, 'x [km]': 6378.137+alt, 'y [km]': 0, 'z [km]': 0, 'vx [km/s]': 0, 'vy [km/s]': 0, 'vz [km/s]': orbital_speed(alt*1e-3)} 
        TargetCoords = {'lat [deg]': 0, 'lon [deg]': poi_lon_deg}
        obsv_metrics = o.calc_data_metrics(SpacecraftOrbitState, TargetCoords)
        self.assertAlmostEqual(obsv_metrics["observation range [km]"], range_km, delta = 1)
        self.assertAlmostEqual(obsv_metrics["incidence angle [deg]"], 2*abs(poi_lon_deg), delta = 0.15)
        self.assertAlmostEqual(obsv_metrics["look angle [deg]"], abs(poi_lon_deg), delta = 0.15)

        # test negative look angle, 90 deg orbit inclination  
        poi_lon_deg = -1 * random.uniform(10, 45)
        range_km = RE # fix range to RE, an isosceles triangle forms
        alt = np.sqrt(RE*RE*(2-2*np.cos(np.deg2rad(180-2*abs(poi_lon_deg))))) - RE
        SpacecraftOrbitState = {'time [JDUT1]':epoch_JDUT1, 'x [km]': 6378.137+alt, 'y [km]': 0, 'z [km]': 0, 'vx [km/s]': 0, 'vy [km/s]': 0, 'vz [km/s]': -1*orbital_speed(alt*1e-3)} 
        TargetCoords = {'lat [deg]': 0, 'lon [deg]': poi_lon_deg}
        obsv_metrics = o.calc_data_metrics(SpacecraftOrbitState, TargetCoords)
        self.assertAlmostEqual(obsv_metrics["observation range [km]"], range_km, delta = 1)
        self.assertAlmostEqual(obsv_metrics["incidence angle [deg]"], 2*abs(poi_lon_deg), delta = 0.15)
        self.assertAlmostEqual(obsv_metrics["look angle [deg]"], -1*abs(poi_lon_deg), delta = 0.15)

    def test_calc_data_metrics_4_1(self):
        """ Test with satellite above POI at (lat = 0,lon = 0), and making observation to South.
            Date chosen so that ECEF and ECI frames are aligned.
            Sensor specs do not influence the below calcs. They do however shall influence the coverage calcs (which is not covered by this test).
            Velocity vector do not influence the calcs.
            Test with reference model which is good for small angles.
        """
        o = BasicSensorModel.from_json('{}')
        epoch_JDUT1 =  2458543.06088 # 2019 Feb 28 13:27:40 is time at which the ECEF and ECI frames approximately align, hence ECEF to ECI rotation is identity. See <https://www.celnav.de/longterm.htm> online calculator of GMST.
        
        # test positive look angle, 0 deg orbit inclination 
        poi_lat_deg = -1 * random.uniform(0.01, 0.1)
        nadir_angle_deg = abs(np.rad2deg(RE*np.deg2rad(poi_lat_deg)/ 500))
        range_km = 500/np.cos(np.deg2rad(nadir_angle_deg))
        SpacecraftOrbitState = {'time [JDUT1]':epoch_JDUT1, 'x [km]': 6878.137, 'y [km]': 0, 'z [km]': 0, 'vx [km/s]': 0, 'vy [km/s]': -7.6126, 'vz [km/s]': 0} # altitude 500 km
        TargetCoords = {'lat [deg]': poi_lat_deg, 'lon [deg]': 0}
        obsv_metrics = o.calc_data_metrics(SpacecraftOrbitState, TargetCoords)
        self.assertAlmostEqual(obsv_metrics["observation range [km]"], range_km, delta = 1)
        self.assertAlmostEqual(obsv_metrics["incidence angle [deg]"], nadir_angle_deg, delta = 0.15)
        self.assertAlmostEqual(obsv_metrics["look angle [deg]"], nadir_angle_deg, delta = 0.15)
        self.assertAlmostEqual(obsv_metrics["solar zenith [deg]"], 20.335, delta = 0.1) # precomputed value at the epoch and near the (lat=0, lon=0) position

        # test negative look angle, 0 deg orbit inclination 
        poi_lat_deg = -1 * random.uniform(0.01, 0.1)
        nadir_angle_deg = abs(np.rad2deg(RE*np.deg2rad(poi_lat_deg)/ 500))
        range_km = 500/np.cos(np.deg2rad(nadir_angle_deg))
        SpacecraftOrbitState = {'time [JDUT1]':epoch_JDUT1, 'x [km]': 6878.137, 'y [km]': 0, 'z [km]': 0, 'vx [km/s]': 0, 'vy [km/s]': 7.6126, 'vz [km/s]': 0} # altitude 500 km
        TargetCoords = {'lat [deg]': poi_lat_deg, 'lon [deg]': 0}
        obsv_metrics = o.calc_data_metrics(SpacecraftOrbitState, TargetCoords)
        self.assertAlmostEqual(obsv_metrics["observation range [km]"], range_km, delta = 1)
        self.assertAlmostEqual(obsv_metrics["incidence angle [deg]"], nadir_angle_deg, delta = 0.15)
        self.assertAlmostEqual(obsv_metrics["look angle [deg]"], -1*nadir_angle_deg, delta = 0.15)
        self.assertAlmostEqual(obsv_metrics["solar zenith [deg]"], 20.335, delta = 0.1) # precomputed value at the epoch and near the (lat=0, lon=0) position

        # test 90 deg orbit inclination 
        poi_lat_deg = -1 * random.uniform(0.01, 0.1)
        nadir_angle_deg = abs(np.rad2deg(RE*np.deg2rad(poi_lat_deg)/ 500))
        range_km = 500/np.cos(np.deg2rad(nadir_angle_deg))
        SpacecraftOrbitState = {'time [JDUT1]':epoch_JDUT1, 'x [km]': 6878.137, 'y [km]': 0, 'z [km]': 0, 'vx [km/s]': 0, 'vy [km/s]': 0, 'vz [km/s]': 7.6126} # altitude 500 km
        TargetCoords = {'lat [deg]': poi_lat_deg, 'lon [deg]': 0}
        obsv_metrics = o.calc_data_metrics(SpacecraftOrbitState, TargetCoords)
        self.assertAlmostEqual(obsv_metrics["observation range [km]"], range_km, delta = 1)
        self.assertAlmostEqual(obsv_metrics["incidence angle [deg]"], nadir_angle_deg, delta = 0.15)
        self.assertAlmostEqual(obsv_metrics["look angle [deg]"], nadir_angle_deg, delta = 0.15)
        self.assertAlmostEqual(obsv_metrics["solar zenith [deg]"], 20.335, delta = 0.1) # precomputed value at the epoch and near the (lat=0, lon=0) position

    def test_calc_data_metrics_4_2(self):
        """ Test with satellite above POI at (lat = 0,lon = 0), and making observation to South.
            Date chosen so that ECEF and ECI frames are aligned.
            Sensor specs do not influence the below calcs. They do however shall influence the coverage calcs (which is not covered by this test).
            Velocity vector do not influence the calcs.
            Test with alternate reference model which works for larger angles for a special case with range = RE.
        """
        o = BasicSensorModel.from_json('{}')
        epoch_JDUT1 =  2458543.06088 # 2019 Feb 28 13:27:40 is time at which the ECEF and ECI frames approximately align, hence ECEF to ECI rotation is identity. See <https://www.celnav.de/longterm.htm> online calculator of GMST.
        
        # test positive look angle, 0 deg orbit inclination 
        poi_lat_deg = -1*random.uniform(10, 45) 
        range_km = RE # fix range to RE, an isosceles triangle forms
        alt = np.sqrt(RE*RE*(2-2*np.cos(np.deg2rad(180-2*poi_lat_deg)))) - RE
        SpacecraftOrbitState = {'time [JDUT1]':epoch_JDUT1, 'x [km]': 6378.137+alt, 'y [km]': 0, 'z [km]': 0, 'vx [km/s]': 0, 'vy [km/s]': -1*orbital_speed(alt*1e-3), 'vz [km/s]': 0}
        TargetCoords = {'lat [deg]': poi_lat_deg, 'lon [deg]': 0}
        obsv_metrics = o.calc_data_metrics(SpacecraftOrbitState, TargetCoords)
        self.assertAlmostEqual(obsv_metrics["observation range [km]"], range_km, delta = 1)
        self.assertAlmostEqual(obsv_metrics["incidence angle [deg]"], 2*abs(poi_lat_deg), delta = 0.15)
        self.assertAlmostEqual(obsv_metrics["look angle [deg]"], abs(poi_lat_deg), delta = 0.15)

        # test negative look angle, 0 deg orbit inclination 
        poi_lat_deg = -1*random.uniform(10, 45) 
        range_km = RE # fix range to RE, an isosceles triangle forms
        alt = np.sqrt(RE*RE*(2-2*np.cos(np.deg2rad(180-2*poi_lat_deg)))) - RE
        SpacecraftOrbitState = {'time [JDUT1]':epoch_JDUT1, 'x [km]': 6378.137+alt, 'y [km]': 0, 'z [km]': 0, 'vx [km/s]': 0, 'vy [km/s]': orbital_speed(alt*1e-3), 'vz [km/s]': 0} 
        TargetCoords = {'lat [deg]': poi_lat_deg, 'lon [deg]': 0}
        obsv_metrics = o.calc_data_metrics(SpacecraftOrbitState, TargetCoords)
        self.assertAlmostEqual(obsv_metrics["observation range [km]"], range_km, delta = 1)
        self.assertAlmostEqual(obsv_metrics["incidence angle [deg]"], 2*abs(poi_lat_deg), delta = 0.15)
        self.assertAlmostEqual(obsv_metrics["look angle [deg]"], -1*abs(poi_lat_deg), delta = 0.15)

        # test 90 deg orbit inclination 
        poi_lat_deg = -1*random.uniform(10, 45) 
        range_km = RE # fix range to RE, an isosceles triangle forms
        alt = np.sqrt(RE*RE*(2-2*np.cos(np.deg2rad(180-2*poi_lat_deg)))) - RE
        SpacecraftOrbitState = {'time [JDUT1]':epoch_JDUT1, 'x [km]': 6378.137+alt, 'y [km]': 0, 'z [km]': 0, 'vx [km/s]': 0, 'vy [km/s]': 0, 'vz [km/s]': orbital_speed(alt*1e-3)} 
        TargetCoords = {'lat [deg]': poi_lat_deg, 'lon [deg]': 0}
        obsv_metrics = o.calc_data_metrics(SpacecraftOrbitState, TargetCoords)
        self.assertAlmostEqual(obsv_metrics["observation range [km]"], range_km, delta = 1)
        self.assertAlmostEqual(obsv_metrics["incidence angle [deg]"], 2*abs(poi_lat_deg), delta = 0.15)
        self.assertAlmostEqual(obsv_metrics["look angle [deg]"], abs(poi_lat_deg), delta = 0.15)
        
    def test_calc_data_metrics_5_1(self):
        """ Test with satellite above POI at (lat = 0,lon = 0), and making observation to North.
            Date chosen so that ECEF and ECI frames are aligned.
            Sensor specs do not influence the below calcs. They do however shall influence the coverage calcs (which is not covered by this test).
            Velocity vector do not influence the calcs.
            Test with reference model which is good for small angles
        """
        o = BasicSensorModel.from_json('{}')
        epoch_JDUT1 =  2458543.06088 # 2019 Feb 28 13:27:40 is time at which the ECEF and ECI frames approximately align, hence ECEF to ECI rotation is identity. See <https://www.celnav.de/longterm.htm> online calculator of GMST.
        
        # test positive look angle, 0 deg orbit inclination 
        poi_lat_deg = random.uniform(0.01, 0.1)
        nadir_angle_deg = abs(np.rad2deg(RE*np.deg2rad(poi_lat_deg)/ 500))
        range_km = 500/np.cos(np.deg2rad(nadir_angle_deg))
        SpacecraftOrbitState = {'time [JDUT1]':epoch_JDUT1, 'x [km]': 6878.137, 'y [km]': 0, 'z [km]': 0, 'vx [km/s]': 0, 'vy [km/s]': 7.6126, 'vz [km/s]': 0} # altitude 500 km
        TargetCoords = {'lat [deg]': poi_lat_deg, 'lon [deg]': 0}
        obsv_metrics = o.calc_data_metrics(SpacecraftOrbitState, TargetCoords)
        self.assertAlmostEqual(obsv_metrics["observation range [km]"], range_km, delta = 1)
        self.assertAlmostEqual(obsv_metrics["incidence angle [deg]"], nadir_angle_deg, delta =  0.15)
        self.assertAlmostEqual(obsv_metrics["look angle [deg]"], nadir_angle_deg, delta =  0.15)
        self.assertAlmostEqual(obsv_metrics["solar zenith [deg]"], 20.335, delta = 0.1) # precomputed value at the epoch and near the (lat=0, lon=0) position

        # test negative look angle, 0 deg orbit inclination 
        poi_lat_deg = random.uniform(0.01, 0.1)
        nadir_angle_deg = abs(np.rad2deg(RE*np.deg2rad(poi_lat_deg)/ 500))
        range_km = 500/np.cos(np.deg2rad(nadir_angle_deg))
        SpacecraftOrbitState = {'time [JDUT1]':epoch_JDUT1, 'x [km]': 6878.137, 'y [km]': 0, 'z [km]': 0, 'vx [km/s]': 0, 'vy [km/s]': -7.6126, 'vz [km/s]': 0} # altitude 500 km
        TargetCoords = {'lat [deg]': poi_lat_deg, 'lon [deg]': 0}
        obsv_metrics = o.calc_data_metrics(SpacecraftOrbitState, TargetCoords)
        self.assertAlmostEqual(obsv_metrics["observation range [km]"], range_km, delta = 1)
        self.assertAlmostEqual(obsv_metrics["incidence angle [deg]"], nadir_angle_deg, delta =  0.15)
        self.assertAlmostEqual(obsv_metrics["look angle [deg]"], -1*nadir_angle_deg, delta =  0.15)
        self.assertAlmostEqual(obsv_metrics["solar zenith [deg]"], 20.335, delta = 0.1) # precomputed value at the epoch and near the (lat=0, lon=0) position

        # test 90 deg orbit inclination
        poi_lat_deg = random.uniform(0.01, 0.1)
        nadir_angle_deg = abs(np.rad2deg(RE*np.deg2rad(poi_lat_deg)/ 500))
        range_km = 500/np.cos(np.deg2rad(nadir_angle_deg))
        SpacecraftOrbitState = {'time [JDUT1]':epoch_JDUT1, 'x [km]': 6878.137, 'y [km]': 0, 'z [km]': 0, 'vx [km/s]': 0, 'vy [km/s]': 0, 'vz [km/s]': 7.6126} # altitude 500 km
        TargetCoords = {'lat [deg]': poi_lat_deg, 'lon [deg]': 0}
        obsv_metrics = o.calc_data_metrics(SpacecraftOrbitState, TargetCoords)
        self.assertAlmostEqual(obsv_metrics["observation range [km]"], range_km, delta = 1)
        self.assertAlmostEqual(obsv_metrics["incidence angle [deg]"], nadir_angle_deg, delta =  0.15)
        self.assertAlmostEqual(obsv_metrics["look angle [deg]"], nadir_angle_deg, delta =  0.15)
        self.assertAlmostEqual(obsv_metrics["solar zenith [deg]"], 20.335, delta = 0.1) # precomputed value at the epoch and near the (lat=0, lon=0) position

    def test_calc_data_metrics_5_2(self):
        """ Test with satellite above POI at (lat = 0,lon = 0), and making observation to North.
            Date chosen so that ECEF and ECI frames are aligned.
            Sensor specs do not influence the below calcs. They do however shall influence the coverage calcs (which is not covered by this test).
            Velocity vector do not influence the calcs.
            Test with alternate reference model which works for larger angles for a special case with range = RE
        """
        o = BasicSensorModel.from_json('{}')
        epoch_JDUT1 =  2458543.06088 # 2019 Feb 28 13:27:40 is time at which the ECEF and ECI frames approximately align, hence ECEF to ECI rotation is identity. See <https://www.celnav.de/longterm.htm> online calculator of GMST.
        
        # test positive look angle, 0 deg orbit inclination 
        poi_lat_deg = random.uniform(10, 45)
        range_km = RE # fix range to RE, an isosceles triangle forms
        alt = np.sqrt(RE*RE*(2-2*np.cos(np.deg2rad(180-2*abs(poi_lat_deg))))) - RE
        SpacecraftOrbitState = {'time [JDUT1]':epoch_JDUT1, 'x [km]': 6378.137+alt, 'y [km]': 0, 'z [km]': 0, 'vx [km/s]': 0, 'vy [km/s]': orbital_speed(alt*1e-3), 'vz [km/s]': 0} # altitude 500 km
        TargetCoords = {'lat [deg]': poi_lat_deg, 'lon [deg]': 0}
        obsv_metrics = o.calc_data_metrics(SpacecraftOrbitState, TargetCoords)
        self.assertAlmostEqual(obsv_metrics["observation range [km]"], range_km, delta = 1)
        self.assertAlmostEqual(obsv_metrics["incidence angle [deg]"], 2*abs(poi_lat_deg), delta = 0.15)
        self.assertAlmostEqual(obsv_metrics["look angle [deg]"], abs(poi_lat_deg), delta = 0.15)

        # test negative look angle, 0 deg orbit inclination 
        poi_lat_deg = random.uniform(10, 45)
        range_km = RE # fix range to RE, an isosceles triangle forms
        alt = np.sqrt(RE*RE*(2-2*np.cos(np.deg2rad(180-2*abs(poi_lat_deg))))) - RE
        SpacecraftOrbitState = {'time [JDUT1]':epoch_JDUT1, 'x [km]': 6378.137+alt, 'y [km]': 0, 'z [km]': 0, 'vx [km/s]': 0, 'vy [km/s]': -1*orbital_speed(alt*1e-3), 'vz [km/s]': 0} # altitude 500 km
        TargetCoords = {'lat [deg]': poi_lat_deg, 'lon [deg]': 0}
        obsv_metrics = o.calc_data_metrics(SpacecraftOrbitState, TargetCoords)
        self.assertAlmostEqual(obsv_metrics["observation range [km]"], range_km, delta = 1)
        self.assertAlmostEqual(obsv_metrics["incidence angle [deg]"], 2*abs(poi_lat_deg), delta = 0.15)
        self.assertAlmostEqual(obsv_metrics["look angle [deg]"], -1*abs(poi_lat_deg), delta = 0.15)

        # test 90 deg orbit inclination 
        poi_lat_deg = random.uniform(10, 45)
        range_km = RE # fix range to RE, an isosceles triangle forms
        alt = np.sqrt(RE*RE*(2-2*np.cos(np.deg2rad(180-2*abs(poi_lat_deg))))) - RE
        SpacecraftOrbitState = {'time [JDUT1]':epoch_JDUT1, 'x [km]': 6378.137+alt, 'y [km]': 0, 'z [km]': 0, 'vx [km/s]': 0, 'vy [km/s]': 0, 'vz [km/s]': orbital_speed(alt*1e-3)} # altitude 500 km
        TargetCoords = {'lat [deg]': poi_lat_deg, 'lon [deg]': 0}
        obsv_metrics = o.calc_data_metrics(SpacecraftOrbitState, TargetCoords)
        self.assertAlmostEqual(obsv_metrics["observation range [km]"], range_km, delta = 1)
        self.assertAlmostEqual(obsv_metrics["incidence angle [deg]"], 2*abs(poi_lat_deg), delta = 0.15)
        self.assertAlmostEqual(obsv_metrics["look angle [deg]"], abs(poi_lat_deg), delta = 0.15)

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

    def test_get_pixel_config(self): #@TODO
        o = BasicSensorModel.from_json('{}')
        self.assertEqual(o.get_pixel_config().numberDetectorRows,4)
        self.assertEqual(o.get_pixel_config().numberDetectorCols,4)

        o = BasicSensorModel.from_json('{"numberDetectorRows":5, "numberDetectorCols":10}')
        self.assertEqual(o.get_pixel_config().numberDetectorRows,5)
        self.assertEqual(o.get_pixel_config().numberDetectorCols,10)

        o = BasicSensorModel.from_json('{"numberDetectorCols":10}')
        self.assertEqual(o.get_pixel_config().numberDetectorRows,4)
        self.assertEqual(o.get_pixel_config().numberDetectorCols,10)
    
    def test_synthesize_observation(self): #@TODO
        pass
    
    def test_get_pointing_option(self): #TODO
        pass
