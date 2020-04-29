"""Unit tests for instrupy.basic_sensor module.
"""

import unittest
import json
import numpy
import sys, os


from instrupy.basic_sensor import *
from instrupy.util import Orientation, FieldOfView

class TestBasicSensor(unittest.TestCase):

    def test_from_json_basic(self):
        
        # Test: Typical case
        o = BasicSensor.from_json('{"name": "Atom","acronym":"At","mass":10,"volume":12.45, "dataRate": 40, "bitsPerPixel": 8, "power": 12, "orientation":{"convention": "XYZ","xRotation":0,"yRotation":0,"zRotation":0}, "fieldOfView": {"sensorGeometry": "CUSTOM", "customConeAnglesVector": 10 }}')
        self.assertEqual(o._type, "Basic Sensor")
        self.assertEqual(o.name, "Atom")
        self.assertIsInstance(o.name, str)
        self.assertEqual(o.acronym, "At")
        self.assertIsInstance(o.acronym, str)
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
        self.assertEqual(o.orientation.euler_angle2, 0)
        self.assertIsInstance(o.fieldOfView, FieldOfView)
        self.assertEqual(o.fieldOfView._coneAngleVec_deg, [10])
        self.assertIsInstance(o, BasicSensor)
        self.assertIsNone(o._id)
        self.assertEqual(o._type, "Basic Sensor")

        # Test: Case with missing fields
        o = BasicSensor.from_json('{"name": "Atom","mass":10,"volume":12.45, "fieldOfView": {"sensorGeometry": "CONICAL", "fullConeAngle": 10 }}')
        self.assertEqual(o.acronym, "Atom")
        self.assertIsNone(o.power)
        self.assertIsNone(o.dataRate)
        self.assertIsInstance(o.orientation, Orientation)
        self.assertIsInstance(o.fieldOfView, FieldOfView)
        self.assertEqual(o.fieldOfView._coneAngleVec_deg, [5])
        self.assertIsNone(o.fieldOfView._clockAngleVec_deg)

        # Test: Incomplete field-of-view specification, test that Exception is raised
        with self.assertRaises(Exception):
            BasicSensor.from_json('{"name": "Atom","mass":10,"volume":12.45, "fieldOfView": {"sensorGeometry": "RECTANGULAR", "alongTrackFieldOfView": 10 }}')


    def test_calc_typ_data_metrics(self):

        o = BasicSensor.from_json('{"name": "Atom","orientation":{"convention": "SIDE_LOOK", "sideLookAngle":22.5}, "fieldOfView": {"sensorGeometry": "CONICAL", "fullConeAngle": 60 }}')
        epoch_JDUT1 =  2458543.06088 # 2019 Feb 28 13:27:40 is time at which the ECEF and ECI frames approximately align, hence ECEF to ECI rotation is identity. See <https://www.celnav.de/longterm.htm> online calculator of GMST.
        
        # Test: Check a simple intuitive scenario when satellite is 500 km above POI at (lat = 0,lon = 0)
        SpacecraftOrbitState = {'Time[JDUT1]':epoch_JDUT1, 'x[km]': 6878.137, 'y[km]': 0, 'z[km]': 0, 'vx[km/s]': 0, 'vy[km/s]': 7.613, 'vz[km/s]': 0} # equatorial orbit, altitude about 500 km
        TargetCoords = {'Lat [deg]': 0, 'Lon [deg]': 0}
        obsv_metrics = o.calc_typ_data_metrics(SpacecraftOrbitState, TargetCoords)
        self.assertTrue(obsv_metrics["Coverage [T/F]"])
        self.assertAlmostEqual(obsv_metrics["Observation Range [km]"], 500, delta = 1)
        self.assertAlmostEqual(obsv_metrics["Incidence angle [deg]"], 0, delta = 0.1)
        self.assertAlmostEqual(obsv_metrics["Look angle [deg]"], 0, delta = 0.1)

        # Test: Check a simple intuitive scenario when satellite is at altitude of 500 km and makes observation of point "ahead" of it in direction of velocity vector.
        poi_lon_deg = 0.1
        nadir_angle_deg = abs(numpy.rad2deg(6378.137*numpy.deg2rad(poi_lon_deg)/ 500))
        range_km = 500/numpy.cos(numpy.deg2rad(nadir_angle_deg))
        SpacecraftOrbitState = {'Time[JDUT1]':epoch_JDUT1, 'x[km]': 6878.137, 'y[km]': 0, 'z[km]': 0, 'vx[km/s]': 0, 'vy[km/s]': 7.613, 'vz[km/s]': 0} # equatorial orbit, altitude about 500 km
        TargetCoords = {'Lat [deg]': 0, 'Lon [deg]': 0.1}
        obsv_metrics = o.calc_typ_data_metrics(SpacecraftOrbitState, TargetCoords)
        self.assertTrue(obsv_metrics["Coverage [T/F]"])
        self.assertAlmostEqual(obsv_metrics["Observation Range [km]"], range_km, delta = 1)
        self.assertAlmostEqual(obsv_metrics["Incidence angle [deg]"], nadir_angle_deg, delta = 0.1)
        self.assertAlmostEqual(obsv_metrics["Look angle [deg]"], nadir_angle_deg, delta = 0.1)

        # Test: Check a simple intuitive scenario when satellite is at altitude of 500 km and makes observation of point "behind" it opposite to direction of velocity vector.
        poi_lon_deg = -1.2
        nadir_angle_deg = abs(numpy.rad2deg(6378.137*numpy.deg2rad(poi_lon_deg)/ 500))
        range_km = 500/numpy.cos(numpy.deg2rad(nadir_angle_deg))
        SpacecraftOrbitState = {'Time[JDUT1]':epoch_JDUT1, 'x[km]': 6878.137, 'y[km]': 0, 'z[km]': 0, 'vx[km/s]': 0, 'vy[km/s]': 7.613, 'vz[km/s]': 0} # equatorial orbit, altitude about 500 km
        TargetCoords = {'Lat [deg]': 0, 'Lon [deg]': -1.2}
        obsv_metrics = o.calc_typ_data_metrics(SpacecraftOrbitState, TargetCoords)
        self.assertTrue(obsv_metrics["Coverage [T/F]"])
        self.assertAlmostEqual(obsv_metrics["Observation Range [km]"], range_km, delta = 1)
        self.assertAlmostEqual(obsv_metrics["Incidence angle [deg]"], nadir_angle_deg, delta = 1) # larger deltas since truth data is approximate calculation (assumes flat-Earth)
        self.assertAlmostEqual(obsv_metrics["Look angle [deg]"], nadir_angle_deg, delta = 1)

        # Test: Check a simple intuitive scenario when satellite is at altitude of 500 km and makes observation of point "right" of the direction of velocity vector.
        poi_lat_deg = -0.1
        nadir_angle_deg = abs(numpy.rad2deg(6378.137*numpy.deg2rad(poi_lat_deg)/ 500))
        range_km = 500/numpy.cos(numpy.deg2rad(nadir_angle_deg))
        SpacecraftOrbitState = {'Time[JDUT1]':epoch_JDUT1, 'x[km]': 6878.137, 'y[km]': 0, 'z[km]': 0, 'vx[km/s]': 0, 'vy[km/s]': 7.613, 'vz[km/s]': 0} # equatorial orbit, altitude about 500 km
        TargetCoords = {'Lat [deg]': -0.1, 'Lon [deg]': 0}
        obsv_metrics = o.calc_typ_data_metrics(SpacecraftOrbitState, TargetCoords)
        self.assertTrue(obsv_metrics["Coverage [T/F]"])
        self.assertAlmostEqual(obsv_metrics["Observation Range [km]"], range_km, delta = 1)
        self.assertAlmostEqual(obsv_metrics["Incidence angle [deg]"], nadir_angle_deg, delta = 0.15)
        self.assertAlmostEqual(obsv_metrics["Look angle [deg]"], nadir_angle_deg, delta = 0.15)

        # Test: Check a simple intuitive scenario when satellite is at altitude of 500 km and makes observation of point "left" of the direction of velocity vector.
        poi_lat_deg = 1.2
        nadir_angle_deg = abs(numpy.rad2deg(6378.137*numpy.deg2rad(poi_lat_deg)/ 500))
        range_km = 500/numpy.cos(numpy.deg2rad(nadir_angle_deg))
        SpacecraftOrbitState = {'Time[JDUT1]':epoch_JDUT1, 'x[km]': 6878.137, 'y[km]': 0, 'z[km]': 0, 'vx[km/s]': 0, 'vy[km/s]': 7.613, 'vz[km/s]': 0} # equatorial orbit, altitude about 500 km
        TargetCoords = {'Lat [deg]': 1.2, 'Lon [deg]': 0}
        obsv_metrics = o.calc_typ_data_metrics(SpacecraftOrbitState, TargetCoords)
        self.assertTrue(obsv_metrics["Coverage [T/F]"])
        self.assertAlmostEqual(obsv_metrics["Observation Range [km]"], range_km, delta = 1)
        self.assertAlmostEqual(obsv_metrics["Incidence angle [deg]"], nadir_angle_deg, delta = 1)
        self.assertAlmostEqual(obsv_metrics["Look angle [deg]"], nadir_angle_deg, delta = 1)
