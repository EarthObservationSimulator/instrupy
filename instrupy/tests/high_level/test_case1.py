"""Unit test for InstruPy end-to-end execution
"""

import unittest
import json
import numpy
import os
import sys
from instrupy.public_library import Instrument


class HighLevelTestCase1(unittest.TestCase):
    """ Test for all three types of instruments. Currently the only test is to
        see that there is execution without any errors or exceptions being
        thrown.
    """

    def setUp(self):
        self.dir_path = os.path.dirname(os.path.realpath(__file__))
        self.instrument_json_dir = os.path.join(self.dir_path, 'data',
                                                'instrument_json_files')
        ORBIT_DATA_PATH = os.path.join(self.dir_path, 'data',
                                       'orbit_data', 'case1')
        self.poi_filepath = os.path.join(ORBIT_DATA_PATH, 'poi.csv')
        self.AccessInfo_filepath = [
                os.path.join(ORBIT_DATA_PATH, 'sat1_accessInfo.csv'),
                os.path.join(ORBIT_DATA_PATH, 'sat2_accessInfo.csv'),
                os.path.join(ORBIT_DATA_PATH, 'sat3_accessInfo.csv')
        ]
        self.level0dataMetrics_filepath = [
            os.path.join(ORBIT_DATA_PATH, 'sat1_level0_data_metrics.csv'),
            os.path.join(ORBIT_DATA_PATH, 'sat2_level0_data_metrics.csv'),
            os.path.join(ORBIT_DATA_PATH, 'sat3_level0_data_metrics.csv')
        ]
        self.level1dataMetrics_filepath = os.path.join(
                ORBIT_DATA_PATH, 'level1_data_metrics.csv')

        self.level2dataMetrics_filepath = os.path.join(
                ORBIT_DATA_PATH, 'level2_data_metrics.csv')

    def tearDown(self):
        os.remove(self.level0dataMetrics_filepath[0])
        os.remove(self.level0dataMetrics_filepath[1])
        os.remove(self.level0dataMetrics_filepath[2])
        os.remove(self.level1dataMetrics_filepath)
        os.remove(self.level2dataMetrics_filepath)

    def produce_data_metrics(self, instru_specs_filepath):
        with open(instru_specs_filepath, 'r') as instru_specs:
            _sensor = Instrument.from_json(instru_specs)

        _sensor.generate_level0_data_metrics(
                self.poi_filepath,
                self.AccessInfo_filepath[0],
                self.level0dataMetrics_filepath[0]
        )
        _sensor.generate_level0_data_metrics(
                self.poi_filepath,
                self.AccessInfo_filepath[1],
                self.level0dataMetrics_filepath[1]
        )
        _sensor.generate_level0_data_metrics(
                self.poi_filepath,
                self.AccessInfo_filepath[2],
                self.level0dataMetrics_filepath[2]
        )

        _sensor.generate_level1_data_metrics(self.level0dataMetrics_filepath,
                                             self.level1dataMetrics_filepath)

        _sensor.generate_level2_data_metrics(self.level1dataMetrics_filepath,
                                             self.level2dataMetrics_filepath)

    def test_Basic_Sensor(self):
        self.produce_data_metrics(os.path.join(self.instrument_json_dir,
                                  'atom.json'))

    def test_PassiveOpticalSensor(self):
        self.produce_data_metrics(os.path.join(self.instrument_json_dir,
                                  'firesat.json'))

    def test_SyntheticApertureRadar(self):
        self.produce_data_metrics(os.path.join(self.instrument_json_dir,
                                  'microXsar.json'))
