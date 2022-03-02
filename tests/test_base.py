"""Unit tests for instrupy.base.
"""
import unittest
import numpy as np
import random
from deepdiff import DeepDiff

from instrupy import InstrumentModelFactory, Instrument
from instrupy.basic_sensor_model import BasicSensorModel
from instrupy.passive_optical_scanner_model import PassiveOpticalScannerModel
from instrupy.synthetic_aperture_radar_model import SyntheticApertureRadarModel
from instrupy.radiometer_model import RadiometerModel
from instrupy.util import SphericalGeometry, Orientation, ViewGeometry, Maneuver 

class TestInstrumentModelFactory(unittest.TestCase):
    
    class DummyNewInstrument:
        def __init__(self, *args, **kwargs):
            pass
            
        def from_dict(self):
            return TestInstrumentModelFactory.DummyNewInstrument()

    def test___init__(self):
        factory = InstrumentModelFactory()

        # test the built-in instrumnet models are registered
        self.assertIn('Basic Sensor', factory._creators)
        self.assertEqual(factory._creators['Basic Sensor'], BasicSensorModel)
        self.assertIn('Passive Optical Scanner', factory._creators)
        self.assertEqual(factory._creators['Passive Optical Scanner'], PassiveOpticalScannerModel)
        self.assertIn('Synthetic Aperture Radar', factory._creators)
        self.assertEqual(factory._creators['Synthetic Aperture Radar'], SyntheticApertureRadarModel)
        self.assertIn('Radiometer', factory._creators)
        self.assertEqual(factory._creators['Radiometer'], RadiometerModel)
    
    def test_register_instrument_model(self):
        factory = InstrumentModelFactory()
        factory.register_instrument_model('New Instrument 2021', TestInstrumentModelFactory.DummyNewInstrument)
        self.assertIn('New Instrument 2021', factory._creators)
        self.assertEqual(factory._creators['New Instrument 2021'], TestInstrumentModelFactory.DummyNewInstrument)
        # test the built-in instrumnet models remain registered
        self.assertIn('Basic Sensor', factory._creators)
        self.assertEqual(factory._creators['Basic Sensor'], BasicSensorModel)
        self.assertIn('Passive Optical Scanner', factory._creators)
        self.assertEqual(factory._creators['Passive Optical Scanner'], PassiveOpticalScannerModel)
        self.assertIn('Synthetic Aperture Radar', factory._creators)
        self.assertEqual(factory._creators['Synthetic Aperture Radar'], SyntheticApertureRadarModel)
        self.assertIn('Radiometer', factory._creators)
        self.assertEqual(factory._creators['Radiometer'], RadiometerModel)

    def test_get_instrument_model(self):
        factory = InstrumentModelFactory()
        # register a dummy instrument model
        factory.register_instrument_model('New Instrument 2021', TestInstrumentModelFactory.DummyNewInstrument)
        # test the instrument model classes can be obtained depending on the input specifications
        # basic sensor model
        specs = {"@type": 'Basic Sensor'} # in practice additional instrument specs shall be present in the dictionary
        bs_model = factory.get_instrument_model(specs)
        self.assertIsInstance(bs_model, BasicSensorModel)
        
        # PassiveOpticalScannerModel, 
        specs = {"@type": 'Passive Optical Scanner', 'scanTechnique': 'PUSHBROOM', "numberDetectorRows":1, "numberDetectorCols":500,  "fieldOfViewGeometry":{"shape":"rectangular", "angleWidth": 5, "angleHeight": 5}} # in practice additional instrument specs shall be present in the dictionary
        ps_model = factory.get_instrument_model(specs)
        self.assertIsInstance(ps_model, PassiveOpticalScannerModel)

        # SyntheticApertureRadarModel
        specs = {"@type": 'Synthetic Aperture Radar', "minimumPRF": 2000, "maximumPRF": 8000, "operatingFrequency": 9.6e9, "antennaHeight":5, "antennaWidth":0.5} # in practice additional instrument specs shall be present in the dictionary
        sar_model = factory.get_instrument_model(specs)
        self.assertIsInstance(sar_model, SyntheticApertureRadarModel)

        # RadiometerModel
        specs = {"@type": 'Radiometer'} # in practice additional instrument specs shall be present in the dictionary
        rad_model = factory.get_instrument_model(specs)
        self.assertIsInstance(rad_model, RadiometerModel)

        # DummyNewInstrument
        specs = {"@type": 'New Instrument 2021'} # in practice additional instrument specs shall be present in the dictionary
        di_model = factory.get_instrument_model(specs)
        self.assertIsInstance(di_model, TestInstrumentModelFactory.DummyNewInstrument)

class TestInstrument(unittest.TestCase):

    bs1 = Instrument.from_json('{"name": "Alpha", "mass":10, "volume":12.45, "dataRate": 40, "bitsPerPixel": 8, "power": 12, \
                                  "orientation": {"referenceFrame": "SC_BODY_FIXED", "convention": "REF_FRAME_ALIGNED"}, \
                                  "fieldOfViewGeometry": {"shape": "CIRCULAR", "diameter":2.5 }, \
                                  "sceneFieldOfViewGeometry": {"shape": "CIRCULAR", "diameter":5 }, \
                                  "maneuver":{"maneuverType": "CIRCULAR", "diameter":10}, \
                                  "pointingOption": [{"referenceFrame": "NADIR_POINTING", "convention": "XYZ", "xRotation":0, "yRotation":2.5, "zRotation":0}, \
                                                     {"referenceFrame": "NADIR_POINTING", "convention": "XYZ", "xRotation":0, "yRotation":-2.5, "zRotation":0}  \
                                                    ], \
                                  "numberDetectorRows":5, "numberDetectorCols":10, "@id":"bs1", "@type":"Basic Sensor" \
                                  }')

    bs2 = Instrument.from_json('{"name": "Beta", "mass":10, "volume":12.45, "dataRate": 40, "bitsPerPixel": 8, "power": 12, \
                                  "fieldOfViewGeometry": {"shape": "CIRCULAR", "diameter":5 }, \
                                  "maneuver":{"maneuverType": "SINGLE_ROLL_ONLY", "A_rollMin":10, "A_rollMax":15}, \
                                  "mode": [{"@id":101, "orientation": {"referenceFrame": "SC_BODY_FIXED", "convention": "REF_FRAME_ALIGNED"}} \
                                          ], \
                                  "numberDetectorRows":5, "numberDetectorCols":10, "@type":"Basic Sensor" \
                                  }')

    bs3 = Instrument.from_json('{"name": "Gamma", "mass":10, "volume":12.45, "dataRate": 40, "bitsPerPixel": 8, "power": 12, \
                                  "fieldOfViewGeometry": {"shape": "RECTANGULAR", "angleHeight":0.25, "angleWidth":10 }, \
                                  "sceneFieldOfViewGeometry": {"shape": "RECTANGULAR", "angleHeight":5, "angleWidth":10 }, \
                                  "maneuver":{"maneuverType": "Double_Roll_Only", "A_rollMin":10, "A_rollMax":15, "B_rollMin":-15, "B_rollMax":-10}, \
                                  "pointingOption": [{"referenceFrame": "NADIR_POINTING", "convention": "SIDE_LOOK", "sideLookAngle":10}, \
                                                     {"referenceFrame": "NADIR_POINTING", "convention": "SIDE_LOOK", "sideLookAngle":15}  \
                                                    ], \
                                  "mode": [{"@id":0, "orientation": {"referenceFrame": "SC_BODY_FIXED", "convention": "REF_FRAME_ALIGNED"}}, \
                                           {"@id":1, "orientation": {"referenceFrame": "SC_BODY_FIXED", "convention": "SIDE_LOOK", "sideLookAngle": 25}}, \
                                           { "orientation": {"referenceFrame": "SC_BODY_FIXED", "convention": "SIDE_LOOK", "sideLookAngle": -25}}  \
                                          ], \
                                  "numberDetectorRows":5, "numberDetectorCols":10, "@id": "bs3", "@type":"Basic Sensor" \
                                  }')

    def test_from_json_basic_sensor(self):
        # test initialization with no mode specification
        self.assertEqual(TestInstrument.bs1.name, "Alpha")
        self.assertEqual(TestInstrument.bs1._id, "bs1")
        self.assertEqual(TestInstrument.bs1._type, "Basic Sensor")
        self.assertIsInstance(TestInstrument.bs1, Instrument)
        self.assertEqual(len(TestInstrument.bs1.mode), 1)
        self.assertIsInstance(TestInstrument.bs1.mode[0], BasicSensorModel)
        mode0 = TestInstrument.bs1.mode[0]
        self.assertEqual(mode0._id, "0")
        self.assertEqual(mode0.mass, 10)
        self.assertEqual(mode0.volume, 12.45)
        self.assertEqual(mode0.dataRate, 40)
        self.assertEqual(mode0.bitsPerPixel, 8)
        self.assertEqual(mode0.power, 12)
        self.assertEqual(mode0.numberDetectorRows, 5)
        self.assertEqual(mode0.numberDetectorCols, 10)
        self.assertEqual(mode0.orientation, Orientation.from_json('{"referenceFrame": "SC_BODY_FIXED", "convention": "REF_FRAME_ALIGNED"}'))
        self.assertEqual(mode0.fieldOfView, ViewGeometry.from_dict({"orientation": {"referenceFrame": "SC_BODY_FIXED", "convention": "REF_FRAME_ALIGNED"}, 
                                                                    "sphericalGeometry": {"shape": "CIRCULAR", "diameter":2.5 }}))
        self.assertEqual(mode0.sceneFieldOfView, ViewGeometry.from_dict({"orientation": {"referenceFrame": "SC_BODY_FIXED", "convention": "REF_FRAME_ALIGNED"}, 
                                                                         "sphericalGeometry": {"shape": "CIRCULAR", "diameter":5 }}))
        self.assertEqual(mode0.maneuver, Maneuver.from_json('{"maneuverType": "CIRCULAR", "diameter": 10}'))
        self.assertEqual(mode0.fieldOfRegard, [ViewGeometry.from_dict({"orientation": {"referenceFrame": "NADIR_POINTING", "convention": "REF_FRAME_ALIGNED"}, 
                                                                       "sphericalGeometry": {"shape": "CIRCULAR", "diameter":15 }})])        
        self.assertEqual(mode0.pointingOption, [Orientation.from_dict({"referenceFrame": "NADIR_POINTING", "convention": "XYZ", "xRotation":0, "yRotation":2.5, "zRotation":0}),
                                                 Orientation.from_dict({"referenceFrame": "NADIR_POINTING", "convention": "XYZ", "xRotation":0, "yRotation":-2.5, "zRotation":0})])

        # test initialization with single mode specification
        self.assertEqual(TestInstrument.bs2.name, "Beta")
        self.assertIsNotNone(TestInstrument.bs2._id) # a random id shall be assigned
        self.assertEqual(TestInstrument.bs2._type, "Basic Sensor")
        self.assertIsInstance(TestInstrument.bs2, Instrument)
        self.assertEqual(len(TestInstrument.bs2.mode), 1)
        self.assertIsInstance(TestInstrument.bs2.mode[0], BasicSensorModel)
        mode0 = TestInstrument.bs2.mode[0]
        self.assertEqual(mode0._id, 101)
        self.assertEqual(mode0.mass, 10)
        self.assertEqual(mode0.volume, 12.45)
        self.assertEqual(mode0.dataRate, 40)
        self.assertEqual(mode0.bitsPerPixel, 8)
        self.assertEqual(mode0.power, 12)
        self.assertEqual(mode0.numberDetectorRows, 5)
        self.assertEqual(mode0.numberDetectorCols, 10)
        self.assertEqual(mode0.orientation, Orientation.from_json('{"referenceFrame": "SC_BODY_FIXED", "convention": "REF_FRAME_ALIGNED"}'))
        self.assertEqual(mode0.fieldOfView, ViewGeometry.from_dict({"orientation":{"referenceFrame": "SC_BODY_FIXED", "convention": "REF_FRAME_ALIGNED"},
                                                                     "sphericalGeometry":{"shape": "CIRCULAR", "diameter": 5}}))
        self.assertEqual(mode0.sceneFieldOfView, mode0.fieldOfView)
        self.assertEqual(mode0.maneuver, Maneuver.from_json('{"maneuverType": "single_ROLL_ONLY", "A_rollMin": 10, "A_rollMax":15}'))
        self.assertIsNone(mode0.pointingOption)   

        # test initialization with multiple mode specifications
        self.assertEqual(TestInstrument.bs3.name, "Gamma")
        self.assertEqual(TestInstrument.bs3._id, "bs3")
        self.assertEqual(TestInstrument.bs3._type, "Basic Sensor")
        self.assertIsInstance(TestInstrument.bs3, Instrument)
        self.assertEqual(len(TestInstrument.bs3.mode), 3)
        self.assertIsInstance(TestInstrument.bs3.mode[0], BasicSensorModel)
        # mode0
        mode0 = TestInstrument.bs3.mode[0]
        self.assertEqual(mode0._id, 0)
        self.assertEqual(mode0.mass, 10)
        self.assertEqual(mode0.volume, 12.45)
        self.assertEqual(mode0.dataRate, 40)
        self.assertEqual(mode0.bitsPerPixel, 8)
        self.assertEqual(mode0.power, 12)
        self.assertEqual(mode0.numberDetectorRows, 5)
        self.assertEqual(mode0.numberDetectorCols, 10)
        self.assertEqual(mode0.orientation, Orientation.from_json('{"referenceFrame": "SC_BODY_FIXED", "convention": "REF_FRAME_ALIGNED"}'))
        self.assertEqual(mode0.fieldOfView, ViewGeometry.from_dict({"orientation":{"referenceFrame": "SC_BODY_FIXED", "convention": "REF_FRAME_ALIGNED"},
                                                                    "sphericalGeometry":{"shape": "RECTANGULAR", "angleHeight":0.25, "angleWidth":10}}))
        self.assertEqual(mode0.sceneFieldOfView, ViewGeometry.from_dict({"orientation":{"referenceFrame": "SC_BODY_FIXED", "convention": "REF_FRAME_ALIGNED"},
                                                                        "sphericalGeometry":{"shape": "RECTANGULAR", "angleHeight":5, "angleWidth":10}}))
        self.assertEqual(mode0.maneuver, Maneuver.from_json('{"maneuverType": "double_roll_only", "A_rollMin": 10, "A_rollMax":15, "B_rollMin":-15, "B_rollMax":-10}'))

        ddiff = DeepDiff(mode0.fieldOfRegard, 
                         [ViewGeometry.from_dict({"orientation":{"referenceFrame": "NADIR_POINTING", "convention": "SIDE_Look", "sideLookAngle":12.5},
                                                                      "sphericalGeometry":{"shape": "RECTANGULAR", "angleHeight":5, "angleWidth":15}}),
                          ViewGeometry.from_dict({"orientation":{"referenceFrame": "NADIR_POINTING", "convention": "SIDE_Look", "sideLookAngle":-12.5},
                                                                      "sphericalGeometry":{"shape": "RECTANGULAR", "angleHeight":5, "angleWidth":15}})
                         ], 
                         significant_digits=7, ignore_numeric_type_changes=True)
        self.assertEqual(ddiff, {}, msg=ddiff)

        self.assertEqual(mode0.pointingOption, [Orientation.from_dict({"referenceFrame": "NADIR_POINTING", "convention": "SIDE_LOOK", "sideLookAngle":10}),
                                                Orientation.from_dict({"referenceFrame": "NADIR_POINTING", "convention": "SIDE_LOOK", "sideLookAngle":15})]) 
        # mode1
        mode1 = TestInstrument.bs3.mode[1]
        self.assertEqual(mode1._id, 1)
        self.assertEqual(mode1.mass, 10)
        self.assertEqual(mode1.volume, 12.45)
        self.assertEqual(mode1.dataRate, 40)
        self.assertEqual(mode1.bitsPerPixel, 8)
        self.assertEqual(mode1.power, 12)
        self.assertEqual(mode1.numberDetectorRows, 5)
        self.assertEqual(mode1.numberDetectorCols, 10)
        self.assertEqual(mode0.fieldOfView, ViewGeometry.from_dict({"orientation":{"referenceFrame": "SC_BODY_FIXED", "convention": "REF_FRAME_ALIGNED"},
                                                                   "sphericalGeometry":{"shape": "RECTANGULAR", "angleHeight":0.25, "angleWidth":10}}))
        self.assertEqual(mode0.sceneFieldOfView, ViewGeometry.from_dict({"orientation":{"referenceFrame": "SC_BODY_FIXED", "convention": "REF_FRAME_ALIGNED"},
                                                                        "sphericalGeometry":{"shape": "RECTANGULAR", "angleHeight":5, "angleWidth":10}}))
        self.assertEqual(mode0.maneuver, Maneuver.from_json('{"maneuverType": "double_roll_only", "A_rollMin": 10, "A_rollMax":15, "B_rollMin":-15, "B_rollMax":-10}'))
        ddiff = DeepDiff(mode0.fieldOfRegard, 
                         [ViewGeometry.from_dict({"orientation":{"referenceFrame": "NADIR_POINTING", "convention": "SIDE_Look", "sideLookAngle":12.5},
                                                                      "sphericalGeometry":{"shape": "RECTANGULAR", "angleHeight":5, "angleWidth":15}}),
                          ViewGeometry.from_dict({"orientation":{"referenceFrame": "NADIR_POINTING", "convention": "SIDE_Look", "sideLookAngle":-12.5},
                                                                      "sphericalGeometry":{"shape": "RECTANGULAR", "angleHeight":5, "angleWidth":15}})
                         ], 
                         significant_digits=7, ignore_numeric_type_changes=True)
        self.assertEqual(ddiff, {}, msg=ddiff)

        self.assertEqual(mode1.pointingOption, mode0.pointingOption)        
        # mode2
        mode2 = TestInstrument.bs3.mode[2]
        self.assertIsNotNone(mode2._id)
        self.assertEqual(mode2.mass, 10)
        self.assertEqual(mode2.volume, 12.45)
        self.assertEqual(mode2.dataRate, 40)
        self.assertEqual(mode2.bitsPerPixel, 8)
        self.assertEqual(mode2.power, 12)
        self.assertEqual(mode2.numberDetectorRows, 5)
        self.assertEqual(mode2.numberDetectorCols, 10)
        self.assertEqual(mode2.orientation, Orientation.from_json('{"referenceFrame": "SC_BODY_FIXED", "convention": "SIDE_LOOK", "sideLookAngle": -25}'))
        self.assertEqual(mode0.fieldOfView, ViewGeometry.from_dict({"orientation":{"referenceFrame": "SC_BODY_FIXED", "convention": "REF_FRAME_ALIGNED"},
                                                                   "sphericalGeometry":{"shape": "RECTANGULAR", "angleHeight":0.25, "angleWidth":10}}))
        self.assertEqual(mode0.sceneFieldOfView, ViewGeometry.from_dict({"orientation":{"referenceFrame": "SC_BODY_FIXED", "convention": "REF_FRAME_ALIGNED"},
                                                                        "sphericalGeometry":{"shape": "RECTANGULAR", "angleHeight":5, "angleWidth":10}}))
        self.assertEqual(mode0.maneuver, Maneuver.from_json('{"maneuverType": "double_roll_only", "A_rollMin": 10, "A_rollMax":15, "B_rollMin":-15, "B_rollMax":-10}'))
        ddiff = DeepDiff(mode0.fieldOfRegard, 
                         [ViewGeometry.from_dict({"orientation":{"referenceFrame": "NADIR_POINTING", "convention": "SIDE_Look", "sideLookAngle":12.5},
                                                                      "sphericalGeometry":{"shape": "RECTANGULAR", "angleHeight":5, "angleWidth":15}}),
                          ViewGeometry.from_dict({"orientation":{"referenceFrame": "NADIR_POINTING", "convention": "SIDE_Look", "sideLookAngle":-12.5},
                                                                      "sphericalGeometry":{"shape": "RECTANGULAR", "angleHeight":5, "angleWidth":15}})
                         ], 
                         significant_digits=7, ignore_numeric_type_changes=True)
        self.assertEqual(ddiff, {}, msg=ddiff)

        self.assertEqual(mode2.pointingOption, mode0.pointingOption)   
    
    def test_get_type(self):
        self.assertEqual(TestInstrument.bs1.get_type(), 'Basic Sensor')
        self.assertEqual(TestInstrument.bs2.get_type(), 'Basic Sensor')
        self.assertEqual(TestInstrument.bs3.get_type(), 'Basic Sensor')

    def test_get_id(self):
        self.assertEqual(TestInstrument.bs1.get_id(), "bs1")
        self.assertIsNotNone(TestInstrument.bs2.get_id())
        self.assertEqual(TestInstrument.bs3.get_id(), "bs3")

    def test_get_mode_id(self): #@TODO
        pass

    def test_get_mode(self): #@TODO
        pass    
    
    def test_get_field_of_regard(self): #@TODO
        # bs1
        # no input mode-id
        self.assertEqual(TestInstrument.bs1.get_field_of_regard(), [ViewGeometry.from_dict({"orientation":{"referenceFrame": "NADIR_POINTING", "convention": "REF_FRAME_ALIGNED"}, 
                                                                                            "sphericalGeometry": {"shape": "CIRCULAR", "diameter": 15}})])
        # input correct mode-id
        self.assertEqual(TestInstrument.bs1.get_field_of_regard(mode_id="0"), [ViewGeometry.from_dict({"orientation":{"referenceFrame": "NADIR_POINTING", "convention": "REF_FRAME_ALIGNED"},
                                                                                                       "sphericalGeometry": {"shape": "CIRCULAR", "diameter": 15}})])
        # input incorrect mode-id, should default to first mode
        self.assertEqual(TestInstrument.bs1.get_field_of_regard(mode_id="abc"), [ViewGeometry.from_dict({"orientation":{"referenceFrame": "NADIR_POINTING", "convention": "REF_FRAME_ALIGNED"},
                                                                                                         "sphericalGeometry": {"shape": "CIRCULAR", "diameter": 15}})])

        # bs2
        # no input mode-id
        self.assertEqual(TestInstrument.bs2.get_field_of_regard(), [ViewGeometry.from_dict({"orientation":{"referenceFrame": "NADIR_POINTING", "convention": "SIDE_LOOK", "sideLookAngle": 12.5},"sphericalGeometry": {"shape": "rectangular", "angleHeight": 5, "angleWidth":10}})])
        # input correct mode-id
        ddiff = DeepDiff(TestInstrument.bs2.get_field_of_regard(mode_id=101), [ViewGeometry.from_dict({"orientation":{"referenceFrame": "NADIR_POINTING", "convention": "SIDE_LOOK", "sideLookAngle": 12.5},"sphericalGeometry": {"shape": "rectangular", "angleHeight": 5, "angleWidth":10}})], ignore_numeric_type_changes=True)
        self.assertEqual(ddiff, {})

        # input incorrect mode-id, should default to first mode
        self.assertEqual(TestInstrument.bs2.get_field_of_regard(mode_id="abc"), [ViewGeometry.from_dict({"orientation":{"referenceFrame": "NADIR_POINTING", "convention": "SIDE_LOOK", "sideLookAngle": 12.5},"sphericalGeometry": {"shape": "rectangular", "angleHeight": 5, "angleWidth":10}})])

        # bs3, all modes have the same field of regard
        # no input mode-id
        ddiff = DeepDiff(TestInstrument.bs3.get_field_of_regard()[0], 
                         ViewGeometry.from_dict({"orientation":{"referenceFrame": "NADIR_POINTING", "convention": "SIDE_LOOK", "sideLookAngle": 12.5},"sphericalGeometry": {"shape": "RECTANGULAR", "angleHeight":5, "angleWidth":15 }}), 
                         significant_digits=7, ignore_numeric_type_changes=True)
        self.assertEqual(ddiff, {}, msg=ddiff)
        
        ddiff = DeepDiff(TestInstrument.bs3.get_field_of_regard()[1], 
                         ViewGeometry.from_dict({"orientation":{"referenceFrame": "NADIR_POINTING", "convention": "SIDE_LOOK", "sideLookAngle": -12.5},"sphericalGeometry": {"shape": "RECTANGULAR", "angleHeight":5, "angleWidth":15 }}), 
                         significant_digits=7, ignore_numeric_type_changes=True)
        self.assertEqual(ddiff, {}, msg=ddiff)

        # input correct mode-id
        ddiff = DeepDiff(TestInstrument.bs3.get_field_of_regard(mode_id=0)[0], 
                         ViewGeometry.from_dict({"orientation":{"referenceFrame": "NADIR_POINTING", "convention": "SIDE_LOOK", "sideLookAngle": 12.5},"sphericalGeometry": {"shape": "RECTANGULAR", "angleHeight":5, "angleWidth":15 }}), 
                         significant_digits=7, ignore_numeric_type_changes=True)
        self.assertEqual(ddiff, {}, msg=ddiff)

        ddiff = DeepDiff(TestInstrument.bs3.get_field_of_regard(mode_id=0)[1], 
                         ViewGeometry.from_dict({"orientation":{"referenceFrame": "NADIR_POINTING", "convention": "SIDE_LOOK", "sideLookAngle": -12.5},"sphericalGeometry": {"shape": "RECTANGULAR", "angleHeight":5, "angleWidth":15 }}), 
                         significant_digits=7, ignore_numeric_type_changes=True)
        self.assertEqual(ddiff, {}, msg=ddiff)

        # input incorrect mode-id, should default to first mode
        self.assertEqual(TestInstrument.bs3.get_field_of_regard(mode_id='abc'), TestInstrument.bs3.get_field_of_regard(mode_id=0))
        # next mode
        self.assertEqual(TestInstrument.bs3.get_field_of_regard(mode_id=1), TestInstrument.bs3.get_field_of_regard(mode_id=0))
        # next mode,
        mode_id = TestInstrument.bs3.mode_id[2]
        self.assertEqual(TestInstrument.bs3.get_field_of_regard(mode_id=mode_id), TestInstrument.bs3.get_field_of_regard(mode_id=0))
                
    def test_get_field_of_view(self):
        # bs1
        # no input mode-id
        self.assertEqual(TestInstrument.bs1.get_field_of_view(), ViewGeometry.from_dict({"orientation":{"referenceFrame": "SC_BODY_FIXED", "convention": "REF_FRAME_ALIGNED"},"sphericalGeometry": {"shape": "CIRCULAR", "diameter": 2.5}}))
        # input correct mode-id
        self.assertEqual(TestInstrument.bs1.get_field_of_view(mode_id="0"), ViewGeometry.from_dict({"orientation":{"referenceFrame": "SC_BODY_FIXED", "convention": "REF_FRAME_ALIGNED"},"sphericalGeometry": {"shape": "CIRCULAR", "diameter": 2.5}}))
        # input incorrect mode-id, should default to first mode
        self.assertEqual(TestInstrument.bs1.get_field_of_view(mode_id="abc"), ViewGeometry.from_dict({"orientation":{"referenceFrame": "SC_BODY_FIXED", "convention": "REF_FRAME_ALIGNED"},"sphericalGeometry": {"shape": "CIRCULAR", "diameter": 2.5}}))

        # bs2
        # no input mode-id
        self.assertEqual(TestInstrument.bs2.get_field_of_view(), ViewGeometry.from_dict({"orientation":{"referenceFrame": "SC_BODY_FIXED", "convention": "REF_FRAME_ALIGNED"},"sphericalGeometry": {"shape": "CIRCULAR", "diameter": 5}}))
        # input correct mode-id
        self.assertEqual(TestInstrument.bs2.get_field_of_view(mode_id=101), ViewGeometry.from_dict({"orientation":{"referenceFrame": "SC_BODY_FIXED", "convention": "REF_FRAME_ALIGNED"},"sphericalGeometry": {"shape": "CIRCULAR", "diameter": 5}}))
        # input incorrect mode-id, should default to first mode
        self.assertEqual(TestInstrument.bs2.get_field_of_view(mode_id="abc"), ViewGeometry.from_dict({"orientation":{"referenceFrame": "SC_BODY_FIXED", "convention": "REF_FRAME_ALIGNED"},"sphericalGeometry": {"shape": "CIRCULAR", "diameter": 5}}))

        # bs3
        # no input mode-id
        self.assertEqual(TestInstrument.bs3.get_field_of_view(), ViewGeometry.from_dict({"orientation":{"referenceFrame": "SC_BODY_FIXED", "convention": "REF_FRAME_ALIGNED"},"sphericalGeometry": {"shape": "RECTANGULAR", "angleHeight":0.25, "angleWidth":10 }}))
        # input correct mode-id
        self.assertEqual(TestInstrument.bs3.get_field_of_view(mode_id=0), ViewGeometry.from_dict({"orientation":{"referenceFrame": "SC_BODY_FIXED", "convention": "REF_FRAME_ALIGNED"},"sphericalGeometry": {"shape": "RECTANGULAR", "angleHeight":0.25, "angleWidth":10 }}))
        # input incorrect mode-id, should default to first mode
        self.assertEqual(TestInstrument.bs3.get_field_of_view(mode_id="abc"), ViewGeometry.from_dict({"orientation":{"referenceFrame": "SC_BODY_FIXED", "convention": "REF_FRAME_ALIGNED"},"sphericalGeometry": {"shape": "RECTANGULAR", "angleHeight":0.25, "angleWidth":10 }}))
        # next mode
        self.assertEqual(TestInstrument.bs3.get_field_of_view(mode_id=1), ViewGeometry.from_dict({"orientation":{"referenceFrame": "SC_BODY_FIXED", "convention": "SIDE_LOOK", "sideLookAngle": 25},"sphericalGeometry": {"shape": "RECTANGULAR", "angleHeight":0.25, "angleWidth":10 }}))
        # next mode
        mode_id = TestInstrument.bs3.mode_id[2]
        self.assertEqual(TestInstrument.bs3.get_field_of_view(mode_id=mode_id), ViewGeometry.from_dict({"orientation":{"referenceFrame": "SC_BODY_FIXED", "convention": "SIDE_LOOK", "sideLookAngle": -25},"sphericalGeometry": {"shape": "RECTANGULAR", "angleHeight":0.25, "angleWidth":10 }}))
        
    def test_get_scene_field_of_view(self):
        # bs1
        # no input mode-id
        self.assertEqual(TestInstrument.bs1.get_scene_field_of_view(), ViewGeometry.from_dict({"orientation":{"referenceFrame": "SC_BODY_FIXED", "convention": "REF_FRAME_ALIGNED"},"sphericalGeometry": {"shape": "CIRCULAR", "diameter": 5}}))
        # input correct mode-id
        self.assertEqual(TestInstrument.bs1.get_scene_field_of_view(mode_id="0"), ViewGeometry.from_dict({"orientation":{"referenceFrame": "SC_BODY_FIXED", "convention": "REF_FRAME_ALIGNED"},"sphericalGeometry": {"shape": "CIRCULAR", "diameter": 5}}))
        # input incorrect mode-id, should default to first mode
        self.assertEqual(TestInstrument.bs1.get_scene_field_of_view(mode_id="abc"), ViewGeometry.from_dict({"orientation":{"referenceFrame": "SC_BODY_FIXED", "convention": "REF_FRAME_ALIGNED"},"sphericalGeometry": {"shape": "CIRCULAR", "diameter": 5}}))

        # bs2
        # no input mode-id
        self.assertEqual(TestInstrument.bs2.get_scene_field_of_view(), ViewGeometry.from_dict({"orientation":{"referenceFrame": "SC_BODY_FIXED", "convention": "REF_FRAME_ALIGNED"},"sphericalGeometry": {"shape": "CIRCULAR", "diameter": 5}}))
        # input correct mode-id
        self.assertEqual(TestInstrument.bs2.get_scene_field_of_view(mode_id=101), ViewGeometry.from_dict({"orientation":{"referenceFrame": "SC_BODY_FIXED", "convention": "REF_FRAME_ALIGNED"},"sphericalGeometry": {"shape": "CIRCULAR", "diameter": 5}}))
        # input incorrect mode-id, should default to first mode
        self.assertEqual(TestInstrument.bs2.get_scene_field_of_view(mode_id="abc"), ViewGeometry.from_dict({"orientation":{"referenceFrame": "SC_BODY_FIXED", "convention": "REF_FRAME_ALIGNED"},"sphericalGeometry": {"shape": "CIRCULAR", "diameter": 5}}))

        # bs3
        # no input mode-id
        self.assertEqual(TestInstrument.bs3.get_scene_field_of_view(), ViewGeometry.from_dict({"orientation":{"referenceFrame": "SC_BODY_FIXED", "convention": "REF_FRAME_ALIGNED"},"sphericalGeometry": {"shape": "RECTANGULAR", "angleHeight":5, "angleWidth":10 }}))
        # input correct mode-id
        self.assertEqual(TestInstrument.bs3.get_scene_field_of_view(mode_id=0), ViewGeometry.from_dict({"orientation":{"referenceFrame": "SC_BODY_FIXED", "convention": "REF_FRAME_ALIGNED"},"sphericalGeometry": {"shape": "RECTANGULAR", "angleHeight":5, "angleWidth":10 }}))
        # input incorrect mode-id, should default to first mode
        self.assertEqual(TestInstrument.bs3.get_scene_field_of_view(mode_id="abc"), ViewGeometry.from_dict({"orientation":{"referenceFrame": "SC_BODY_FIXED", "convention": "REF_FRAME_ALIGNED"},"sphericalGeometry": {"shape": "RECTANGULAR", "angleHeight":5, "angleWidth":10 }}))
        # next mode
        self.assertEqual(TestInstrument.bs3.get_scene_field_of_view(mode_id=1), ViewGeometry.from_dict({"orientation":{"referenceFrame": "SC_BODY_FIXED", "convention": "SIDE_LOOK", "sideLookAngle": 25},"sphericalGeometry": {"shape": "RECTANGULAR", "angleHeight":5, "angleWidth":10 }}))
        # next mode
        mode_id = TestInstrument.bs3.mode_id[2]
        self.assertEqual(TestInstrument.bs3.get_scene_field_of_view(mode_id=mode_id), ViewGeometry.from_dict({"orientation":{"referenceFrame": "SC_BODY_FIXED", "convention": "SIDE_LOOK", "sideLookAngle": -25},"sphericalGeometry": {"shape": "RECTANGULAR", "angleHeight":5, "angleWidth":10 }}))
       

    def test_get_orientation(self):
        # bs1
        # no input mode-id
        self.assertEqual(TestInstrument.bs1.get_orientation(), Orientation.from_json('{"referenceFrame": "SC_BODY_FIXED", "convention": "REF_FRAME_ALIGNED"}'))
        # input correct mode-id
        self.assertEqual(TestInstrument.bs1.get_orientation(mode_id="0"), Orientation.from_json('{"referenceFrame": "SC_BODY_FIXED", "convention": "REF_FRAME_ALIGNED"}'))
        # input incorrect mode-id, should default to first mode
        self.assertEqual(TestInstrument.bs1.get_orientation(mode_id="abc"), Orientation.from_json('{"referenceFrame": "SC_BODY_FIXED", "convention": "REF_FRAME_ALIGNED"}'))

        # bs2
        # no input mode-id
        self.assertEqual(TestInstrument.bs2.get_orientation(), Orientation.from_json('{"referenceFrame": "SC_BODY_FIXED", "convention": "REF_FRAME_ALIGNED"}'))
        # input correct mode-id
        self.assertEqual(TestInstrument.bs2.get_orientation(mode_id=101), Orientation.from_json('{"referenceFrame": "SC_BODY_FIXED", "convention": "REF_FRAME_ALIGNED"}'))
        # input incorrect mode-id, should default to first mode
        self.assertEqual(TestInstrument.bs2.get_orientation(mode_id="abc"), Orientation.from_json('{"referenceFrame": "SC_BODY_FIXED", "convention": "REF_FRAME_ALIGNED"}'))

        # bs3
        # no input mode-id
        self.assertEqual(TestInstrument.bs3.get_orientation(), Orientation.from_json('{"referenceFrame": "SC_BODY_FIXED", "convention": "REF_FRAME_ALIGNED"}'))
        # input correct mode-id
        self.assertEqual(TestInstrument.bs3.get_orientation(mode_id=0), Orientation.from_json('{"referenceFrame": "SC_BODY_FIXED", "convention": "REF_FRAME_ALIGNED"}'))
        # input incorrect mode-id, should default to first mode
        self.assertEqual(TestInstrument.bs3.get_orientation(mode_id="abc"), Orientation.from_json('{"referenceFrame": "SC_BODY_FIXED", "convention": "REF_FRAME_ALIGNED"}'))
        # next mode
        self.assertEqual(TestInstrument.bs3.get_orientation(mode_id=1), Orientation.from_json('{"referenceFrame": "SC_BODY_FIXED", "convention": "SIDE_LOOK", "sideLookAngle": 25}'))
        # next mode
        mode_id = TestInstrument.bs3.mode_id[2]
        self.assertEqual(TestInstrument.bs3.get_orientation(mode_id=mode_id), Orientation.from_json('{"referenceFrame": "SC_BODY_FIXED", "convention": "SIDE_LOOK", "sideLookAngle": -25}'))
        

    def test_pixel_config(self): #@TODO
        pass

    def test_calc_data_metrics_bs1(self): #@TODO
        """ Simple test involving satellite above POI at (lat = 0,lon = 0). Date chosen so that ECEF and ECI frames are aligned.
            Sensor specs do not influence the below calcs. They do however shall influence the coverage calcs (which is not covered by this test).
        """
        epoch_JDUT1 =  2458543.06088 # 2019 Feb 28 13:27:40 is time at which the ECEF and ECI frames approximately align, hence ECEF to ECI rotation is identity. See <https://www.celnav.de/longterm.htm> online calculator of GMST.
        
        SpacecraftOrbitState = {'time [JDUT1]':epoch_JDUT1, 'x [km]': 6878.137, 'y [km]': 0, 'z [km]': 0, 'vx [km/s]': 0, 'vy [km/s]': 7.6126, 'vz [km/s]': 0} # altitude 500 km
        TargetCoords = {'lat [deg]': 0, 'lon [deg]': 0}
        
        # no input mode-id
        obsv_metrics = TestInstrument.bs1.calc_data_metrics(None, SpacecraftOrbitState, TargetCoords)
        self.assertAlmostEqual(obsv_metrics["observation range [km]"], 500, delta = 1)
        self.assertAlmostEqual(obsv_metrics["incidence angle [deg]"], 0, delta = 0.1)
        self.assertAlmostEqual(obsv_metrics["look angle [deg]"], 0, delta = 0.1)
        self.assertAlmostEqual(obsv_metrics["solar zenith [deg]"], 20.335, delta = 0.1) # precomputed value at the epoch and (lat=0, lon=0) position

        # correct mode-id
        obsv_metrics = TestInstrument.bs1.calc_data_metrics("0", SpacecraftOrbitState, TargetCoords)
        self.assertAlmostEqual(obsv_metrics["observation range [km]"], 500, delta = 1)
        self.assertAlmostEqual(obsv_metrics["incidence angle [deg]"], 0, delta = 0.1)
        self.assertAlmostEqual(obsv_metrics["look angle [deg]"], 0, delta = 0.1)
        self.assertAlmostEqual(obsv_metrics["solar zenith [deg]"], 20.335, delta = 0.1) # precomputed value at the epoch and (lat=0, lon=0) position

        # incorrect mode-id
        obsv_metrics = TestInstrument.bs1.calc_data_metrics("abc", SpacecraftOrbitState, TargetCoords)
        self.assertAlmostEqual(obsv_metrics["observation range [km]"], 500, delta = 1)
        self.assertAlmostEqual(obsv_metrics["incidence angle [deg]"], 0, delta = 0.1)
        self.assertAlmostEqual(obsv_metrics["look angle [deg]"], 0, delta = 0.1)
        self.assertAlmostEqual(obsv_metrics["solar zenith [deg]"], 20.335, delta = 0.1) # precomputed value at the epoch and (lat=0, lon=0) position

    def test_synthesize_observation(self): #@TODO
        pass

    def test_get_pointing_option(self): #@TODO
        pass  
    



