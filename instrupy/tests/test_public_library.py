"""Unit tests for instrupy.public_library module.
"""

import unittest
import json
import numpy
import sys, os


from instrupy.public_library import Instrument

class TestPublicLibrary(unittest.TestCase):

    def test_get_coverage_specs(self):
        o = Instrument.from_json('{"@type": "Basic Sensor", "name": "Atom","acronym":"At","mass":10,"volume":12.45, "dataRate": 40, "power": 12, "orientation":{"convention": "XYZ","xRotation":10,"yRotation":20.65,"zRotation":-20}, "fieldOfView": {"sensorGeometry": "CONICAL", "fullConeAngle": 10 }}')
        cov = json.loads(o.get_coverage_specs())
        self.assertEqual(cov['Orientation']['eulerAngle1'], 10)
        self.assertEqual(cov['Orientation']['eulerAngle2'], 20.65)
        self.assertEqual(cov['Orientation']['eulerAngle3'], 340)
        self.assertEqual(cov['Orientation']['eulerSeq1'], 1)
        self.assertEqual(cov['Orientation']['eulerSeq2'], 2)
        self.assertEqual(cov['Orientation']['eulerSeq3'], 3)
        self.assertEqual(cov['fieldOfView']['geometry'], "CONICAL")
        self.assertEqual(cov['fieldOfView']['coneAnglesVector'], [5])
        self.assertEqual(cov['fieldOfView']['clockAnglesVector'], None)
        self.assertEqual(cov['fieldOfView']['AlongTrackFov'], 10)
        self.assertEqual(cov['fieldOfView']['CrossTrackFov'], 10)

        firesat_json = ('{"@type": "Passive Optical Scanner",'
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
                        '   "sideLookAngle": 20'
                        ' },'
                        '"dataRate": 85,'
                        '"numberDetectorRowsAT": 256,'
                        '"numberDetectorColsCT": 1,'
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
                        '"snrThreshold": 10}')
        o = Instrument.from_json(firesat_json)
        cov = json.loads(o.get_coverage_specs())
        self.assertEqual(cov['Orientation']['eulerAngle1'], 0)
        self.assertEqual(cov['Orientation']['eulerAngle2'], 20)
        self.assertEqual(cov['Orientation']['eulerAngle3'], 0)
        self.assertEqual(cov['Orientation']['eulerSeq1'], 1)
        self.assertEqual(cov['Orientation']['eulerSeq2'], 2)
        self.assertEqual(cov['Orientation']['eulerSeq3'], 3)
        self.assertEqual(cov['fieldOfView']['geometry'], "RECTANGULAR")
        numpy.testing.assert_almost_equal(cov['fieldOfView']['coneAnglesVector'], [57.90053973274467, 57.90053973274467, 57.90053973274467, 57.90053973274467], 7)
        numpy.testing.assert_almost_equal(cov['fieldOfView']['clockAnglesVector'], [0.37066537305415925, 179.62933462694585, 180.37066537305415, 359.6293346269458], 7)
        self.assertAlmostEqual(cov['fieldOfView']['AlongTrackFov'], 0.628, places = 7)
        self.assertAlmostEqual(cov['fieldOfView']['CrossTrackFov'], 115.8, places = 7)

        microxsar_json =  ( '{"@type": "Synthetic Aperture Radar",'
                            '"name": "MiroXSAR",'  
                            '"mass": 130,' 
                            '"volume": 0.343,' 
                            '"power": 1100,' 
                            '"orientation": {'
                            '    "convention": "SIDE_LOOK",'
                            '    "sideLookAngle": 30'
                            '},'
                            '"dataRate": 2000,'
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

        o = Instrument.from_json(microxsar_json)
        cov = json.loads(o.get_coverage_specs())
        self.assertEqual(cov['Orientation']['eulerAngle1'], 0)
        self.assertEqual(cov['Orientation']['eulerAngle2'], 30)
        self.assertEqual(cov['Orientation']['eulerAngle3'], 0)
        self.assertEqual(cov['Orientation']['eulerSeq1'], 1)
        self.assertEqual(cov['Orientation']['eulerSeq2'], 2)
        self.assertEqual(cov['Orientation']['eulerSeq3'], 3)
        self.assertEqual(cov['fieldOfView']['geometry'], "RECTANGULAR")
        self.assertEqual(cov['fieldOfView']['coneAnglesVector'], [1.2843229275932238, 1.2843229275932238, 1.2843229275932238, 1.2843229275932238])
        self.assertEqual(cov['fieldOfView']['clockAnglesVector'], [8.130787572398706, 171.8692124276013, 188.1307875723987, 351.86921242760127])
        self.assertAlmostEqual(cov['fieldOfView']['AlongTrackFov'], 0.36326197, places = 7)
        self.assertAlmostEqual(cov['fieldOfView']['CrossTrackFov'], 2.54283383, places = 7)