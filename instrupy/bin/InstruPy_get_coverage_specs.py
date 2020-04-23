import argparse
import os, sys
import csv
import glob
import json
import instrupy

from instrupy.public_library import Instrument

def main(input_json):
    '''
        ##### Deprecated #####
        This module allows command line execution of the python function to extract from json input instrument specifications, OC understandable
        coverage specifications (also a json string).

        example usage: 
        `python bin/instrupy_get_coverage_specs.py "{\"@type\": \"Basic Sensor\", \"name\": \"Atom\",\"acronym\":\"At\",\"mass\":10,\"volume\":12.45, \"dataRate\": 40, \"power\": 12, \"orientation\":{\"convention\": \"XYZ\",\"xRotation\":0,\"yRotation\":0,\"zRotation\":0}, \"fieldOfView\": {\"sensorGeometry\": \"CONICAL\", \"fullConeAngle\": 10 }}"`

        .. note:: The character '\' needs to be prefixed before every double-quote in the json formatted string input (as tested in Windows).
    
    '''

    o = Instrument.from_json(input_json)
    print(o.get_coverage_specs())
    return o.get_coverage_specs()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Translate Instrument high-level specifications to OC understandable Coverage specifications (orientation and field-of-view).'
    )
    parser.add_argument(
        'input_json',
        help = "Instrument specifications in JSON format.",
        type = json.loads
    )
    args = parser.parse_args()
    main(args.input_json)





   