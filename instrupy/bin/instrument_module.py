import argparse
import os
import sys
import csv
import glob
import json
import instrupy
import pandas
from instrupy.public_library import Instrument


def main(arch_dir):

    poi_filepath = os.path.join(arch_dir, 'poi.csv')
    _AccessInfo_files = glob.glob(arch_dir+"/*_accessInfo.csv")
    _instruSpecs_files = glob.glob(arch_dir+"/*_accessInfo.json")

    indx = 0
    level0dataMetrics_filepath = []
    for AccessInfo_filepath in _AccessInfo_files:

        with open(_instruSpecs_files[indx], 'r') as instru_specs_file:
            instru_specs = instrupy.util.FileUtilityFunctions.from_json(instru_specs_file)

        x = Instrument.from_json(instru_specs)

        ''' Extract the satellite index as written in the Access filename. '''
        temp_AccessInfo_filepath = AccessInfo_filepath.split(os.path.sep)
        temp_last_index = len(temp_AccessInfo_filepath) - 1
        satIndx = temp_AccessInfo_filepath[temp_last_index].split('_')[0]
        sat_filename = str(satIndx) + '_level0_data_metrics.csv'
        level0dataMetrics_filepath.append(os.path.join(arch_dir, sat_filename))
        x.generate_level0_data_metrics(poi_filepath, AccessInfo_filepath,
                                       level0dataMetrics_filepath[indx])
        indx = indx + 1

    level1dataMetrics_filepath = os.path.join(arch_dir,
                                              'level1_data_metrics.csv')
    x.generate_level1_data_metrics(level0dataMetrics_filepath,
                                   level1dataMetrics_filepath)

    level2dataMetrics_filepath = os.path.join(arch_dir,
                                              'level2_data_metrics.csv')
    x.generate_level2_data_metrics(level1dataMetrics_filepath,
                                   level2dataMetrics_filepath)

    """ ''' Generate coverage metrics '''
    # Get mission duration from one of the access files
    # Note that 4th row contains the mission duration
    missionDuration_days = pandas.read_csv(_AccessInfo_files[0],
                                           skiprows=[0, 1, 2],
                                           nrows=1,
                                           header=None)

    missionDuration_days = float(missionDuration_days[0][0].split()[4])

    level1CoverageMetrics_filepath = os.path.join(
            arch_dir,
            'level1_coverage_metrics.csv')
    level2CoverageMetrics_filepath = os.path.join(
            arch_dir,
            'level2_coverage_metrics.csv')
    x.generate_level1_coverage_metrics(missionDuration_days,
                                       level0dataMetrics_filepath,
                                       level1CoverageMetrics_filepath)
    x.generate_level2_data_metrics(level1CoverageMetrics_filepath,
                                   level2CoverageMetrics_filepath) """


class readable_dir(argparse.Action):
    """Defines a custom argparse Action to identify a readable directory."""
    def __call__(self, parser, namespace, values, option_string=None):
        prospective_dir = values
        if not os.path.isdir(prospective_dir):
            raise argparse.ArgumentTypeError(
                '{0} is not a valid path'.format(prospective_dir)
            )
        if os.access(prospective_dir, os.R_OK):
            setattr(namespace, self.dest, prospective_dir)
        else:
            raise argparse.ArgumentTypeError(
                '{0} is not a readable dir'.format(prospective_dir)
            )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Run instrument analysis'
    )
    parser.add_argument(
        'archdir',
        action=readable_dir,
        help="Architecture directory to read inputs/write outputs"
    )
    args = parser.parse_args()
    main(args.archdir)
