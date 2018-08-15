#!/usr/bin/env python3

import argparse

from bedsim.btools.statprop import Statprop


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate statistics of apyBDsim h5 file.')
    simgroup = parser.add_argument_group('Statistics settings')
    simgroup.add_argument('--filename', metavar='path', help='File to analyze.', required=True)
    args = parser.parse_args()

    Statprop(input_filename=args.filename)
