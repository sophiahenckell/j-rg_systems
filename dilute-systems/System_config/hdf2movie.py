#!/usr/bin/env python3

import argparse

from bedsim.btools.animate import Animate 


if __name__ == '__main__':
#     parser = argparse.ArgumentParser(description='Turn pyBDsim h5 files into movies.')
#     simgroup = parser.add_argument_group('Movie settings')
#     simgroup.add_argument('--filename', metavar='path', help='File generate movie of.', required=True)
#     args = parser.parse_args()
#     
    Animate("output/Ellipse-Benchmark-n1024-k1.8-phi0.3-id1-equilibrate.h5", output_filename="output/testframes/test.mp4")

#     Animate(input_filename=args.filename, output_filename=args.filename+".mp4")
