import argparse
import glob
import numpy as np
import re
import sys


def calc(path):
    pattern = path + "/logs/*.o*"
    filenames = glob.glob(pattern)

    dat_coll = []
    dat_corr = []
    for fname in filenames:
        fh = open(fname, "r")
        for line in fh:
            if re.search("^Executed correction events", line):
                num = ''.join(filter(lambda x: x.isdigit(), line))
                dat_corr.append(int(num))
            elif re.search("^Executed collision events", line):
                num = ''.join(filter(lambda x: x.isdigit(), line))
                dat_coll.append(int(num))

    print("path: ", path)
    print("coll: ", dat_coll)
    print("corr: ", dat_corr)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate mean msd of various files.')
    simgroup = parser.add_argument_group('path specific')
    simgroup.add_argument('--path', metavar='path', help='Analyze all appropriate files in given path.', required=True)
    args = parser.parse_args()

    calc(path=args.path)
