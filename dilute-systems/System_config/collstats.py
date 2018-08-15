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

    mean_coll = np.mean(dat_coll)
    std_coll  = np.std(dat_coll, ddof=1)

    mean_corr = np.mean(dat_corr)
    std_corr  = np.std(dat_corr, ddof=1)

    print("path: ", path)
    print("mean coll: ", mean_coll, " with stddev ", std_coll)
    print("mean corr: ", mean_corr, " with stddev ", std_corr)
        #MSDs.append(np.loadtxt(fname))
    #mean = np.mean(MSDs, axis=0)
    #std = np.std(MSDs, axis=0, ddof=1)
    #np.savetxt(path + '/mean-msd.dat', mean)
    #np.savetxt(path + '/mean-msd-std.dat', std)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate mean msd of various files.')
    simgroup = parser.add_argument_group('path specific')
    simgroup.add_argument('--path', metavar='path', help='Analyze all appropriate files in given path.', required=True)
    args = parser.parse_args()

    calc(path=args.path)
