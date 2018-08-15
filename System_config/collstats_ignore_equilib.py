import argparse
import glob
import numpy as np
import re
import sys


def calc(path):
    pattern = path + "/logs/*.o*"
    filenames = glob.glob(pattern)

    dat_coll = {}
    dat_corr = {}
    for fname in filenames:
        found_coll = 0
        found_corr = 0
        fh = open(fname, "r")
        for line in fh:
            if re.search("^Executed correction events", line):
                if str(found_corr) not in dat_corr:
                    dat_corr[str(found_corr)] = []
                num = ''.join(filter(lambda x: x.isdigit(), line))
                dat_corr[str(found_corr)].append(int(num))
                found_corr += 1
            elif re.search("^Executed collision events", line):
                if str(found_coll) not in dat_coll:
                    dat_coll[str(found_coll)] = []
                num = ''.join(filter(lambda x: x.isdigit(), line))
                dat_coll[str(found_coll)].append(int(num))
                found_coll += 1
    #print("coll: ", dat_coll)

    mean_coll = [np.mean(dat_coll[k]) for k in dat_coll]
    std_coll  = [np.std(dat_coll[k], ddof=1) for k in dat_coll]

    mean_corr = [np.mean(dat_corr[k]) for k in dat_corr]
    std_corr  = [np.std(dat_corr[k], ddof=1) for k in dat_corr]

    print("path: ", path)
    print("mean coll: ", mean_coll, " with stddev ", std_coll)
    print("mean corr: ", mean_corr, " with stddev ", std_corr)
    print("last mean coll: ", mean_coll[-1], " with stddev ", std_coll[-1])
    print("last mean corr: ", mean_corr[-1], " with stddev ", std_corr[-1])

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
