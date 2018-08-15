import argparse
import glob
import numpy as np


def calc(path):
    pattern = path + "/*.h5-msd.dat"
    filenames = glob.glob(pattern)
    MSDs = []
    for fname in filenames:
        MSDs.append(np.loadtxt(fname))
    mean = np.mean(MSDs, axis=0)
    std = np.std(MSDs, axis=0, ddof=1)
    np.savetxt(path + '/mean-msd.dat', mean)
    np.savetxt(path + '/mean-msd-std.dat', std)

def calc_msr(path):
    pattern = path + "/*.h5-msr.dat"
    filenames = glob.glob(pattern)
    MSRs = []
    for fname in filenames:
        MSRs.append(np.loadtxt(fname))
    mean = np.mean(MSRs, axis=0)
    std = np.std(MSRs, axis=0, ddof=1)
    np.savetxt(path + '/mean-msr.dat', mean)
    np.savetxt(path + '/mean-msr-std.dat', std)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate mean msd of various files.')
    simgroup = parser.add_argument_group('path specific')
    simgroup.add_argument('--path', metavar='path', help='Analyze all appropriate files in given path.', required=True)
    args = parser.parse_args()

    calc(path=args.path)
    calc_msr(path=args.path)
