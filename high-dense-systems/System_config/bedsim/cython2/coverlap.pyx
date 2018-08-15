import numpy as np
from matplotlib.mathtext import DELTA
cimport numpy as np
cimport cython

cdef extern from "coverlap.h": 
    long double f_AB(long double aA, long double aB, long double bA, long double bB, long double cA, long double cB, long double qA0, long double qA1, long double qA2, long double qA3, long double qB0, long double qB1, long double qB2, long double qB3, long double rA0, long double rA1, long double rA2, long double rB0, long double rB1, long double rB2, long double lambda0);
    long double fAB(long double aA, long double aB, long double bA, long double bB, long double cA, long double cB, long double qA0, long double qA1, long double qA2, long double qA3, long double qB0, long double qB1, long double qB2, long double qB3, long double rA0, long double rA1, long double rA2, long double rB0, long double rB1, long double rB2, long double lambda0);
    long double hAB(long double aA, long double aB, long double bA, long double bB, long double cA, long double cB, long double qA0, long double qA1, long double qA2, long double qA3, long double qB0, long double qB1, long double qB2, long double qB3, long double rA0, long double rA1, long double rA2, long double rB0, long double rB1, long double rB2, long double lambda0);
    long double dhAB(long double aA, long double aB, long double bA, long double bB, long double cA, long double cB, long double qA0, long double qA1, long double qA2, long double qA3, long double qB0, long double qB1, long double qB2, long double qB3, long double rA0, long double rA1, long double rA2, long double rB0, long double rB1, long double rB2, long double lambda0);
    long double rCx(long double aA, long double aB, long double bA, long double bB, long double cA, long double cB, long double qA0, long double qA1, long double qA2, long double qA3, long double qB0, long double qB1, long double qB2, long double qB3, long double rA0, long double rA1, long double rA2, long double rB0, long double rB1, long double rB2, long double lambda0);
    long double rCy(long double aA, long double aB, long double bA, long double bB, long double cA, long double cB, long double qA0, long double qA1, long double qA2, long double qA3, long double qB0, long double qB1, long double qB2, long double qB3, long double rA0, long double rA1, long double rA2, long double rB0, long double rB1, long double rB2, long double lambda0);
    long double rCz(long double aA, long double aB, long double bA, long double bB, long double cA, long double cB, long double qA0, long double qA1, long double qA2, long double qA3, long double qB0, long double qB1, long double qB2, long double qB3, long double rA0, long double rA1, long double rA2, long double rB0, long double rB1, long double rB2, long double lambda0);
    long double normx(long double aA, long double aB, long double bA, long double bB, long double cA, long double cB, long double qA0, long double qA1, long double qA2, long double qA3, long double qB0, long double qB1, long double qB2, long double qB3, long double rA0, long double rA1, long double rA2, long double rB0, long double rB1, long double rB2, long double lambda0);
    long double normy(long double aA, long double aB, long double bA, long double bB, long double cA, long double cB, long double qA0, long double qA1, long double qA2, long double qA3, long double qB0, long double qB1, long double qB2, long double qB3, long double rA0, long double rA1, long double rA2, long double rB0, long double rB1, long double rB2, long double lambda0);
    long double normz(long double aA, long double aB, long double bA, long double bB, long double cA, long double cB, long double qA0, long double qA1, long double qA2, long double qA3, long double qB0, long double qB1, long double qB2, long double qB3, long double rA0, long double rA1, long double rA2, long double rB0, long double rB1, long double rB2, long double lambda0);
    long double vCx(long double aA, long double aB, long double bA, long double bB, long double cA, long double cB, long double qA0, long double qA1, long double qA2, long double qA3, long double qB0, long double qB1, long double qB2, long double qB3, long double rA0, long double rA1, long double rA2, long double rB0, long double rB1, long double rB2, long double vA0, long double vA1, long double vA2, long double vB0, long double vB1, long double vB2, long double omegaA0, long double omegaA1, long double omegaA2, long double omegaB0, long double omegaB1, long double omegaB2, long double lambda0);
    long double vCy(long double aA, long double aB, long double bA, long double bB, long double cA, long double cB, long double qA0, long double qA1, long double qA2, long double qA3, long double qB0, long double qB1, long double qB2, long double qB3, long double rA0, long double rA1, long double rA2, long double rB0, long double rB1, long double rB2, long double vA0, long double vA1, long double vA2, long double vB0, long double vB1, long double vB2, long double omegaA0, long double omegaA1, long double omegaA2, long double omegaB0, long double omegaB1, long double omegaB2, long double lambda0);
    long double vCz(long double aA, long double aB, long double bA, long double bB, long double cA, long double cB, long double qA0, long double qA1, long double qA2, long double qA3, long double qB0, long double qB1, long double qB2, long double qB3, long double rA0, long double rA1, long double rA2, long double rB0, long double rB1, long double rB2, long double vA0, long double vA1, long double vA2, long double vB0, long double vB1, long double vB2, long double omegaA0, long double omegaA1, long double omegaA2, long double omegaB0, long double omegaB1, long double omegaB2, long double lambda0);
    long double fdotAB(long double aA, long double aB, long double bA, long double bB, long double cA, long double cB, long double qA0, long double qA1, long double qA2, long double qA3, long double qB0, long double qB1, long double qB2, long double qB3, long double rA0, long double rA1, long double rA2, long double rB0, long double rB1, long double rB2, long double vA0, long double vA1, long double vA2, long double vB0, long double vB1, long double vB2, long double omegaA0, long double omegaA1, long double omegaA2, long double omegaB0, long double omegaB1, long double omegaB2, long double lambda0);
    long double normdenom(long double q0, long double q1, long double q2, long double q3, long double omega0, long double omega1, long double omega2, long double time);
    long double qtime0(long double q0, long double q1, long double q2, long double q3, long double omega0, long double omega1, long double omega2, long double time);
    long double qtime1(long double q0, long double q1, long double q2, long double q3, long double omega0, long double omega1, long double omega2, long double time);
    long double qtime2(long double q0, long double q1, long double q2, long double q3, long double omega0, long double omega1, long double omega2, long double time);
    long double qtime3(long double q0, long double q1, long double q2, long double q3, long double omega0, long double omega1, long double omega2, long double time);
    long double tobody0( long double q0, long double q1, long double q2, long double q3, long double v0, long double v1, long double v2);
    long double tobody1( long double q0, long double q1, long double q2, long double q3, long double v0, long double v1, long double v2);
    long double tobody2( long double q0, long double q1, long double q2, long double q3, long double v0, long double v1, long double v2); 
    long double tospace0( long double q0, long double q1, long double q2, long double q3, long double v0, long double v1, long double v2);
    long double tospace1( long double q0, long double q1, long double q2, long double q3, long double v0, long double v1, long double v2);
    long double tospace2( long double q0, long double q1, long double q2, long double q3, long double v0, long double v1, long double v2);
#     int collision_time_lambda(long double majorA, long double majorB, long double minorA, long double minorB, long double minor2A, long double minor2B, long double quatAx, long double quatAy, long double quatAz, long double quatAw, long double quatBx, long double quatBy, long double quatBz, long double quatBw, long double posAx, long double posAy, long double posAz, long double posBx, long double posBy, long double posBz, long double velAx, long double velAy, long double velAz, long double velBx, long double velBy, long double velBz, long double angvelAx, long double angvelAy, long double angvelAz, long double angvelBx, long double angvelBy, long double angvelBz, long double *T, long double *lambda0, long double *tc);

def pyf_AB( long double lambda0, long double majorA, long double majorB, long double minorA, long double minorB, long double minor2A, long double minor2B, np.ndarray[np.double_t,ndim=1] qA, np.ndarray[np.double_t,ndim=1] qB, np.ndarray[np.double_t,ndim=1] rA, np.ndarray[np.double_t,ndim=1] rB):
    [qA0, qA1, qA2, qA3] = qA
    [qB0, qB1, qB2, qB3] = qB
    [rA0, rA1, rA2] = rA
    [rB0, rB1, rB2] = rB
    ret = f_AB(<long double> majorA, <long double> majorB, <long double> minorA, <long double> minorB, <long double> minor2A, <long double> minor2B, <long double> qA0, <long double> qA1, <long double> qA2, <long double> qA3, <long double> qB0, <long double> qB1, <long double> qB2, <long double> qB3, <long double> rA0, <long double> rA1, <long double> rA2, <long double> rB0, <long double> rB1, <long double> rB2, <long double> lambda0)
    return ret -1


def pyfAB( long double lambda0, long double majorA, long double majorB, long double minorA, long double minorB, long double minor2A, long double minor2B, np.ndarray[np.double_t,ndim=1] qA, np.ndarray[np.double_t,ndim=1] qB, np.ndarray[np.double_t,ndim=1] rA, np.ndarray[np.double_t,ndim=1] rB):
    [qA0, qA1, qA2, qA3] = qA
    [qB0, qB1, qB2, qB3] = qB
    [rA0, rA1, rA2] = rA
    [rB0, rB1, rB2] = rB
    ret = fAB(<long double> majorA, <long double> majorB, <long double> minorA, <long double> minorB, <long double> minor2A, <long double> minor2B, <long double> qA0, <long double> qA1, <long double> qA2, <long double> qA3, <long double> qB0, <long double> qB1, <long double> qB2, <long double> qB3, <long double> rA0, <long double> rA1, <long double> rA2, <long double> rB0, <long double> rB1, <long double> rB2, <long double> lambda0)
    return ret - 1

# 3D overlap potential derivative WITHOUT denominator (for root solving)
def pyhAB( long double lambda0, long double majorA, long double majorB, long double minorA, long double minorB, long double minor2A, long double minor2B, np.ndarray[np.double_t,ndim=1] qA, np.ndarray[np.double_t,ndim=1] qB, np.ndarray[np.double_t,ndim=1] rA, np.ndarray[np.double_t,ndim=1] rB):
    [qA0, qA1, qA2, qA3] = qA
    [qB0, qB1, qB2, qB3] = qB
    [rA0, rA1, rA2] = rA
    [rB0, rB1, rB2] = rB
    ret = hAB(<long double> majorA, <long double> majorB, <long double> minorA, <long double> minorB, <long double> minor2A, <long double> minor2B, <long double> qA0, <long double> qA1, <long double> qA2, <long double> qA3, <long double> qB0, <long double> qB1, <long double> qB2, <long double> qB3, <long double> rA0, <long double> rA1, <long double> rA2, <long double> rB0, <long double> rB1, <long double> rB2, <long double> lambda0)
    return ret

# 3D overlap potential derivative WITHOUT denominator (for root solving)
def pydhAB( long double lambda0, long double majorA, long double majorB, long double minorA, long double minorB, long double minor2A, long double minor2B, np.ndarray[np.double_t,ndim=1] qA, np.ndarray[np.double_t,ndim=1] qB, np.ndarray[np.double_t,ndim=1] rA, np.ndarray[np.double_t,ndim=1] rB):
    [qA0, qA1, qA2, qA3] = qA
    [qB0, qB1, qB2, qB3] = qB
    [rA0, rA1, rA2] = rA
    [rB0, rB1, rB2] = rB
    ret = dhAB(<long double> majorA, <long double> majorB, <long double> minorA, <long double> minorB, <long double> minor2A, <long double> minor2B, <long double> qA0, <long double> qA1, <long double> qA2, <long double> qA3, <long double> qB0, <long double> qB1, <long double> qB2, <long double> qB3, <long double> rA0, <long double> rA1, <long double> rA2, <long double> rB0, <long double> rB1, <long double> rB2, <long double> lambda0)
    return ret

def pyrC( long double lambda0, long double majorA, long double majorB, long double minorA, long double minorB, long double minor2A, long double minor2B, np.ndarray[np.double_t,ndim=1] qA, np.ndarray[np.double_t,ndim=1] qB, np.ndarray[np.double_t,ndim=1] rA, np.ndarray[np.double_t,ndim=1] rB):
    [qA0, qA1, qA2, qA3] = qA
    [qB0, qB1, qB2, qB3] = qB
    [rA0, rA1, rA2] = rA
    [rB0, rB1, rB2] = rB
    retx = rCx(<long double> majorA, <long double> majorB, <long double> minorA, <long double> minorB, <long double> minor2A, <long double> minor2B, <long double> qA0, <long double> qA1, <long double> qA2, <long double> qA3, <long double> qB0, <long double> qB1, <long double> qB2, <long double> qB3, <long double> rA0, <long double> rA1, <long double> rA2, <long double> rB0, <long double> rB1, <long double> rB2, <long double> lambda0)
    rety = rCy(<long double> majorA, <long double> majorB, <long double> minorA, <long double> minorB, <long double> minor2A, <long double> minor2B, <long double> qA0, <long double> qA1, <long double> qA2, <long double> qA3, <long double> qB0, <long double> qB1, <long double> qB2, <long double> qB3, <long double> rA0, <long double> rA1, <long double> rA2, <long double> rB0, <long double> rB1, <long double> rB2, <long double> lambda0)
    retz = rCz(<long double> majorA, <long double> majorB, <long double> minorA, <long double> minorB, <long double> minor2A, <long double> minor2B, <long double> qA0, <long double> qA1, <long double> qA2, <long double> qA3, <long double> qB0, <long double> qB1, <long double> qB2, <long double> qB3, <long double> rA0, <long double> rA1, <long double> rA2, <long double> rB0, <long double> rB1, <long double> rB2, <long double> lambda0)
    return np.array([retx, rety, retz])

def pynC( long double lambda0, long double majorA, long double majorB, long double minorA, long double minorB, long double minor2A, long double minor2B, np.ndarray[np.double_t,ndim=1] qA, np.ndarray[np.double_t,ndim=1] qB, np.ndarray[np.double_t,ndim=1] rA, np.ndarray[np.double_t,ndim=1] rB):
    [qA0, qA1, qA2, qA3] = qA
    [qB0, qB1, qB2, qB3] = qB
    [rA0, rA1, rA2] = rA
    [rB0, rB1, rB2] = rB
    nCx = normx(<long double> majorA, <long double> majorB, <long double> minorA, <long double> minorB, <long double> minor2A, <long double> minor2B, <long double> qA0, <long double> qA1, <long double> qA2, <long double> qA3, <long double> qB0, <long double> qB1, <long double> qB2, <long double> qB3, <long double> rA0, <long double> rA1, <long double> rA2, <long double> rB0, <long double> rB1, <long double> rB2, <long double> lambda0)
    nCy = normy(<long double> majorA, <long double> majorB, <long double> minorA, <long double> minorB, <long double> minor2A, <long double> minor2B, <long double> qA0, <long double> qA1, <long double> qA2, <long double> qA3, <long double> qB0, <long double> qB1, <long double> qB2, <long double> qB3, <long double> rA0, <long double> rA1, <long double> rA2, <long double> rB0, <long double> rB1, <long double> rB2, <long double> lambda0)
    nCz = normz(<long double> majorA, <long double> majorB, <long double> minorA, <long double> minorB, <long double> minor2A, <long double> minor2B, <long double> qA0, <long double> qA1, <long double> qA2, <long double> qA3, <long double> qB0, <long double> qB1, <long double> qB2, <long double> qB3, <long double> rA0, <long double> rA1, <long double> rA2, <long double> rB0, <long double> rB1, <long double> rB2, <long double> lambda0)
    return np.array([nCx, nCy, nCz])

def pyvC( long double lambda0, long double majorA, long double majorB, long double minorA, long double minorB, long double minor2A, long double minor2B, np.ndarray[np.double_t,ndim=1] qA, np.ndarray[np.double_t,ndim=1] qB, np.ndarray[np.double_t,ndim=1] rA, np.ndarray[np.double_t,ndim=1] rB, np.ndarray[np.double_t,ndim=1] vA, np.ndarray[np.double_t,ndim=1] vB, np.ndarray[np.double_t,ndim=1] omegaA, np.ndarray[np.double_t,ndim=1] omegaB):
    [qA0, qA1, qA2, qA3] = qA
    [qB0, qB1, qB2, qB3] = qB
    [rA0, rA1, rA2] = rA
    [rB0, rB1, rB2] = rB
    [vA0, vA1, vA2] = vA
    [vB0, vB1, vB2] = vB
    [omegaA0, omegaA1, omegaA2] = omegaA
    [omegaB0, omegaB1, omegaB2] = omegaB
    vx = vCx(<long double> majorA, <long double> majorB, <long double> minorA, <long double> minorB, <long double> minor2A, <long double> minor2B, <long double> qA0, <long double> qA1, <long double> qA2, <long double> qA3, <long double> qB0, <long double> qB1, <long double> qB2, <long double> qB3, <long double> rA0, <long double> rA1, <long double> rA2, <long double> rB0, <long double> rB1, <long double> rB2, <long double> vA0, <long double> vA1, <long double> vA2, <long double> vB0, <long double> vB1, <long double> vB2, <long double> omegaA0, <long double> omegaA1, <long double> omegaA2, <long double> omegaB0, <long double> omegaB1, <long double> omegaB2, <long double> lambda0)
    vy = vCy(<long double> majorA, <long double> majorB, <long double> minorA, <long double> minorB, <long double> minor2A, <long double> minor2B, <long double> qA0, <long double> qA1, <long double> qA2, <long double> qA3, <long double> qB0, <long double> qB1, <long double> qB2, <long double> qB3, <long double> rA0, <long double> rA1, <long double> rA2, <long double> rB0, <long double> rB1, <long double> rB2, <long double> vA0, <long double> vA1, <long double> vA2, <long double> vB0, <long double> vB1, <long double> vB2, <long double> omegaA0, <long double> omegaA1, <long double> omegaA2, <long double> omegaB0, <long double> omegaB1, <long double> omegaB2, <long double> lambda0)
    vz = vCz(<long double> majorA, <long double> majorB, <long double> minorA, <long double> minorB, <long double> minor2A, <long double> minor2B, <long double> qA0, <long double> qA1, <long double> qA2, <long double> qA3, <long double> qB0, <long double> qB1, <long double> qB2, <long double> qB3, <long double> rA0, <long double> rA1, <long double> rA2, <long double> rB0, <long double> rB1, <long double> rB2, <long double> vA0, <long double> vA1, <long double> vA2, <long double> vB0, <long double> vB1, <long double> vB2, <long double> omegaA0, <long double> omegaA1, <long double> omegaA2, <long double> omegaB0, <long double> omegaB1, <long double> omegaB2, <long double> lambda0)
    return np.array([vx, vy, vz])

# 3D overlap potential time derivative
def pyfdotAB( long double lambda0, long double majorA, long double majorB, long double minorA, long double minorB, long double minor2A, long double minor2B, np.ndarray[np.double_t,ndim=1] qA, np.ndarray[np.double_t,ndim=1] qB, np.ndarray[np.double_t,ndim=1] rA, np.ndarray[np.double_t,ndim=1] rB, np.ndarray[np.double_t,ndim=1] vA, np.ndarray[np.double_t,ndim=1] vB, np.ndarray[np.double_t,ndim=1] omegaA, np.ndarray[np.double_t,ndim=1] omegaB):
    [qA0, qA1, qA2, qA3] = qA
    [qB0, qB1, qB2, qB3] = qB
    [rA0, rA1, rA2] = rA
    [rB0, rB1, rB2] = rB
    [vA0, vA1, vA2] = vA
    [vB0, vB1, vB2] = vB
    [omegaA0, omegaA1, omegaA2] = omegaA
    [omegaB0, omegaB1, omegaB2] = omegaB
    ret = fdotAB(<long double> majorA, <long double> majorB, <long double> minorA, <long double> minorB, <long double> minor2A, <long double> minor2B, <long double> qA0, <long double> qA1, <long double> qA2, <long double> qA3, <long double> qB0, <long double> qB1, <long double> qB2, <long double> qB3, <long double> rA0, <long double> rA1, <long double> rA2, <long double> rB0, <long double> rB1, <long double> rB2, <long double> vA0, <long double> vA1, <long double> vA2, <long double> vB0, <long double> vB1, <long double> vB2, <long double> omegaA0, <long double> omegaA1, <long double> omegaA2, <long double> omegaB0, <long double> omegaB1, <long double> omegaB2, <long double> lambda0)
    return ret

def pyquat(np.ndarray[np.double_t,ndim=1] q, np.ndarray[np.double_t,ndim=1] omega, long double time):
    [q0, q1, q2, q3] = q
    [omega0, omega1, omega2] = omega
    if omega.all() == 0:
        quatx = q0
        quaty = q1
        quatz = q2
        quatw = q3
    else:
        quatx = qtime0(<long double> q0, <long double> q1, <long double> q2, <long double> q3, <long double> omega0, <long double> omega1, <long double> omega2, <long double> time)
        quaty = qtime1(<long double> q0, <long double> q1, <long double> q2, <long double> q3, <long double> omega0, <long double> omega1, <long double> omega2, <long double> time)
        quatz = qtime2(<long double> q0, <long double> q1, <long double> q2, <long double> q3, <long double> omega0, <long double> omega1, <long double> omega2, <long double> time)
        quatw = qtime3(<long double> q0, <long double> q1, <long double> q2, <long double> q3, <long double> omega0, <long double> omega1, <long double> omega2, <long double> time)
        norms = normdenom(<long double> q0, <long double> q1, <long double> q2, <long double> q3, <long double> omega0, <long double> omega1, <long double> omega2, <long double> time)
    
        quatx = quatx # Sophia FIXME
        quaty = quaty
        quatz = quatz
        quatw = quatw

    return np.array([quatx, quaty, quatz, quatw])

def pytobody(np.ndarray[np.double_t,ndim=1] q, np.ndarray[np.double_t,ndim=1] vec):
    [q0, q1, q2, q3] = q
    [v0, v1, v2] = vec
    tobodyx = tobody0(<long double> q0, <long double> q1, <long double> q2, <long double> q3, <long double> v0, <long double> v1, <long double> v2)
    tobodyy = tobody1(<long double> q0, <long double> q1, <long double> q2, <long double> q3, <long double> v0, <long double> v1, <long double> v2)
    tobodyz = tobody2(<long double> q0, <long double> q1, <long double> q2, <long double> q3, <long double> v0, <long double> v1, <long double> v2)
    return np.array([tobodyx, tobodyy, tobodyz])

def pytospace(np.ndarray[np.double_t,ndim=1] q, np.ndarray[np.double_t,ndim=1] vec):
    [q0, q1, q2, q3] = q
    [v0, v1, v2] = vec
    tospacex = tospace0(<long double> q0, <long double> q1, <long double> q2, <long double> q3, <long double> v0, <long double> v1, <long double> v2)
    tospacey = tospace1(<long double> q0, <long double> q1, <long double> q2, <long double> q3, <long double> v0, <long double> v1, <long double> v2)
    tospacez = tospace2(<long double> q0, <long double> q1, <long double> q2, <long double> q3, <long double> v0, <long double> v1, <long double> v2)
    return np.array([tospacex, tospacey, tospacez])

#cdef double* collision_time_c(ellipse1, ellipse2, double T):
#cdef double* collision_time_c_wrapper(ellipse1, ellipse2, double T):
from libc.stdlib cimport malloc, free
# def collision_time_c(ellipse1, ellipse2, long double T):
#     cdef long double *tmax = &T
#     cdef long double *tc      = <long double *>malloc(sizeof(double))
#     cdef long double *lambda0 = <long double *>malloc(sizeof(double))
#     cdef long double tcx
#     cdef long double lambda0x
#     
#     majorA = ellipse1.major
#     minorA = ellipse1.minor
#     minor2A = minorA
#     [quatAx, quatAy, quatAz, quatAw] = ellipse1.angle
#     [posAx, posAy, posAz] = ellipse1.position
#     [velAx, velAy, velAz] = ellipse1.velocity
#     [angvelAx, angvelAy, angvelAz] = ellipse1.angvel
# 
#     majorB = ellipse2.major
#     minorB = ellipse2.minor
#     minor2B = minorB
#     [quatBx, quatBy, quatBz, quatBw] = ellipse2.angle
#     [posBx, posBy, posBz] = ellipse2.position
#     [velBx, velBy, velBz] = ellipse2.velocity
#     [angvelBx, angvelBy, angvelBz] = ellipse2.angvel
# 
#     status = collision_time_lambda(<long double> majorA, <long double> majorB, <long double> minorA, <long double> minorB, <long double> minor2A, <long double> minor2B, <long double> quatAx, <long double> quatAy, <long double> quatAz, <long double> quatAw, <long double> quatBx, <long double> quatBy, <long double> quatBz, <long double> quatBw, <long double> posAx, <long double> posAy, <long double> posAz, <long double> posBx, <long double> posBy, <long double> posBz, <long double> velAx, <long double> velAy, <long double> velAz, <long double> velBx, <long double> velBy, <long double> velBz, <long double> angvelAx, <long double> angvelAy, <long double> angvelAz, <long double> angvelBx, <long double> angvelBy, <long double> angvelBz, <long double*> tmax, <long double*> lambda0, <long double*> tc)
#     tcx = tc[0]
#     lambda0x = lambda0[0] 
#     free(tc)
#     free(lambda0)
#     
#     if status == 0: # collision at time t with lambda0
#         rc = pyrC( lambda0x, majorA, majorB, minorA, minorB, minor2A, minor2B, ellipse1.angle, ellipse2.angle, ellipse1.position + ellipse1.velocity * tcx, ellipse2.position + ellipse2.velocity * tcx)
#         #return (tcx, lambda0x)
#         return (0, [tcx] + rc.tolist())
#     elif status == -1: # no collision until T
#         return None
#     elif status == -2: # overlap
#         return (-2, 0)
#     elif status == -3: # repredict
#         return (-3, tcx)
        
"""
############################## EXPERIMENTAL
ct helper
"""

from cpython cimport bool

from math import isfinite
from collections import deque
from functools import partial

from scipy.optimize import newton, brentq
from scipy.interpolate import splrep, splev, sproot, splder


#cdef double overlap_potential(ellipse1, ellipse2, double lambda0, double time):
cpdef overlap_potential(double time, double lambda0, ellipse1, ellipse2, double dtswell, double aspectratio):

    qnew1 = pyquat(ellipse1.angle,ellipse1.angvel,time)
    qnew2 = pyquat(ellipse2.angle,ellipse2.angvel,time)
    
    majorA = ellipse1.major + dtswell * time
    minorA = majorA/aspectratio
    majorB = ellipse2.major + dtswell * time
    minorB = majorB/aspectratio

    return pyfAB(lambda0, majorA, majorB, minorA, minorB, minorA, minorB, qnew1, qnew2, ellipse1.position + ellipse1.velocity * time, ellipse2.position + ellipse2.velocity * time)

#cdef double overlap_potential_dlambda0(ellipse1, ellipse2, double lambda0, double time):
cpdef double overlap_potential_dlambda0(double lambda0, double time, ellipse1, ellipse2, double dtswell, double aspectratio):

    qnew1 = pyquat(ellipse1.angle,ellipse1.angvel,time)
    qnew2 = pyquat(ellipse2.angle,ellipse2.angvel,time)
    
    majorA = ellipse1.major + dtswell * time
    minorA = majorA/aspectratio
    majorB = ellipse2.major + dtswell * time
    minorB = majorB/aspectratio
     
    return pyhAB(lambda0, majorA, majorB, minorA, minorB, minorA, minorB, qnew1, qnew2, ellipse1.position + ellipse1.velocity * time, ellipse2.position + ellipse2.velocity * time)

cpdef overlap_potential_d2lambda0(double lambda0, double time, ellipse1, ellipse2, double dtswell, double aspectratio):

    qnew1 = pyquat(ellipse1.angle,ellipse1.angvel,time)
    qnew2 = pyquat(ellipse2.angle,ellipse2.angvel,time)
    majorA = ellipse1.major + dtswell * time
    minorA = majorA/aspectratio
    majorB = ellipse2.major + dtswell * time
    minorB = majorB/aspectratio
    
    return pydhAB(lambda0, majorA, majorB, minorA, minorB, minorA, minorB, qnew1, qnew2, ellipse1.position + ellipse1.velocity * time, ellipse2.position + ellipse2.velocity * time)

cpdef double overlap_potential_dt(double lambda0, double time, ellipse1, ellipse2, double dtswell, double aspectratio):
    majorA = ellipse1.major + dtswell * time
    minorA = majorA/aspectratio
    majorB = ellipse2.major + dtswell * time
    minorB = majorB/aspectratio

    qA = pyquat(ellipse1.angle,ellipse1.angvel,time)
    qB = pyquat(ellipse2.angle,ellipse2.angvel,time)
    
    rA = ellipse1.position + ellipse1.velocity * time
    rB = ellipse2.position + ellipse2.velocity * time
    
    vA = ellipse1.velocity
    vB = ellipse2.velocity
    
    omegaA = ellipse1.angvel
    omegaB = ellipse2.angvel
    
    return pyfdotAB(lambda0, majorA, majorB, minorA, minorB, minorA, minorB, qA, qB, rA, rB,  vA, vB, omegaA, omegaB) # FIXME: CHECK!



"""
ctypedef double (*fABdlambda_ptr)(double, double, double, double, double, double, double, np.ndarray, np.ndarray)
cdef double pyfABdlambda_cdef( double lambda0, double majorA, double minorA, double majorB, double minorB, double thetaA, double thetaB, np.ndarray[np.double_t,ndim=1] rA, np.ndarray[np.double_t,ndim=1] rB):
    [rA0, rA1] = rA
    [rB0, rB1] = rB
    ret = fABdlambla(<double> majorA, <double> majorB, <double> minorA, <double> minorB, <double> thetaA, <double> thetaB, <double> rA0, <double> rA1, <double> rB0, <double> rB1, <double> lambda0)
    return ret
"""

ctypedef double (*fABdlambda_ptr)(double, double, double, double, double, double[::1], double[::1], double[::1], double[::1])
cdef double pyfABdlambda_cdef( double lambda0, double majorA, double minorA, double majorB, double minorB, double[::1] qA, double[::1] qB, double[::1] rA, double[::1] rB):
    [qA0, qA1, qA2, qA3] = qA
    [qB0, qB1, qB2, qB3] = qB
    [rA0, rA1, rA2] = rA
    [rB0, rB1, rB2] = rB
    ret = hAB(<double> majorA, <double> majorB, <double> minorA, <double> minorB, <double> minorA, <double> minorB, <double> qA0, <double> qA1, <double> qA2, <double> qA3, <double> qB0, <double> qB1, <double> qB2, <double> qB3, <double> rA0, <double> rA1, <double> rA2, <double> rB0, <double> rB1, <double> rB2, <double> lambda0)
    return ret


# def __find_time_of_initial_roots__(ellipse1, ellipse2, double lambda0, double tstart, double tend): # FIXME: maybe unsigned double tstart/tend?
#     """
#     TODO: write doc
#     """
#     x = np.linspace(tstart, tend, 10) # FIXME: adjust step n
#     y = np.array([overlap_potential(time, lambda0, ellipse1, ellipse2) for time in x])
#     #y = np.asarray([pyfAB(lambda0, major1, minor1, major2, minor2, angle1 + angvel1*xx, angle2 + angvel2 * xx, position1 + velocity1 * xx, position2 + velocity2 * xx) for xx in x])
#     ##y = np.asarray(pyfABseries(lambda0, major1, minor1, major2, minor2, angle1, angle2, position1, position2, velocity1, velocity2, angvel1, angvel2, x))
#     
#     try:
#         spl = splrep(x,y)
#         roots = sproot(tck=spl, mest=5) # FIXME: mest is number of expected roots... maybe tune this parameter empirically
#         
#         #roots = deque(roots.tolist())
#         return (roots[::-1]).tolist() # reverse to get descending list => pop() gets lowest first
#     except TypeError:
#         pass
#     except ValueError:
#         pass
#     return None
# 
# 
# 
# def __find_time2__(ellipse1, ellipse2, double t0, double lambda0, double tstart, double tend): # FIXME: maybe unsigned double tstart/tend?
#     """
#     Find the time where the contact function has a root for fixed lambda0. 
#     
#     NOTE: may rise TypeError if data to interpolate is crappy or ValueError
#     @param lambda0 lambda0 \in [0,1]
#     """
#     #fx = lambda t: pyfAB(lambda0, major1, minor1, major2, minor2, angle1 + angvel1*t, angle2 + angvel2 * t, position1 + velocity1 * t, position2 + velocity2 * t)
# 
#     # FIXME: check if (fABdot at t=root) << 0
#     #print("overlap_potential_dt(root)=",overlap_potential_dt(t0, lambda0, ellipse1, ellipse2))
#     #fx = lambda t: pyfAB(lambda0, major1, minor1, major2, minor2, angle1 + angvel1*t, angle2 + angvel2 * t, position1 + velocity1 * t, position2 + velocity2 * t)
#     fx = lambda t: overlap_potential(t, lambda0, ellipse1, ellipse2)
#     
#     try:
#         tc = newton(func=fx, x0=t0)                    
# 
#         #if tstart <= b2 <= tend and abs(overlap_potential(b2, lambda0, ellipse1, ellipse2)) < 1e-4: # newton safeguard
#         if tstart <= tc <= tend: # newton safeguard
#             return tc
#     except RuntimeError: # root 'roots' of interpolated function seems to be no root of function F
#         pass # FIXME: extract testcases here or implement additional torquato stuff (i.e. schedule collision detection at last known non-convergent point)
#     
#     """
#     absdiff = (tend-tstart)/100
#     #if overlap_potential_dt(root, lambda0, ellipse1, ellipse2) < -0.1:                    
#     try:
#         (b2, r) = brentq(f=fx, a=t0-absdiff, b=t0+absdiff, full_output=True)
#     except ValueError as ve:
#         print(ve.args)
#     else:
#         print("b2=",b2)
#         if r.converged:
#             return b2
#     """
# 
#     return None
# 
# 
# def tc_finder4(ellipse1, ellipse2, double tstart, double tend):
#     cdef unsigned int counter
#     #cdef double epsilon, lambda0, b, b2, tc, delta_t, delta_lambda # some of these variables can get 'None' so we need a double/None type
#     
#     epsilon = 1e-4
#     #lambda0 = ellipse1.major/(ellipse1.major + ellipse2.major) # guess of the initial lambda
#     lambda0 = min(ellipse1.major, ellipse2.major)/(ellipse1.major + ellipse2.major) # guess of the initial lambda
# 
#     counter = 0
#     tc = tstart
#     delta_t = float('inf') # initial value
#     delta_lambda = float('inf') # initial value
# 
#     b = __find_lambda__(ellipse1, ellipse2, tc, lambda0)
#     if b is not None:
#         lambda0 = b
#     
#     maxit = 10
#     
#     init_roots = __find_time_of_initial_roots__(ellipse1, ellipse2, lambda0, tstart, tend)
#     while init_roots:
#         tc = init_roots.pop()
#         #print("root=",tc)
#         
#         #while counter < maxit and (delta_t > epsilon or delta_lambda > epsilon):
#         while counter < maxit and (delta_t > 1e-4 and delta_lambda > 1e-4):
#             counter += 1
#             
#             b = __find_lambda__(ellipse1, ellipse2, tc, lambda0)
#             if b is not None:
#                 delta_lambda = abs(lambda0 - b)
#                 lambda0 = b
#             else:
#                 pass
#                 #return None
#             
#             b2 = __find_time2__(ellipse1, ellipse2, tc, lambda0, tstart, tend)
#             if b2 is not None:
#                 delta_t = abs(tc - b2)
#                 tc = b2
#             else:
#                 #delta_t = float('inf')
#                 #tend = (tstart + tend)/2
#                 #return None # FIXME: or raise future prediction...
#                 pass
#             
#             #print("IT",counter,": l0= ",lambda0, "; tc=", tc)
#             
#             if b is None and b2 is None: # both not changed so next iterations won't change anything as well => no collision for this spl root
#                 # break while loop
#                 #delta_t = float('inf')
#                 #delta_lambda0 = float('inf')
#                 break
#     
#         """
#         #if counter == maxit and (delta_t > (tend-tstart)*epsilon or delta_lambda > epsilon):
#         if counter == maxit and (delta_t > 1e-2 or delta_lambda > 1e-4):
#             if delta_t > epsilon:
#                 print("T NOT CONVERGING ", delta_t)
#             if delta_lambda > epsilon:
#                 print("lambda NOT CONVERGING ", delta_lambda)
#             if abs(overlap_potential(tc, lambda0, ellipse1, ellipse2)) < .1:
#                 print("NOT CONVERGING BUT VALID SOLUTION, LOL :D")
#             
#             if tstart <= tc <= tend:
#                 raise RuntimeError(tc)
#         
#         return (tc, lambda0, counter, -1)
#         """
#     
#         if abs(overlap_potential(tc, lambda0, ellipse1, ellipse2)) < .1:
#             return (tc, lambda0, counter, -1)
#         elif delta_t == float('inf'): # FIXME: check this
#             pass
#         elif counter == maxit and (delta_t > 1e-4 or delta_lambda > 1e-4):
#     
#             if delta_t > epsilon:
#                 print("T NOT CONVERGING ", delta_t)
#             if delta_lambda > epsilon:
#                 print("lambda NOT CONVERGING ", delta_lambda)
#             if abs(overlap_potential(tc, lambda0, ellipse1, ellipse2)) < .1:
#                 print("NOT CONVERGING BUT VALID SOLUTION, LOL :D")
#     
#             if tstart <= tc <= tend:
#                 raise RuntimeError(tc)
#                 #return (tc, lambda0, counter, -1)
#             else:
#                 #raise RuntimeError(tstart)
#                 #return None
#                 pass
#         else:
#             #return None
#             pass
#     return None




"""
#
# tc_finder3
#
"""

def __find_lambda__(ellipse1, ellipse2, double time, double lambda0, double dtswell, double aspectratio):
    """
    Checked - Aiyin
    bla... FIXME: doc # Simplified version NOT OPTIMIZED
    
    NOTE: this functions may throw a RuntimeError if the newton solver does not converge!
    """
    cdef double b, fdl0, fdl1

    fdl0 = overlap_potential_dlambda0(0, time, ellipse1, ellipse2, dtswell, aspectratio)
    fdl1 = overlap_potential_dlambda0(1, time, ellipse1, ellipse2, dtswell, aspectratio)
    if np.sign(fdl0) == np.sign(fdl1):
        return None
 
    try:
        b = newton(func=overlap_potential_dlambda0, x0=lambda0, args=(time, ellipse1, ellipse2, dtswell, aspectratio)) # secant method
        if 0 < b < 1:
            return b
    except RuntimeError: # should not happen as we have a concave function and checked signs above, but just leave it in case of numerical failure
        pass
    return None



def __find_time__(ellipse1, ellipse2, double lambda0, double tstart, double tend, double dtswell, double aspectratio): # FIXME: maybe unsigned double tstart/tend?
    """
    -checked Aiyin
    Find the time where the contact function has a root for fixed lambda0. 
    
    NOTE: may rise TypeError if data to interpolate is crappy or ValueError
    @param lambda0 lambda0 \in [0,1]
    """

    subdiv = 10  # FIXME: adjust step ########### IMPORTANT
    x = np.linspace(tstart, tend, subdiv)
#     if (ellipse1._id,ellipse2._id) == (23,26):
#         print("#3 checking overlap potentials")
    y = np.array([overlap_potential(time, lambda0, ellipse1, ellipse2, dtswell, aspectratio) for time in x])
#     if (ellipse1._id,ellipse2._id) == (23,26):
#         for time in x:
#             print(time,overlap_potential(time, lambda0, ellipse1, ellipse2))
    try:
        spl = splrep(x,y)
        roots = sproot(tck=spl, mest=5) # FIXME: mest is number of expected roots... maybe tune this parameter empirically
        dspl = splder(spl)
        roots = deque(roots.tolist())
        
        
#         if (ellipse1._id,ellipse2._id) == (23,26):
#             print("#3 lambda0",lambda0,"Fstart",overlap_potential(tstart, lambda0, ellipse1, ellipse2),"Fend",overlap_potential(tend, lambda0, ellipse1, ellipse2))
#             for time in x:
#                 print("> time",time)
#                 print("> position",ellipse1.position + ellipse1.velocity*time,ellipse2.position + ellipse2.velocity*time)
#                 print("> angle",pyquat(ellipse1.angle,ellipse1.angvel,time) ,pyquat(ellipse2.angle,ellipse2.angvel,time))
#                 print("> F0",overlap_potential(time, lambda0, ellipse1, ellipse2))
                
        if not roots: # no root found... no collision in current interval with provided lambda0
#             try:
#                 # maybe insert C algorithm here to find different roots
#             except RuntimeError:
#                 return None
            return None
        else:
            ######### IMPORTANT
            absdiff = (tend-tstart)/(2 * subdiv) # FIXME: zuerst /100, dann /10 ### => dafür brauchen wir eine richtige Abschätzung!!
            while roots:
                root = roots.popleft()
                fx = lambda t: overlap_potential(t, lambda0, ellipse1, ellipse2, dtswell, aspectratio)  
                
                try:
#                     if(ellipse1._id,ellipse2._id) == (23,26):
#                         print("calculating brent result with root",root,"and absdiff",absdiff,"root-absdiff",root-absdiff)
#                         print("fx at a",fx(root-absdiff),"fx at b",fx(root+absdiff),"fx at root", fx(root),"fx at tstart",fx(tstart))
                    bracket_start = root - absdiff
                    bracket_end = root + absdiff
                    if bracket_start < 0: # because the potential function cannot go to negative times
                        bracket_start = tstart # you also assume that bracket_start is always positive
                    (b2, r) = brentq(f = fx,a=bracket_start,b=bracket_end, full_output=True)
#                     (b2, r) = brentq(f=fx, a=root-absdiff, b=root+absdiff, full_output=True)
#                     if(ellipse1._id,ellipse2._id) == (23,26):
#                         print("b2,r",b2,r)
#                         print("fx at root",fx(b2),fx(1e-25))

                except ValueError as ve:
#                     if (ellipse1._id,ellipse2._id) == (23,26):
#                         z = newton(func=fx, x0=root)
#                         print("brentq failed, here is newton",z,fx(z))
                    pass
                else:
#                     if(ellipse1._id,ellipse2._id) == (23,26):
#                         print("r is not converging",r)
                    if r.converged:
                        return b2
                
    except TypeError:
        pass
    except ValueError:
        pass
    return None


def tc_finder3(ellipse1, ellipse2, double tstart, double tend, double dtswell, double aspectratio):
    """
    -checked -Aiyin
    """
    cdef unsigned int counter
    #cdef double epsilon, lambda0, b, b2, tc, delta_t, delta_lambda # some of these variables can get 'None' so we need a double/None type
    
#     if (ellipse1._id,ellipse2._id) == (1,5):
#         print("inside tc_finder3 before computation",tstart,tend)

    epsilon = 1e-4
    lambda0 = ellipse1.major/(ellipse1.major + ellipse2.major) # guess of the initial lambda

    counter = 0
    tc = tstart
    
    delta_t = float('inf') # initial value
    delta_lambda = float('inf') # initial value

    maxit = 10 # FIXME: wähle das dynamisch je nach Intervall und Geschwindigkeit
    #while counter < maxit and (delta_t > 1e-4 or  delta_lambda > 1e-4):
    while counter < maxit and (delta_t > 1e-4 and delta_lambda > 1e-5):
        counter += 1
        
        b = __find_lambda__(ellipse1, ellipse2, tc, lambda0, dtswell, aspectratio)
#         print("from __find_lambda b",b,"tc",tc)
        if b is not None:
#             if (ellipse1._id,ellipse2._id) == (1,5):
#                 print(":) lambda0",lambda0)
            if counter == 1 and lambda0 == b:
                delta_lambda = 1
            else:
                delta_lambda = abs(lambda0 - b)
            
            #if delta_lambda > 0.1:
            #    lambda0 += np.sign(b-lambda0) * 0.1
                
            lambda0 = b

        else:
            pass
        
#         if (ellipse1._id,ellipse2._id) == (23,26):
#             print("#2 lambda0",lambda0)

        b2 = __find_time__(ellipse1, ellipse2, lambda0, tstart, tend, dtswell, aspectratio)
        
        if b2 is not None:
            delta_t = abs(tc - b2)
            tc = b2

        else:
            pass
        
#         if (ellipse1._id,ellipse2._id) == (23,26):
#             print("#4","counter",counter,"lambda",lambda0,"tc",tc)
#             print("Fstart",overlap_potential(tstart, lambda0, ellipse1, ellipse2),"Fend",overlap_potential(tend, lambda0, ellipse1, ellipse2))
#             print("b",b,"b2",b2)
            
        if b is None and b2 is None: # both not changed so next iterations won't change anything as well
            # FIXME: return that no collision will happen at least until tstart or tc?
            return None
    


    if delta_t == float('inf'): # FIXME: check this
        return None
    elif counter == maxit and (delta_t > 1e-4 or delta_lambda > 1e-4):

        if tstart <= tc <= tend:
            raise RuntimeError(tc)
            #return (tc, lambda0, counter, -1)
        else:
            #raise RuntimeError(tstart)
            return None
    elif abs(overlap_potential(tc, lambda0, ellipse1, ellipse2, dtswell, aspectratio)) < 1e-4: # or 1e-5
        return (tc, lambda0, counter, -1)
        
    else:
        return None


# def tc_finder(ellipse1, ellipse2, double tstart, double tend):
#     cdef unsigned int counter
#     cdef double epsilon, lambda0, b, b2, tc, delta_t, delta_lambda, root
#     cdef bool found
#     cdef np.ndarray x
#     
#     epsilon = 1e-4
#     lambda0 = ellipse1.major/(ellipse1.major + ellipse2.major) # guess of the initial lambda
# 
#     counter = 0
#     #tc = 0
#     tc = tstart
#     #tc = float('inf')
#     delta_t = float('inf') # initial value
#     delta_lambda = float('inf') # initial value
#     
#     while counter < 100 and delta_t > epsilon and delta_lambda > epsilon:
#         counter += 1
#         """ find lambda0, so that f_AB(\lambda) gets maximal. This lambda0 can be found by determining the first root of h_{AB} """
#         fdl0 = overlap_potential_dlambda0(0, tc, ellipse1, ellipse2)
#         fdl1 = overlap_potential_dlambda0(1, tc, ellipse1, ellipse2)
#         if np.sign(fdl0) == np.sign(fdl1):
#             return None
#         # FIXME: use proper newton method, not secant method!
#         try: 
#             b = newton(func=overlap_potential_dlambda0, x0=lambda0, args=(tc, ellipse1, ellipse2)) # secant method
#             delta_lambda = abs(b - lambda0)
#             lambda0 = b
#             """ 2. now find root of F0 (NOTE: proof of concept implementation => do this nicer and more efficient in the future) """
#             
#             x = np.linspace(tstart, tend, 10) # FIXME: adjust step n
#             y = np.array([overlap_potential(xx, lambda0, ellipse1, ellipse2) for xx in x])
#             try:
#                 spl = splrep(x,y)
#                 #pol = partial(splev, tck=spl)
#                 roots = sproot(tck=spl, mest=5) # FIXME: mest is number of expected roots... maybe tune this parameter empirically
#             except TypeError: # this error arises if you are interpolating crappy data, e.g. if tstart=tend
#                 print("Provided data not suited for interpolation. Check input!")
#             except ValueError:
#                 return None
#             
#             roots = deque(roots.tolist())
#             if roots is None or not roots: # no root found... no collision in current interval with provided lambda0
#                 tend = (tstart + tend)/2
#             else:
#                 found = False
#                 while not found and roots:
#                     root = roots.popleft()
#                     try:
#                         b2 = newton(func=overlap_potential, x0=root, args=(lambda0, ellipse1, ellipse2)) # FIXME: provide fprime for speedup!
# 
#                         if tstart <= b2 <= tend: # newton safeguard
#                             delta_t = abs(b2 - tc)
#                             tc = b2
#                             found = True
#                         else:
#                             pass
#                     except RuntimeError: # root 'roots' of interpolated function seems to be no root of function F
#                         pass # FIXME: extract testcases here or implement additional torquato stuff (i.e. schedule collision detection at last known non-convergent point)
#  
#         except RuntimeError:
#             pass
# 
#     return (tc, lambda0)
