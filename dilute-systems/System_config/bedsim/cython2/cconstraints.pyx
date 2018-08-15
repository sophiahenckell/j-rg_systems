import numpy as np
cimport numpy as np

cdef extern from "cconstraints.h": 
    int A(double angle, double major, double minor, double *ret);
    double E(double t, double x, double y, double px, double py, double vx, double vy, double angle, double angvel, double major, double minor, double* grad);

def pyA( np.double_t angle, np.double_t major, np.double_t minor):
    # np.ndarray[np.double_t,ndim=1] A,
    """ wrap np arrays to fc( a.data ... ) """
    #assert N <= len(A) == len(B) == len(Z)
    #fcret = fc( N, <double*> A.data, <double*> B.data, <double*> Z.data )
    #    # fcret = fc( N, A.data, B.data, Z.data )  grr char*
    cdef np.ndarray[np.double_t, ndim=2] aret
    aret = np.array([[0,0],[0,0]], dtype=np.float64)
    i = A(<double> angle, <double> major, <double> minor, <double*> aret.data)
    return aret

def pyE(double t, double x, double y, double px, double py, double vx, double vy, double angle, double angvel, double major, double minor, np.ndarray[np.double_t,ndim=1] grad):
    if grad is not None:
        ret = E(<double> t, <double> x, <double> y, <double> px, <double> py, <double> vx, <double> vy, <double> angle, <double> angvel, <double> major, <double> minor, <double*> grad.data)
        return (ret, grad)
    else:
        grad = np.array([0,0,0], dtype=np.float64)
        ret = E(<double> t, <double> x, <double> y, <double> px, <double> py, <double> vx, <double> vy, <double> angle, <double> angvel, <double> major, <double> minor, <double*> grad.data)
        return (ret, None)
