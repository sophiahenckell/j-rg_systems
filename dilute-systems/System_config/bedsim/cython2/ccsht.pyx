import numpy as np
cimport numpy as np

cdef extern from "ccsht.h": 
    int project_on(double proj_angle, double angle, double major, double minor, double *tnpoint);
    int get_normal_projection_interval(double t, double proj_angle, double angle, double angvel, double major, double minor, double *position, double *velocity, double *interval);

def py_project_on(double proj_angle, double angle, double major, double minor, np.ndarray[np.double_t,ndim=1] tnpoint):
    ret = project_on(<double> proj_angle, <double> angle, <double> major, <double> minor, <double*> tnpoint.data)
    return tnpoint

def py_get_normal_projection_interval(double t, double proj_angle, double angle, double angvel, double major, double minor, np.ndarray[np.double_t,ndim=1] position, np.ndarray[np.double_t,ndim=1] velocity, np.ndarray[np.double_t,ndim=1] interval):
    ret = get_normal_projection_interval(<double> t, <double> proj_angle, <double> angle, <double> angvel, <double> major, <double> minor, <double*> position.data, <double*> velocity.data, <double*> interval.data)
    return interval
