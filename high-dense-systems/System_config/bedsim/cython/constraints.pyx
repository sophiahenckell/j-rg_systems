"""
Cython module for high performance ellipse constraints.
Usually particle.py should be amended but this is a long term task.
This module is a quick workaround.
"""

import cython
import numpy as np
cimport numpy as np

#from libc.math cimport sin, cos

DTYPE = np.float64
ctypedef np.float64_t DTYPE_t


"""
cpdef E(np.ndarray x, np.ndarray position, np.ndarray velocity, double angle, double angvel, double major, double minor):
    cdef double t, th
    cdef np.ndarray r, p
    
    [t, p] = [x[0], x[1:]]
    th = angle + angvel * t
    
    r = position + velocity * t
    pr = p-r

    return np.dot(np.dot(np.transpose(pr), A(th, major, minor)), pr) - 1.0


cpdef gradE(np.ndarray x, np.ndarray position, np.ndarray velocity, double angle, double angvel, double major, double minor):
    cpdef double t, th
    cpdef np.ndarray r, p, DxE
    cdef DTYPE_t DtE
    
    [t, p] = [x[0], x[1:]]
    om = angvel
    th = angle + om * t
    
    r = position + velocity * t
    pr = p-r

    DxE = 2.0 * np.dot(A(th, major, minor), pr)
    DtE = om * np.dot(np.transpose(pr), np.dot(DtA(th, major, minor), pr)) - 2.0 * np.dot(np.transpose(velocity), np.dot(A(th, major, minor), pr))
    return np.concatenate(([DtE], DxE), axis=0)




cpdef A(double angle, double major, double minor):
    cdef double axes_pref = (1.0/(major*major) - 1.0/(minor*minor))
    cdef np.ndarray mat = np.array([[np.cos(2.0 * angle), np.sin(2.0 * angle)], [np.sin(2.0 * angle), -np.cos(2.0 * angle)]], dtype=np.float64)
    #return 0.5 * ((1.0/self.major**2 - 1.0/self.minor**2)*mat + (1.0/self.major**2 + 1.0/self.minor**2)*np.identity(2, float))
    return 0.5 * axes_pref * (mat + np.identity(2, float))

cpdef DtA(angle, double major, double minor):
    cdef double axes_pref = (1.0/(major*major) - 1.0/(minor*minor))
    cdef mat = np.array([[-np.sin(2.0 * angle), np.cos(2.0 * angle)], [np.cos(2.0 * angle), np.sin(2.0 * angle)]])
    #return (1.0/self.major**2 - 1.0/self.minor**2)*mat
    return axes_pref*mat
"""


cdef class EllipseConstraint:
    cdef np.ndarray position, velocity
    cdef double angle, angvel, major, minor, axes_prefm, axes_prefp
    cdef DTYPE_t s2a, c2a

    """
    cpdef cp(self, with_ellipse_constraint, time):
        p1 = self.position + self.velocity * time
        p2 = with_ellipse_constraint.position + with_ellipse_constraint.velocity * time
        delta = self.cell.cellspace.system.boundary.delta_dir(p1, p2)
        cp = p1 + self.major * delta/np.linalg.norm(delta)
        return cp
    """

    
    #cpdef E(self, np.ndarray x):
    #cpdef DTYPE_t E(self, np.ndarray[DTYPE_t, ndim=1] x): # NOTE: ndim represents the array dimensions, i.e. ndim([..]) = 1, ndim([[..],[..]])=2, etc.
    @cython.profile(True)
    cpdef DTYPE_t E(self, np.ndarray x):
        cdef DTYPE_t t, th
        cdef np.ndarray r, p, pr

        [t, p] = [x[0], x[1:]]
        th = self.angle + self.angvel * t
        
        r = self.position + self.velocity * t
        pr = p-r
    
        return np.dot(np.dot(np.transpose(pr), self.A(th)), pr) - 1.0
    
    
    #cpdef gradE(self, np.ndarray x):
    cpdef np.ndarray gradE(self, np.ndarray[DTYPE_t, ndim=1] x):
        cdef DTYPE_t t, th, om
        cdef DTYPE_t DtE
        cdef np.ndarray r, p, DxE, pr
        """cdef np.ndarray[DTYPE_t, ndim=1] r
        cdef np.ndarray[DTYPE_t, ndim=1] p
        cdef np.ndarray[DTYPE_t, ndim=1] pr
        cdef np.ndarray[DTYPE_t, ndim=1] DxE"""
        
        [t, p] = [x[0], x[1:]]
        om = self.angvel
        th = self.angle + om * t
        
        r = self.position + self.velocity * t
        pr = p-r
    
        DxE = 2.0 * np.dot(self.A(th), pr)
        #DtE = om * np.dot(np.transpose(pr), np.dot(self.DtA(th), pr)) - 2.0 * np.dot(np.transpose(self.velocity), np.dot(self.A(th), pr))
        DtE = om * np.dot(pr, np.dot(self.DtA(th), pr)) - 2.0 * np.dot(self.velocity, np.dot(self.A(th), pr))
        return np.concatenate(([DtE], DxE), axis=0)
        #return np.array([DtE, DxE[0], DxE[1]])
 
    
    #@cython.boundscheck(False) # turn of bounds-checking for entire function
    #cpdef np.ndarray[DTYPE_t, ndim=2] A(self, double angle):
    @cython.profile(True)
    cpdef np.ndarray A(self, double angle):
        cdef DTYPE_t sin2a, cos2a, major2r, minor2r
        cdef np.ndarray mat
        #cdef np.ndarray[DTYPE_t, ndim=2] mat
        
        sin2a = np.sin(2.0 * angle)
        cos2a = np.cos(2.0 * angle)
        
        mat = np.array([[cos2a, sin2a], [sin2a, -cos2a]], dtype=np.float64)
        #mat = np.array([[np.cos(2.0 * angle), np.sin(2.0 * angle)], [np.sin(2.0 * angle), -np.cos(2.0 * angle)]], dtype=np.float64)
        
        major2r = 1/(self.major**2)
        minor2r = 1/(self.minor**2)
        
        return 0.5 * ((major2r - minor2r)*mat + (major2r + minor2r)*np.identity(2, float))
        #return 0.5 * ((1.0/self.major**2 - 1.0/self.minor**2)*mat + (1.0/self.major**2 + 1.0/self.minor**2)*np.identity(2, float))
        #return 0.5 * (self.axes_prefm * mat + self.axes_prefp * np.identity(2, float))
    
    
    #cpdef np.ndarray[DTYPE_t, ndim=2] DtA(self, double angle):
    cpdef np.ndarray DtA(self, double angle):
        #cdef np.ndarray[DTYPE_t, ndim=2] mat
        cdef np.ndarray mat
        cdef DTYPE_t sin2a, cos2a
        sin2a = np.sin(2.0 * angle)
        cos2a = np.cos(2.0 * angle)
        mat = np.array([[-sin2a, cos2a], [cos2a, sin2a]])
        #mat = np.array([[-np.sin(2.0 * angle), np.cos(2.0 * angle)], [np.cos(2.0 * angle), np.sin(2.0 * angle)]])
        #return (1.0/self.major**2 - 1.0/self.minor**2)*mat
        return self.axes_prefm * mat


    #def __init__(self, np.ndarray[DTYPE_t, ndim=1] position, np.ndarray[DTYPE_t, ndim=1] velocity, double angle, double angvel, double major, double minor):
    #def __init__(self, np.ndarray[DTYPE_t, ndim=1] position, np.ndarray[DTYPE_t, ndim=1] velocity, double angle, double angvel, double major, double minor):
    def __init__(self, np.ndarray position, np.ndarray velocity, double angle, double angvel, double major, double minor):
        self.position = position
        self.velocity = velocity
        self.angle = angle
        self.angvel = angvel
        self.major = major
        self.minor = minor
        
        #self.axes_pref = (1.0/major**2 - 1.0/minor**2)
        self.axes_prefm = (1.0/(major*major) - 1.0/(minor*minor))
        self.axes_prefp = (1.0/(major*major) + 1.0/(minor*minor))