import cython
import numpy as np
cimport numpy as np


cdef class Boundary(object):
    #cdef unwrap(self, x)
    #cdef delta(self, x, y)
    #cdef delta_dir(self, x, y)
    pass

    
cdef class PeriodicBox(Boundary):
    cpdef public object system
    cpdef np.ndarray __box_size
    cpdef np.ndarray __box_center
    cpdef public _box_corners
    cpdef public np.ndarray extent
    cpdef public _bordervertices
    cpdef public _bordervertices_coord
    cpdef public _eq_vertices
  

    #@cython.locals(xmin=double, xmax=double, ymin=double, ymax=double)
    cpdef np.ndarray[np.double_t,ndim=1] unwrap(self, np.ndarray x)
    
    cpdef np.ndarray[np.double_t,ndim=1] delta(self, np.ndarray x, np.ndarray y)
    
    @cython.locals(size=np.ndarray, xmin=double, xmax=double, ymin=double, ymax=double, delta=np.ndarray)
    cpdef np.ndarray[np.double_t,ndim=1] delta_dir(self, np.ndarray[np.double_t,ndim=1] x, np.ndarray[np.double_t,ndim=1] y)
    
    cpdef __boundary_simplex_to_cell(self, boundary_simplex)
    cpdef __boundary_edge_assignment(self)
    cpdef __boundary_diagonal_corner_neighbour_assignment(self)
    cpdef __grid_to_box_corners(self)