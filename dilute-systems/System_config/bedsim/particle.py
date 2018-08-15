# -*- coding: utf-8 -*-
"""
Created on 14.01.2015

@author: mahe
"""

import sys
from math import sin, cos, tan, atan, sqrt, pi, isnan, isfinite, fsum, degrees, copysign
from collections import deque
from functools import partial

import numpy as np
from scipy.optimize import minimize, brentq, bisect, newton
from scipy.interpolate import PchipInterpolator, interp1d, splrep, splev, sproot
from numpy.ma.testutils import approx
from sympy.polys.rationaltools import together
from numpy import arccos
from PyQt4.Qt import QNetworkAccessManager

"""
Try to include some optimized high performance modules.
If not present just use the python reference implementations.
NOTE: Especially for the ellipse constraints the speed gain
is approx factor 16 compared to pure n-dimensional numpy python.
"""
try:
    from bedsim.cython2.coverlap import pyfAB, pyfdotAB, pyhAB, pyrC, pynC, pyvC, pydhAB, pyquat, tc_finder3, pytobody, \
        pytospace  # ,  collision_time_c

    USE_CYTHON = True
except ImportError:
    print("CYTHON COULD NOT BE INITIALIZED")
    USE_CYTHON = False


class ParticleOverlap(ValueError):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


class ConvergenceError(ValueError):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


class Particle(object):
    """
    A 'Particle' is an arbitrary geometric object with
    - a vectorial position
    - a vectorial velocity
    - a timestep when the last update of the particle happened (intrinsic property of event-driven simulations)
    """
    _particle_counter = 0

    def npfsum(self, arrays):
        same_coord = np.transpose(arrays)
        coord_sum = [fsum(sc) for sc in same_coord]
        return np.array(coord_sum)

    def update(self):

        localtime = self.cell.cellspace.system.system_properties.localtime()

        #         swellrate = self.cell.cellspace.system.system_properties.swelling_rate
        #         print("PARTICLE UPDATE before", self._id,self.position,localtime)
        deltaT = localtime - self.lastUpdated
        if deltaT != 0:
            newpos = self.position + self.velocity * deltaT  # calculate movement of the particle with velocity 'velocity' during a deltaT time-step
            self.lastUpdated = localtime  # save the time-step when the particle was updated
            self.position = self.cell.cellspace.system.boundary.unwrap(
                newpos)  # consider boundary conditions # FIXME: this should be a cell crossing event!
            self.cumulative_position += self.velocity * deltaT
            # if self.cell.cellspace.system.system_properties.localtime() >= 0.36300000000000027:
            #     sys.exit("inside update particle")
        # print("newpos",self._id,newpos)

        return deltaT  # return this to save the subtraction above if called from child classes

    def mass(self):
        """
        The mass method should return the mass of the particle.
        Calculation and normalization of the mass is left to the user (e.g. constant density, variable density, proportionality to area, ...).
        """
        pass

    def cell_cross(self, to_cell,
                   crossing_point):  # FIXME: add crossing point and update! should be integrated somehow in this manner due to collision time inaccuracies!
        """
        Perform cell crossing for the current particle to to_cell.
        Change self.cell entry and remove particle reference from old cell.
        @param to_cell New cell for the particle.
        """
        # 0. update particle to new position (prevent wrong cell assignment due to numeric errors)
        self.position = crossing_point
        self.lastUpdated = self.cell.cellspace.system.system_properties.localtime()  # save the time-step when the particle was updated
        # 1. remove particle from self.cell's particle list
        #        if self.cell.particles: # This could be a problem if the cell is empty -A
        self.cell.particles = list(filter(lambda x: x is not self, self.cell.particles))

        # 2. set cell reference of the particle to the new cell
        self.cell = to_cell
        # 3. append current particle to to_cell's particle list
        to_cell.particles.append(self)

    """
    Pinning
    """

    def pin(self):  # particle will not move anymore and gets infinite mass
        self.pinned = True
        self.__oldmass = self.mass
        self.velocity = np.array([0, 0, 0], dtype=np.float64)  # FIXME: use [0,0] or 'None'?

        def newmass():
            return float('inf')

        self.mass = newmass

    def unpin(self):
        self.pinned = False
        self.mass = self.__oldmass
        self.velocity = np.array([0, 0, 0], dtype=np.float64)  # in case we decided to use None

    """
    Serialization
    """

    def get_name(self):
        return "%s_%s" % (type(self).__name__, self._id)

    def get_id(self):
        return int(self._id)

    def get_type(self):
        return type(self).__name__

    def set_summary(self, bmf_data):
        self._id = bmf_data.pop('id')
        self.position = bmf_data.pop('position')
        self.cumulative_position = bmf_data.pop('cumulative_position')
        self.velocity = bmf_data.pop('velocity')
        self.pinned = bmf_data.pop('pinned')
        print("pinned: ", self.pinned)
        # Orientable.set_summary(self, bmf_data) # FIXME: is this correct?!

    def get_summary(self):
        a = {}
        # a['name'] = self.get_name()
        a['id'] = self.get_id()
        a['type'] = self.get_type()
        a['position'] = self.position
        a['cumulative_position'] = self.cumulative_position
        a['velocity'] = self.velocity
        a['pinned'] = self.pinned
        return a

    def get_summary_dynamics(self):
        a = {}
        # a['name'] = self.get_name() # this is static, but keep for identification
        a['id'] = self.get_id()
        a['type'] = self.get_type()
        a['position'] = self.position
        a['cumulative_position'] = self.cumulative_position
        a['velocity'] = self.velocity
        return a

    def get_summary_statics(self):
        a = {}
        # a['name'] = self.get_name()
        a['id'] = self.get_id()
        a['type'] = self.get_type()
        a['pinned'] = self.pinned
        return a

    """
    Constructor
    """

    def __init__(self, **kwargs):
        """
        @param position,velocity vector or scalar.
        @param cell Reference to the Cell object which contains the particle.
        """
        position = kwargs.pop('position')
        velocity = kwargs.pop('velocity')
        #        print('hi i am in particle class right now')
        #        print("Velocity is {}".format(velocity))
        self.pinned = (
            kwargs.pop('pinned', 0) == 1)  # if 'pinned' is not set in HDF-file, set particle as not pinned by default
        if self.pinned:
            self.pin()
        self.position = np.array(position, dtype=np.float64)  # coordinates of the particle

        # TODO: think about definition of cumulative position for t=0
        self.cumulative_position = np.array([0, 0, 0], dtype=np.float64)
        self.velocity = np.array(velocity, dtype=np.float64)  # velocity of the particle

        self.lastUpdated = 0.0  # last time the particle was updated
        self.lastCollision = {}
        self._particle_counter += 1
        # self._id = Particle._particle_counter  # FIXME
        self._id = int(kwargs.pop('id',
                                  Particle._particle_counter))  # if key 'id' is not present use internal counter # FIXME: ensure that this won't create duplicate ids

        self.cell = None  # Assign cell after particle is created. system should deal with the assignment.
        # self.cell = kwargs.pop('cell')

        self.particleTime = 0.0

        """ Shortcut to the boundary class instance """
        # self.cell.cellspace.system.boundary = self.cell.cellspace.system.boundary
        """ Define a shortcut for accessing the system_properties (needed here for boundary unwrapping) """
        # self.cell.cellspace.system.system_properties = self.cell.cellspace.system.system_properties

    # FIXME: this does not work in cython
    # def __eq__(self, another): ## FIXME: make sure that _id is unique!
    #    return hasattr(another, '_id') and self._id == another._id

    def __richcmp__(self, another, mode):
        if mode == 2:
            return hasattr(another, '_id') and self._id == another._id
        else:
            raise NotImplementedError()

    def __hash__(self):
        return hash(self._id)


class Orientable(Particle):
    """
    An 'Orientable' is an object of arbitrary dimension which has the property
    of a distinguished spatial orientation.
    """

    def pin(self):
        self.angvel = 0
        Particle.pin(self)

    def update(self):
        deltaT = Particle.update(self)
        omegap = np.array([self.angvel[0], self.angvel[1], self.angvel[2], 0.])
        self.cumulative_angle += omegap * deltaT
        self.angvel = np.array([0.1, 0.1, 0.1])  # Sophia FIXME
        self.angle = pyquat(self.angle, self.angvel, deltaT)

        return deltaT

    """ Serialization """

    def set_summary(self, bmf_data):
        self.angle = bmf_data.pop('angle')
        self.cumulative_angle = bmf_data.pop('cumulative_angle')
        self.angvel = bmf_data.pop('angvel')
        Particle.set_summary(self, bmf_data)

    def get_summary(self):
        a = Particle.get_summary(self)
        a['angle'] = self.angle
        a['cumulative_angle'] = self.cumulative_angle
        a['angvel'] = self.angvel
        return a

    def get_summary_dynamics(self):
        a = Particle.get_summary_dynamics(self)
        a['angle'] = self.angle
        a['cumulative_angle'] = self.cumulative_angle
        a['angvel'] = self.angvel
        return a

    def get_summary_statics(self):
        a = Particle.get_summary_statics(self)
        return a

    """ Constructor """

    def __init__(self, **kwargs):
        """
        @param position,velocity,angle,angvel vector or scalar
        """
        angle = kwargs.pop('angle')
        angvel = kwargs.pop('angvel')
        self.angle = angle  # orientation angles of the particle
        self.cumulative_angle = np.array([0., 0., 0., 0.], dtype=np.float64)
        self.angvel = angvel  # angular velocities of the particle
        Particle.__init__(self, **kwargs)


class Circle(Particle):
    """
    Minimal reference implementation of a circle.
    It it intentionally (for comparison of the algorithms) not derived from the Ellipse object.
    """

    def mass(self):
        return np.float64(1.0)  # FIMXE: pi*radius**2

    def collision(self, with_circle, collision_point):

        self.lastCollision[with_circle] = self.cell.cellspace.system.system_properties.localtime()
        with_circle.lastCollision[self] = self.lastCollision[with_circle]

        sigma = self.radius + with_circle.radius  # should be the same as norm(delta_r)
        delta_r = self.cell.cellspace.system.boundary.delta_dir(with_circle.position, self.position)  # r1-r2
        delta_v = self.velocity - with_circle.velocity  # v1-v2
        delta_rv = np.dot(delta_r, delta_v)

        # self.velocity = self.velocity - (delta_r * delta_rv/sigma**2) * 2/(1/self.mass() + 1/with_circle.mass())
        # with_circle.velocity = with_circle.velocity + (delta_r * delta_rv/sigma**2) * 2/(1/self.mass() + 1/with_circle.mass())

        self.velocity = self.velocity - 2 * delta_rv / (
            (1 + self.mass() / with_circle.mass()) * sigma) * delta_r / sigma
        with_circle.velocity = with_circle.velocity + 2 * delta_rv / (
            (1 + with_circle.mass() / self.mass()) * sigma) * delta_r / sigma

    def overlap(self, with_circle):
        delta = self.cell.cellspace.system.boundary.delta_dir(self.position, with_circle.position)
        # return np.linalg.norm(delta) < self.radius + with_circle.radius
        return abs(np.linalg.norm(delta) - self.radius + with_circle.radius) < 10 ** -8

    # FIXME: remove when tpw gets default in Ellipse class 
    def tpw_collision_time(self, with_circle, T):
        return self.collision_time(with_circle)

    def collision_time(self, with_circle):
        r_term = lambda obj: obj.position
        v_term = lambda obj: obj.velocity
        rsum = self.radius + with_circle.radius

        r = -self.cell.cellspace.system.boundary.delta_dir(r_term(self), r_term(
            with_circle))  #### FIXME: for 2x2 systems this may give the wrong distance vector against direction of movement
        # r2shifted = r_term(self) - r
        v = v_term(self) - v_term(with_circle)
        b = np.dot(r, v)
        rr = np.dot(r, r)
        vv = np.dot(v, v)
        arg = b ** 2 - vv * (rr - rsum ** 2)

        if arg >= 0:
            tmin = (-b - sqrt(arg)) / vv

            if tmin < 0 and b < 0:
                tmin = 0
        else:
            return None

        pcp = lambda obj: obj.cell.cellspace.system.boundary.unwrap(
            r_term(obj) + v_term(obj) * tmin)  # THIS IS NOT THE COLLISION POINT BUT THE POSITION OF PARTICLE 'obj'!!!!!
        delta = lambda obj1, obj2: self.cell.cellspace.system.boundary.delta_dir(pcp(obj1), pcp(obj2))
        cp = lambda obj1, obj2: obj1.cell.cellspace.system.boundary.unwrap(
            pcp(obj1) + obj1.radius * delta(obj1, obj2) / np.linalg.norm(delta(obj1, obj2)))

        if tmin >= 0:
            return np.concatenate(([tmin], cp(self, with_circle)), axis=0)
        else:
            return None

    def set_summary(self, bmf_data):
        self.radius = bmf_data.pop('radius')
        Particle.set_summary(self, bmf_data)

    def get_summary(self):
        a = Particle.get_summary(self)
        a['radius'] = self.radius
        return a

    def get_summary_dynamics(self):
        a = Particle.get_summary_dynamics(self)
        return a

    def get_summary_statics(self):
        a = Particle.get_summary_statics(self)
        a['radius'] = self.radius
        return a

    def __init__(self, **kwargs):
        """
        @param position,velocity,radius vector or scalar
        """
        self.radius = kwargs.pop('radius')
        Particle.__init__(self, **kwargs)


class Ellipse(Orientable):
    """
    checked -A
    An 'Ellipse' is a two dimensional orientable object with
    - scalar semi-major and semi-minor axes.
    - one scalar angle for orientation (Orientable2D)
    - one scalar angular velocity (Orientable2D)
    - vector position and velocity (Particle)
    """

    def mass(self):
        return np.float64(
            1)  # np.pi * self.major * self.minor # FIXME: only valid for constant density # FIXME: proportionality constant...

    def update(self):

        deltaT = Orientable.update(self)

        return deltaT

    def collision(self, withEllipse, collisionPoint):  # FIXME: NOTE withEllipse! Add support for other geometries!
        """
        Calculate and apply the results of a collision of the current 'Ellipse' with the 'withEllipse' at a certain 'collisionPoint'.
        @param withEllipse Reference to another Ellipse as collision target
        @param collisionPoint Coordinates of the collision point
        Now collisionPoint also includes normal velctor at point of collision -A
        Note omegas are in body coordinates
        """
        self.lastCollision[withEllipse] = self.cell.cellspace.system.system_properties.localtime()
        withEllipse.lastCollision[self] = self.lastCollision[withEllipse]

        rc = collisionPoint[0:3]  # point of collision, already unwrapped
        normal = collisionPoint[
                 3:6]  # gives the end point of the normal vector wrt to with the origin at rc, this has to be translated to the origin of the simbox
        lambda0 = collisionPoint[-1]

        normal_min = self.cell.cellspace.system.boundary.delta_dir(np.array([0., 0., 0.]),
                                                                   normal)  # translated vector from rc to the simbox
        e_perpendicular = normal_min / np.linalg.norm(normal_min)

        #         e_perpendicular_body = lambda obj: pytobody(obj.angle,e_perpendicular)

        #         print("e_perpendicular",e_perpendicular,"collision point",rc,"va",np.linalg.norm(self.velocity),np.linalg.norm(withEllipse.velocity))
        rcx = lambda obj: self.cell.cellspace.system.boundary.delta_dir(rc, obj.position)
        #         rcx_body = lambda obj: pytobody(obj.angle,rcx(obj))

        e_parallel = lambda obj: self.cell.cellspace.system.boundary.delta_dir(np.array([0., 0., 0.]),
                                                                               np.cross(rcx(obj), e_perpendicular))
        #         e_parallel_body = lambda obj: self.cell.cellspace.system.boundary.delta_dir(np.array([0.,0.,0.]),np.cross(rcx_body(obj),e_perpendicular_body(obj)))

        #         s1 = self.angle[3] # get last element of the array
        #         s2 = withEllipse.angle[3]
        #         p1 = self.angle[:3] # get first three elements of the array
        #         p2 = withEllipse.angle[:3]
        #
        #         P1 = np.array([[0,-p1[2],p1[1]],[p1[2],0,-p1[0]],[-p1[1],p1[0],0]])
        #         P2 = np.array([[0,-p2[2],p2[1]],[p2[2],0,-p2[0]],[-p2[1],p2[0],0]])
        #
        #         Q1 = 2.* (p1*p1.transpose() - s1*P1 + (s1*s1 - 0.5)*np.identity(3) )
        #         Q2 = 2.* (p2*p2.transpose() - s2*P2 + (s2*s2 - 0.5)*np.identity(3) )

        v_term = lambda obj: obj.velocity + np.cross(rcx(obj), obj.angvel)

        #         Oinv1 = np.array([[1./self.major,0,0],[0,1./self.minor,0],[0,0,1./self.minor2]])
        #         Oinv2 = np.array([[1./withEllipse.major,0,0],[0,1./withEllipse.minor,0],[0,0,1./withEllipse.minor2]])

        #         Gamma1 = Q1.transpose()*(Oinv1*swell_matrix)*Q1
        #         Gamma2 = Q2.transpose()*(Oinv2*swell_matrix)*Q2


        #         v_term = lambda obj: obj.velocity + np.cross(rcx(obj),obj.angvel) + np.matmul(Gamma(obj),-rcx(obj))

        Theta_term = lambda obj: obj.mass() * (
            obj.minor ** 2 + obj.minor ** 2) / 5.  # lambda obj: 2.0 * obj.mass() * obj.major**2 / 5.0 #Ia,Ib

        vn = np.matmul(e_perpendicular.transpose(),
                       v_term(withEllipse) - v_term(self))  # np.dot(e_perpendicular,v_term(withEllipse) - v_term(self))

        denom = (
            1. / self.mass() + 1. / withEllipse.mass() + (np.linalg.norm(e_parallel(self)) ** 2) / Theta_term(self) + (
                np.linalg.norm(e_parallel(withEllipse)) ** 2) / Theta_term(withEllipse))

        del_pab = 2.0 * vn / denom

        # just to apply shortest image distance befor F0 computation
        r = -self.cell.cellspace.system.boundary.delta_dir(self.position, withEllipse.position)
        r2shifted = self.position - r
        oldpos2 = np.copy(withEllipse.position)
        withEllipse.position = r2shifted

        #         print("COLLISION between",self._id,withEllipse._id,"del_pab",del_pab,"F0",self.overlap_potential(withEllipse,lambda0,0.),"lambda",lambda0)

        #         withEllipse.position = oldpos2

        old_veloc_a = self.velocity
        old_veloc_b = withEllipse.velocity
        old_angvel_a = self.angvel
        old_angvel_b = withEllipse.angvel

        del_pab_vector = del_pab * e_perpendicular

        if del_pab < 0.:
            #             print("original",self.angvel,withEllipse.angvel)
            self.velocity = old_veloc_a + del_pab * e_perpendicular / self.mass()
            withEllipse.velocity = old_veloc_b - del_pab * e_perpendicular / withEllipse.mass()
            self.angvel = old_angvel_a - del_pab * e_parallel(self) / Theta_term(self)
            withEllipse.angvel = old_angvel_b + del_pab * e_parallel(withEllipse) / Theta_term(withEllipse)

            #             r = -self.cell.cellspace.system.boundary.delta_dir(self.position, withEllipse.position)
            #             r2shifted = self.position - r
            #             oldpos2 = np.copy(withEllipse.position)
            #             withEllipse.position = r2shifted

            lambda1 = self.major / (withEllipse.major + self.major)
            fz = lambda l: self.overlap_potential_dlambda0(withEllipse, l, 1e-8)
            lambda1 = newton(func=fz, x0=lambda1)

            F0dot = self.overlap_potential(withEllipse, lambda1, 1e-8) - self.overlap_potential(withEllipse, lambda0,
                                                                                                0.)

        # withEllipse.position = oldpos2

        #             if F0dot < 0.:

        #                 print("original",F0dot,pytospace(self.angle,self.angvel),pytospace(withEllipse.angle,withEllipse.angvel))
        #                 self.angvel = old_angvel_a + del_pab * e_parallel_body(self) / Theta_term(self)
        #                 withEllipse.angvel = old_angvel_b - del_pab * e_parallel_body(withEllipse) / Theta_term(withEllipse)
        #             print("F0",self.overlap_potential(withEllipse,lambda0,0.),"F0dot",F0dot)


        #         else:
        #             print("delpab > 0")

        withEllipse.position = oldpos2

    def collision_for_swelling(self, withEllipse,
                               collisionPoint):  # FIXME: NOTE withEllipse! Add support for other geometries!
        """
        Calculate and apply the results of a collision of the current 'Ellipse' with the 'withEllipse' at a certain 'collisionPoint'.
        @param withEllipse Reference to another Ellipse as collision target
        @param collisionPoint Coordinates of the collision point
        Now collisionPoint also includes normal velctor at point of collision -A
        Note omegas are in body coordinates
        """
        self.lastCollision[withEllipse] = self.cell.cellspace.system.system_properties.localtime()
        withEllipse.lastCollision[self] = self.lastCollision[withEllipse]

        rc = collisionPoint[0:3]  # point of collision, already unwrapped
        normal = collisionPoint[
                 3:6]  # gives the end point of the normal vector wrt to with the origin at rc, this has to be translated to the origin of the simbox
        lambda0 = collisionPoint[-1]

        normal_min = self.cell.cellspace.system.boundary.delta_dir(np.array([0., 0., 0.]),
                                                                   normal)  # translated vector from rc to the simbox
        e_perpendicular = normal_min / np.linalg.norm(normal_min)

        #         e_perpendicular_body = lambda obj: pytobody(obj.angle,e_perpendicular)

        #         print("e_perpendicular",e_perpendicular,"collision point",rc,"va",np.linalg.norm(self.velocity),np.linalg.norm(withEllipse.velocity))
        rcx = lambda obj: self.cell.cellspace.system.boundary.delta_dir(rc, obj.position)
        #         rcx_body = lambda obj: pytobody(obj.angle,rcx(obj))

        e_parallel = lambda obj: self.cell.cellspace.system.boundary.delta_dir(np.array([0., 0., 0.]),
                                                                               np.cross(rcx(obj), e_perpendicular))
        #         e_parallel_body = lambda obj: self.cell.cellspace.system.boundary.delta_dir(np.array([0.,0.,0.]),np.cross(rcx_body(obj),e_perpendicular_body(obj)))

        #         s1 = self.angle[3] # get last element of the array
        #         s2 = withEllipse.angle[3]
        #         p1 = self.angle[:3] # get first three elements of the array
        #         p2 = withEllipse.angle[:3]
        #
        #         P1 = np.array([[0,-p1[2],p1[1]],[p1[2],0,-p1[0]],[-p1[1],p1[0],0]])
        #         P2 = np.array([[0,-p2[2],p2[1]],[p2[2],0,-p2[0]],[-p2[1],p2[0],0]])
        #
        #         Q1 = 2.* (p1*p1.transpose() - s1*P1 + (s1*s1 - 0.5)*np.identity(3) )
        #         Q2 = 2.* (p2*p2.transpose() - s2*P2 + (s2*s2 - 0.5)*np.identity(3) )

        swellrate = self.cell.cellspace.system.system_properties.swelling_rate

        if swellrate is not None:

            swell_matrix = np.identity(3) * swellrate
            s = lambda obj: obj.angle[3]
            p = lambda obj: obj.angle[:3]

            P = lambda obj: np.array(
                [[0, -p(obj)[2], p(obj)[1]], [p(obj)[2], 0, -p(obj)[0]], [-p(obj)[1], p(obj)[0], 0]])
            Q = lambda obj: 2. * (
                p(obj) * p(obj).transpose() - s(obj) * P(obj) + (s(obj) * s(obj) - 0.5) * np.identity(3))
            Oinv = lambda obj: np.array([[1. / obj.major, 0, 0], [0, 1. / obj.minor, 0], [0, 0, 1. / obj.minor2]])
            Gamma = lambda obj: Q(obj).transpose() * (Oinv(obj) * swell_matrix) * Q(obj)

            v_term = lambda obj: obj.velocity + np.cross(rcx(obj), obj.angvel) + np.matmul(Gamma(obj), -rcx(obj))

        else:

            v_term = lambda obj: obj.velocity + np.cross(rcx(obj), obj.angvel)

            #         Oinv1 = np.array([[1./self.major,0,0],[0,1./self.minor,0],[0,0,1./self.minor2]])
            #         Oinv2 = np.array([[1./withEllipse.major,0,0],[0,1./withEllipse.minor,0],[0,0,1./withEllipse.minor2]])

            #         Gamma1 = Q1.transpose()*(Oinv1*swell_matrix)*Q1
            #         Gamma2 = Q2.transpose()*(Oinv2*swell_matrix)*Q2


            #         v_term = lambda obj: obj.velocity + np.cross(rcx(obj),obj.angvel) + np.matmul(Gamma(obj),-rcx(obj))

        Theta_term = lambda obj: obj.mass() * (
            obj.minor ** 2 + obj.minor ** 2) / 5.  # lambda obj: 2.0 * obj.mass() * obj.major**2 / 5.0 #Ia,Ib

        vn = np.matmul(e_perpendicular.transpose(),
                       v_term(withEllipse) - v_term(self))  # np.dot(e_perpendicular,v_term(withEllipse) - v_term(self))

        denom = (
            1. / self.mass() + 1. / withEllipse.mass() + (np.linalg.norm(e_parallel(self)) ** 2) / Theta_term(self) + (
                np.linalg.norm(e_parallel(withEllipse)) ** 2) / Theta_term(withEllipse))

        del_pab = 2.0 * vn / denom

        # just to apply shortest image distance befor F0 computation
        r = -self.cell.cellspace.system.boundary.delta_dir(self.position, withEllipse.position)
        r2shifted = self.position - r
        oldpos2 = np.copy(withEllipse.position)
        withEllipse.position = r2shifted

        print("COLLISION between", self._id, withEllipse._id, "del_pab", del_pab, "F0",
              self.overlap_potential(withEllipse, lambda0, 0.), "lambda", lambda0)

        #         withEllipse.position = oldpos2

        old_veloc_a = self.velocity
        old_veloc_b = withEllipse.velocity
        old_angvel_a = self.angvel
        old_angvel_b = withEllipse.angvel

        del_pab_vector = del_pab * e_perpendicular

        if del_pab < 0.:

            #             print("original",self.angvel,withEllipse.angvel)
            self.velocity = old_veloc_a + del_pab * e_perpendicular / self.mass()
            withEllipse.velocity = old_veloc_b - del_pab * e_perpendicular / withEllipse.mass()
            self.angvel = old_angvel_a - del_pab * e_parallel(self) / Theta_term(self)
            withEllipse.angvel = old_angvel_b + del_pab * e_parallel(withEllipse) / Theta_term(withEllipse)

            lambda1 = self.major / (withEllipse.major + self.major)
            fz = lambda l: self.overlap_potential_dlambda0(withEllipse, l, 1e-8)
            lambda1 = newton(func=fz, x0=lambda1)

            F0dot = self.overlap_potential(withEllipse, lambda1, 1e-8) - self.overlap_potential(withEllipse, lambda0,
                                                                                                0.)

            print("F0", self.overlap_potential(withEllipse, lambda0, 0.), "F0dot", F0dot)

        else:
            print("delpab > 0")

        withEllipse.position = oldpos2

    def correction(self, withEllipse):
        print("correction move is being applied", self._id, withEllipse._id)
        """
        Calculate and apply a correction move between the current and 'withEllipse' particles.
        @param withEllipse Correction with this ellipse is desired
        """
        # Set currect time step as 'last collision time'
        self.lastCollision[withEllipse] = self.cell.cellspace.system.system_properties.localtime()
        withEllipse.lastCollision[self] = self.lastCollision[withEllipse]

        r = -self.cell.cellspace.system.boundary.delta_dir(self.position, withEllipse.position)

        v1 = self.velocity
        v2 = withEllipse.velocity
        vsum = self.cell.cellspace.system.boundary.unwrap(v1 + v2)
        vdiff = self.cell.cellspace.system.boundary.delta_dir(v1, v2)
        #         unit_vsum = vsum/np.linalg.norm(vsum)

        # let e_parallel = e1, here you just want to construct some 3 random new euclidean coordinates -Aiyin
        e1 = r / np.linalg.norm(r)  # first unit vector parallel to the center of mass separation

        e2 = np.cross(e1, vsum)

        e2 = e2 / np.linalg.norm(e2)  # unit vector 2
        e3 = np.cross(e1, e2)  # unit vector 3

        if not self.pinned and not withEllipse.pinned:  # both not pinned
            m1 = self.mass()
            m2 = withEllipse.mass()
            self.velocity = 1 / (2 * m1) * (
                np.dot(e1, (np.dot(m1 * v1 + m2 * v2, e1) + np.linalg.norm(m1 * v1 - m2 * v2))) + np.dot(e2, np.dot(
                    m1 * v1 + m2 * v2, e2)) + np.dot(e3, np.dot(m1 * v1 + m2 * v2, e3)))
            withEllipse.velocity = 1 / (2 * m2) * (
                np.dot(e1, (np.dot(m1 * v1 + m2 * v2, e1) - np.linalg.norm(m2 * v1 - m2 * v2))) + np.dot(e2, np.dot(
                    m1 * v1 + m2 * v2, e2)) + np.dot(e3, np.dot(m1 * v1 + m2 * v2, e3)))
        elif self.pinned and not withEllipse.pinned:  # only self pinned
            withEllipse.velocity = 0.5 * (
                np.dot(e_parallel, (np.dot(v2, e_parallel) - np.linalg.norm(v2))) + np.dot(e_perpendicular,
                                                                                           np.dot(v2, e_perpendicular)))
        elif not self.pinned and withEllipse.pinned:  # only withEllipse pinned
            self.velocity = 0.5 * (
                np.dot(e_parallel, (np.dot(v1, e_parallel) + np.linalg.norm(v1))) + np.dot(e_perpendicular,
                                                                                           np.dot(v1, e_perpendicular)))
        else:  # both pinned
            pass

    def overlap(self, withEllipse, dt):
        lambda0 = self.major / (self.major + withEllipse.major)
        fdl0 = self.overlap_potential_dlambda0(withEllipse, 0, dt)
        fdl1 = self.overlap_potential_dlambda0(withEllipse, 1, dt)
        if fdl0 == 0 or fdl1 == 0:
            psi_term = -1  # overlap #ugly
        else:
            try:
                #                 print("position",self.position,withEllipse.position)
                #                 print("velocity",self.velocity,withEllipse.velocity)
                #                 print("angle",self.angle,withEllipse.angle)
                #                 print("angvel",self.angvel,withEllipse.angvel)
                fdlz = lambda l: self.overlap_potential_dlambda0(withEllipse, l, dt)
                b = newton(func=fdlz, x0=lambda0)  # secant method

                if 0 < b < 1:
                    lambda0 = b
                    psi_term = self.overlap_potential(withEllipse, lambda0, dt)

            except RuntimeError:  # should not happen as we have a concave function and checked signs above, but just leave it in case of numerical failure
                print("Error")

        return (psi_term < 0.)

    """
    Ellipse overlap functions according to Torquato, Perram & Wertheim
    """

    def overlap_potential(self, withEllipse, lambda0,
                          time):  # set f in constructor to cython version of python version if available
        raise NotImplementedError()

    def __f_cython(self, withEllipse, lambda0, time):
        if isfinite(time):
            qnewA = pyquat(self.angle, self.angvel, time)
            qnewB = pyquat(withEllipse.angle, withEllipse.angvel, time)
            #             if (self._id,withEllipse._id) == (62,126):
            #                 print("time",time,np.linalg.norm(self.position - withEllipse.position))
            #                 print("position",self.position + self.velocity * time,withEllipse.position + withEllipse.velocity*time)
            #                 print("angle new",qnewA,qnewB)
            #                 print("angvel",self.angvel,withEllipse.angvel)
            #                 print("F0",pyfAB(lambda0, self.major, withEllipse.major, self.minor, withEllipse.minor, self.minor2, withEllipse.minor2, qnewA, qnewB, self.position + self.velocity * time, withEllipse.position + withEllipse.velocity * time))
            swellrate = self.cell.cellspace.system.system_properties.swelling_rate
            aspectratio = self.cell.cellspace.system.system_properties.aspect_ratio
            if swellrate is not None:
                majorA = self.major + swellrate * time
                majorB = withEllipse.major + swellrate * time
                minorA = majorA / aspectratio
                minorB = majorB / aspectratio
            else:
                majorA = self.major
                majorB = withEllipse.major
                minorA = self.minor
                minorB = withEllipse.minor

            return pyfAB(lambda0, majorA, majorB, minorA, minorB, minorA, minorB, qnewA, qnewB,
                         self.position + self.velocity * time, withEllipse.position + withEllipse.velocity * time)
        else:
            return float('inf')  # FIXME

    # derivative of overlap_potential by time
    def overlap_potential_dt(self, withEllipse, lambda0,
                             time):  # set h in constructor to cython version of python version if available
        raise NotImplementedError()

    def __fdot_cython(self, withEllipse, lambda0, time):
        qnewA = pyquat(self.angle, self.angvel, time)
        qnewB = pyquat(withEllipse.angle, withEllipse.angvel, time)
        swellrate = self.cell.cellspace.system.system_properties.swelling_rate
        aspectratio = self.cell.cellspace.system.system_properties.aspect_ratio
        if swellrate is not None:
            majorA = self.major + swellrate * time
            majorB = withEllipse.major + swellrate * time
            minorA = majorA / aspectratio
            minorB = majorB / aspectratio
        else:
            majorA = self.major
            majorB = withEllipse.major
            minorA = self.minor
            minorB = withEllipse.minor

        return pyfdotAB(lambda0, majorA, majorB, minorA, minorB, minorA, minorB, qnewA, qnewB,
                        self.position + self.velocity * time, withEllipse.position + withEllipse.velocity * time,
                        self.velocity, withEllipse.velocity, self.angvel, withEllipse.angvel)

    # derivative of overlap_potential by lambda0
    def overlap_potential_dlambda0(self, withEllipse, lambda0,
                                   time):  # set h in constructor to cython version of python version if available
        raise NotImplementedError()

    #     def __h_cython(self, withEllipse, lambda0, time):
    #         qnewA = pyquat(self.angle,self.angvel,time)
    #         qnewB = pyquat(withEllipse.angle, withEllipse.angvel, time)
    #         return pyhAB(lambda0, self.major, withEllipse.major, self.minor, withEllipse.minor, self.minor2, withEllipse.minor2, qnewA, qnewB, self.position + self.velocity * time, withEllipse.position + withEllipse.velocity * time)

    def __fABdlambda_cython(self, withEllipse, lambda0, time):
        qnewA = pyquat(self.angle, self.angvel, time)
        qnewB = pyquat(withEllipse.angle, withEllipse.angvel, time)
        swellrate = self.cell.cellspace.system.system_properties.swelling_rate
        aspectratio = self.cell.cellspace.system.system_properties.aspect_ratio
        if swellrate is not None:
            majorA = self.major + swellrate * time
            majorB = withEllipse.major + swellrate * time
            minorA = majorA / aspectratio
            minorB = majorB / aspectratio
        else:
            majorA = self.major
            majorB = withEllipse.major
            minorA = self.minor
            minorB = withEllipse.minor

        return pyhAB(lambda0, majorA, majorB, minorA, minorB, minorA, minorB, qnewA, qnewB,
                     self.position + self.velocity * time, withEllipse.position + withEllipse.velocity * time)

    def overlap_potential_d2lambda0(self, withEllipse, lambda0,
                                    time):  # set h in constructor to cython version of python version if available
        raise NotImplementedError()

    def __fABd2lambda_cython(self, withEllipse, lambda0, time):
        qnewA = pyquat(self.angle, self.angvel, time)
        qnewB = pyquat(withEllipse.angle, withEllipse.angvel, time)
        swellrate = self.cell.cellspace.system.system_properties.swelling_rate
        aspectratio = self.cell.cellspace.system.system_properties.aspect_ratio
        if swellrate is not None:
            majorA = self.major + swellrate * time
            majorB = withEllipse.major + swellrate * time
            minorA = majorA / aspectratio
            minorB = majorB / aspectratio
        else:
            majorA = self.major
            majorB = withEllipse.major
            minorA = self.minor
            minorB = withEllipse.minor

        return pydhAB(lambda0, majorA, majorB, minorA, minorB, minorA, minorB, qnewA, qnewB,
                      self.position + self.velocity * time, withEllipse.position + withEllipse.velocity * time)

    def overlap_collision_point(self, withEllipse, lambda0, time):
        raise NotImplementedError()

    def __rC_cython(self, withEllipse, lambda0, time):
        qnewA = pyquat(self.angle, self.angvel, time)
        qnewB = pyquat(withEllipse.angle, withEllipse.angvel, time)
        swellrate = self.cell.cellspace.system.system_properties.swelling_rate
        aspectratio = self.cell.cellspace.system.system_properties.aspect_ratio
        if swellrate is not None:
            majorA = self.major + swellrate * time
            majorB = withEllipse.major + swellrate * time
            minorA = majorA / aspectratio
            minorB = majorB / aspectratio
        else:
            majorA = self.major
            majorB = withEllipse.major
            minorA = self.minor
            minorB = withEllipse.minor
        return pyrC(lambda0, majorA, majorB, minorA, minorB, minorA, minorB, qnewA, qnewB,
                    self.position + self.velocity * time, withEllipse.position + withEllipse.velocity * time)

    def __nC_cython(self, withEllipse, lambda0, time):
        qnewA = pyquat(self.angle, self.angvel, time)
        qnewB = pyquat(withEllipse.angle, withEllipse.angvel, time)
        swellrate = self.cell.cellspace.system.system_properties.swelling_rate
        aspectratio = self.cell.cellspace.system.system_properties.aspect_ratio
        if swellrate is not None:
            majorA = self.major + swellrate * time
            majorB = withEllipse.major + swellrate * time
            minorA = majorA / aspectratio
            minorB = majorB / aspectratio
        else:
            majorA = self.major
            majorB = withEllipse.major
            minorA = self.minor
            minorB = withEllipse.minor

        return pynC(lambda0, majorA, majorB, minorA, minorB, minorA, minorB, qnewA, qnewB,
                    self.position + self.velocity * time, withEllipse.position + withEllipse.velocity * time)

    def __vC_cython(self, withEllipse, lambda0, time):
        qnewA = pyquat(self.angle, self.angvel, time)
        qnewB = pyquat(withEllipse.angle, withEllipse.angvel, time)
        swellrate = self.cell.cellspace.system.system_properties.swelling_rate
        aspectratio = self.cell.cellspace.system.system_properties.aspect_ratio
        if swellrate is not None:
            majorA = self.major + swellrate * time
            majorB = withEllipse.major + swellrate * time
            minorA = majorA / aspectratio
            minorB = majorB / aspectratio
        else:
            majorA = self.major
            majorB = withEllipse.major
            minorA = self.minor
            minorB = withEllipse.minor

        return pyvC(lambda0, majorA, majorB, minorA, minorB, minorA, minorB, qnewA, qnewB,
                    self.position + self.velocity * time, withEllipse.position + withEllipse.velocity * time)

    def tpw_collision_time3(self, withEllipse, T=None):

        """
        checked -Aiyin but self.E not included
        Collision time algorithm by Torquato, Perram and Wertheim
        @param withEllipse check collision with 'withEllipse'
        @param T check for collisions in the time interval [0,T]
        Returns collision time and collision point (and maybe normal vector is useful too) -A
        """
        # if self.overlap(withEllipse, 0) and self.overlap(withEllipse, 1e-10):
        #    raise ParticleOverlap("%s and %s are overlapping." % (self.get_name(), withEllipse.get_name()))
        ct_approx = self.circle_approximation(withEllipse)
        if ct_approx is not None:
            (tmin, tmax) = ct_approx
        else:
            return None

        # steps (a) and (b)
        if T is not None:
            if tmin > T:
                return None
            else:
                tend = min(T, tmax)
        else:
            tend = tmax
        tstart = tmin

        r = -self.cell.cellspace.system.boundary.delta_dir(self.position, withEllipse.position)
        r2shifted = self.position - r
        oldpos2 = np.copy(withEllipse.position)
        #         if not np.array_equal(withEllipse.position, r2shifted):
        #             print("check position shifting",withEllipse.position,r2shifted)
        withEllipse.position = r2shifted
        #         print("after",self._id,withEllipse._id,self.position,withEllipse.position)

        # step (c)
        lambda0 = self.major / (self.major + withEllipse.major)  # guess of the initial lambda

        # """ # FIXME: should be activated...
        # initial lambda0 may be the max. of F at time tstart, not the lambda0 initial guess above
        fdl0 = self.overlap_potential_dlambda0(withEllipse, 0, tstart)
        fdl1 = self.overlap_potential_dlambda0(withEllipse, 1, tstart)

        if np.sign(fdl0) == np.sign(fdl1):
            withEllipse.position = oldpos2
            return None
        try:
            fz = lambda l: self.overlap_potential_dlambda0(withEllipse, l, tstart)
            lambda0 = newton(func=fz, x0=lambda0)  # secant method
        # print(self._id,withEllipse._id,lambda0)
        except RuntimeError:
            pass
            # NOTE: if you change lambda0 here, i.e. if you calculate if for another t you'll have to change the overlap checking time-point as well
            # """
            # DEBUG CODE: check if cpt lies on boundary of both ellipses
            # print("collision point lies on boundary of ellipses: ", self.E(cpt), "\tand ", withEllipse.E(cpt), ";\toverlap: ",self.overlap(withEllipse, 0))

        #         #OVERRIDE EVERYTHING
        #         tstart = 0
        #         tend = 0.231481481481

        F0 = self.overlap_potential(withEllipse, lambda0, tstart)
        F0dot2 = self.overlap_potential_dt(withEllipse, lambda0, tstart)

        # catch circles or exact coll.
        if tstart == tend:  # FIXME: maybe move above, maybe simplify, maybe maybe... :D
            cp1 = self.overlap_collision_point(withEllipse, lambda0, tstart)
            cp1_unwrap = self.cell.cellspace.system.boundary.unwrap(cp1)
            cp2 = self.overlap_normal_vector(withEllipse, lambda0, tstart)
            cp2_unwrap = self.cell.cellspace.system.boundary.unwrap(cp2)
            cp3 = np.append(cp1_unwrap, cp2_unwrap)
            cp = np.append(cp3, lambda0)
            cpt = np.concatenate(([tstart], cp), axis=0)
            withEllipse.position = oldpos2
            #             print("calc 1,tstart = tend",self._id,withEllipse._id)
            return cpt

        if tstart < 1e-8 and F0 < 0 and self.overlap_potential(withEllipse, lambda0,
                                                               tstart + 1e-8) < 0:  # in case lambda0 is for t=tstart

            withEllipse.position = oldpos2
            # print("OVERLAP, tmin=",tmin," tmax= ",tmax)
            cp1 = self.overlap_collision_point(withEllipse, lambda0, tstart)
            cp1_unwrap = self.cell.cellspace.system.boundary.unwrap(cp1)
            cp2 = self.overlap_normal_vector(withEllipse, lambda0, tstart)
            cp2_unwrap = self.cell.cellspace.system.boundary.unwrap(cp2)
            cp3 = np.append(cp1_unwrap, cp2_unwrap)
            cp = np.append(cp3, lambda0)
            cpt = np.concatenate(([tstart], cp), axis=0)
            #             print("calc2, OVERLAP!",self._id,withEllipse._id,"tc",tstart,"F0",F0)
            raise ParticleOverlap("%s and %s are overlapping." % (self.get_name(), withEllipse.get_name()))

            # """


            #         epsilon = 1e-4 # 1e-3 - 1e-4 according to Torquato
            #         # step (d)
            #         # FIXME: secant method is used atm, because I suspect that the Torquato time derivative is not correct... at least my calculation of the torquato derivative
            #         F0dot = (self.overlap_potential(withEllipse, lambda0, tstart+1e-4) - F0)/1e-4 # I don't trust the fdotAB implementation
            #         #""" # FIXME: should be activated
            #
            #         if abs(F0) < epsilon: # ellipses are overlapping or nearly touching
            #             #FIXME: Fdot not working properly!!
            #             #print("F0=",F0, "; F0dot=",F0dot)
            #             if F0dot >= 0: # ellipse distance increasing
            #                 F0 = epsilon
            # #                 print("calc3",self._id,withEllipse._id,"ellipse distance increasing")
            #             else: # ellipses are approaching
            #                 cp1 = self.overlap_collision_point(withEllipse, lambda0, tstart)
            #                 cp1_unwrap = self.cell.cellspace.system.boundary.unwrap(cp1)
            #                 cp2 = self.overlap_normal_vector(withEllipse,lambda0,tstart)
            #                 cp2_unwrap = self.cell.cellspace.system.boundary.unwrap(cp2)
            #                 cp3 = np.append(cp1_unwrap,cp2_unwrap)
            #                 cp = np.append(cp3,lambda0)
            # #                 rc = self.position + self.velocity*tstart + (1.0-lambda0)*np.matmul(cp1[0],cp2_unwrap)
            # #                 rc_unwrap = self.cell.cellspace.system.boundary.unwrap(rc)
            # #                 cp = np.append(rc_unwrap,cp2_unwrap)
            #                 cpt = np.concatenate(([tstart],cp), axis=0)
            #                 #print("EXIT2: cpoint lies on boundary of ellipses: ", self.E(cpt), "\tand ", withEllipse.E(cpt), ";\toverlap: ",self.overlap(withEllipse, 0))
            #                 withEllipse.position = oldpos2
            # #                 print("calc4",self._id,withEllipse._id,"tc",tstart,"F0",F0,"F0dot",F0dot)
            #
            #                 return cpt
            # #         print(self._id,withEllipse._id,F0,F0dot,F0dot2,"going inside tc_finder3")
            #         try:
            # #             if (self._id,withEllipse._id) == (23,26):
            # #                 print("#1 going to tc_finder3",self._id,withEllipse._id,"F0",F0,"tstart",tstart,"tend",tend)
            #             aspectratio = self.cell.cellspace.system.system_properties.aspect_ratio
            #             tcf = tc_finder3(ellipse1=self, ellipse2=withEllipse, tstart=tstart, tend=tend, dtswell = 0.,aspectratio = aspectratio)
            #
            #         except RuntimeError as runerr:
            #             val = runerr.args
            #             if val is not None:
            #                 withEllipse.position = oldpos2
            #                 raise ConvergenceError(val)
            #             else:
            #                 withEllipse.position = oldpos2
            #                 return None
            #         else:
            #             if tcf is not None:
            #                 (tc, lambda0, series, rte) = tcf
            #
            #                 if isfinite(tc) and tc >= 0 and self.overlap_potential(withEllipse, lambda0, tc) < epsilon:
            # #                     cp = self.overlap_collision_point(withEllipse, lambda0, tc)
            #                     cp1 = self.overlap_collision_point(withEllipse, lambda0, tc)
            #                     cp1_unwrap = self.cell.cellspace.system.boundary.unwrap(cp1)
            #                     cp2 = self.overlap_normal_vector(withEllipse,lambda0,tc)
            #                     cp2_unwrap = self.cell.cellspace.system.boundary.unwrap(cp2)
            # #                     rc = self.position + self.velocity*tc + (1.0-lambda0)*np.matmul(cp1[0],cp2_unwrap)
            # #                     rc_unwrap = self.cell.cellspace.system.boundary.unwrap(rc)
            #                     cp3 = np.append(cp1_unwrap,cp2_unwrap)
            #                     cp = np.append(cp3,lambda0)
            # #                     cp = np.append(cp1_unwrap,cp2_unwrap)
            #                     cpt = np.concatenate(([tc],cp), axis=0)
            # #                     print("calc5",self.major,self._id,withEllipse._id,"tc",tc,"F0dot",F0dot)#"overlap potential",self.overlap_potential(withEllipse,lambda0,tc))
            #
            #                     withEllipse.position = oldpos2
            #                     return cpt
            #
            #         withEllipse.position = oldpos2
            #         return None # FIXME: return $t_c^tilde$ instead
            #

    def find_swelling_time(self, withEllipse, T=None):
        """
        checked -Aiyin but self.E not included
 
        @param withEllipse check collision with 'withEllipse'
        @param T check for collisions in the time interval [0,T]
        Returns collision time and collision point (and maybe normal vector is useful too) -A
        """
        # if self.overlap(withEllipse, 0) and self.overlap(withEllipse, 1e-10):
        #    raise ParticleOverlap("%s and %s are overlapping." % (self.get_name(), withEllipse.get_name()))
        ct_approx = self.circle_approximation(withEllipse)

        if ct_approx is not None:
            (tmin, tmax) = ct_approx
        else:
            return None

        # steps (a) and (b)
        if T is not None:
            if tmin > T:
                return None
            else:
                tend = min(T, tmax)
        else:
            tend = tmax
        tstart = tmin
        # Move ellipse to avoid discontinuities,
        # NO NEED FOR THIS
        #         print("before",self._id,withEllipse._id,self.position,withEllipse.position)

        r = -self.cell.cellspace.system.boundary.delta_dir(self.position, withEllipse.position)
        r2shifted = self.position - r
        oldpos2 = np.copy(withEllipse.position)
        #         if not np.array_equal(withEllipse.position, r2shifted):
        #             print("check position shifting",withEllipse.position,r2shifted)
        withEllipse.position = r2shifted
        #         print("after",self._id,withEllipse._id,self.position,withEllipse.position)

        # step (c)
        lambda0 = self.major / (self.major + withEllipse.major)  # guess of the initial lambda

        # """ # FIXME: should be activated...
        # initial lambda0 may be the max. of F at time tstart, not the lambda0 initial guess above
        fdl0 = self.overlap_potential_dlambda0(withEllipse, 0, tstart)
        fdl1 = self.overlap_potential_dlambda0(withEllipse, 1, tstart)

        if np.sign(fdl0) == np.sign(fdl1):
            withEllipse.position = oldpos2
            return None
        try:
            fz = lambda l: self.overlap_potential_dlambda0(withEllipse, l, tstart)
            lambda0 = newton(func=fz, x0=lambda0)  # secant method
        # print(self._id,withEllipse._id,lambda0)
        except RuntimeError:
            pass
            # NOTE: if you change lambda0 here, i.e. if you calculate if for another t you'll have to change the overlap checking time-point as well
            # """
            # DEBUG CODE: check if cpt lies on boundary of both ellipses
            # print("collision point lies on boundary of ellipses: ", self.E(cpt), "\tand ", withEllipse.E(cpt), ";\toverlap: ",self.overlap(withEllipse, 0))

        #         #OVERRIDE EVERYTHING
        #         tstart = 0
        #         tend = 0.231481481481
        F0 = self.overlap_potential(withEllipse, lambda0, tstart)
        F0dot2 = self.overlap_potential_dt(withEllipse, lambda0, tstart)

        # catch circles or exact coll.
        if tstart == tend:  # FIXME: maybe move above, maybe simplify, maybe maybe... :D
            cp1 = self.overlap_collision_point(withEllipse, lambda0, tstart)
            cp1_unwrap = self.cell.cellspace.system.boundary.unwrap(cp1)
            cp2 = self.overlap_normal_vector(withEllipse, lambda0, tstart)
            cp2_unwrap = self.cell.cellspace.system.boundary.unwrap(cp2)
            cp3 = np.append(cp1_unwrap, cp2_unwrap)
            cp = np.append(cp3, lambda0)
            cpt = np.concatenate(([tstart], cp), axis=0)
            withEllipse.position = oldpos2
            #             print("swell calc 1,tstart = tend",self._id,withEllipse._id)
            return cpt

        if tstart < 1e-8 and F0 < 0 and self.overlap_potential(withEllipse, lambda0,
                                                               tstart + 1e-8) < 0:  # in case lambda0 is for t=tstart

            withEllipse.position = oldpos2
            # print("OVERLAP, tmin=",tmin," tmax= ",tmax)
            cp1 = self.overlap_collision_point(withEllipse, lambda0, tstart)
            cp1_unwrap = self.cell.cellspace.system.boundary.unwrap(cp1)
            cp2 = self.overlap_normal_vector(withEllipse, lambda0, tstart)
            cp2_unwrap = self.cell.cellspace.system.boundary.unwrap(cp2)
            cp3 = np.append(cp1_unwrap, cp2_unwrap)
            cp = np.append(cp3, lambda0)
            cpt = np.concatenate(([tstart], cp), axis=0)
            print("swell calc2, OVERLAP!", self._id, withEllipse._id, "major", self.major, self.minor, "tc", tstart,
                  "F0", F0)
            #             sys.exit("Overlap detected")
            raise ParticleOverlap("%s and %s are overlapping." % (self.get_name(), withEllipse.get_name()))
            return cpt
            # """

        #         epsilon = 1e-4 # 1e-3 - 1e-4 according to Torquato
        #         # step (d)
        #         # FIXME: secant method is used atm, because I suspect that the Torquato time derivative is not correct... at least my calculation of the torquato derivative
        #         F0dot = (self.overlap_potential(withEllipse, lambda0, tstart+1e-4) - F0)/1e-4 # I don't trust the fdotAB implementation
        #         #""" # FIXME: should be activated
        #
        #         if abs(F0) < epsilon: # ellipses are overlapping or nearly touching
        #             #FIXME: Fdot not working properly!!
        #             #print("F0=",F0, "; F0dot=",F0dot)
        #             if F0dot >= 0: # ellipse distance increasing
        #                 F0 = epsilon
        # #                 print("swell calc3",self._id,withEllipse._id,"ellipse distance increasing")
        #             else: # ellipses are approaching
        #                 cp1 = self.overlap_collision_point(withEllipse, lambda0, tstart)
        #                 cp1_unwrap = self.cell.cellspace.system.boundary.unwrap(cp1)
        #                 cp2 = self.overlap_normal_vector(withEllipse,lambda0,tstart)
        #                 cp2_unwrap = self.cell.cellspace.system.boundary.unwrap(cp2)
        #                 cp3 = np.append(cp1_unwrap,cp2_unwrap)
        #                 cp = np.append(cp3,lambda0)
        # #                 rc = self.position + self.velocity*tstart + (1.0-lambda0)*np.matmul(cp1[0],cp2_unwrap)
        # #                 rc_unwrap = self.cell.cellspace.system.boundary.unwrap(rc)
        # #                 cp = np.append(rc_unwrap,cp2_unwrap)
        #                 cpt = np.concatenate(([tstart],cp), axis=0)
        #                 #print("EXIT2: cpoint lies on boundary of ellipses: ", self.E(cpt), "\tand ", withEllipse.E(cpt), ";\toverlap: ",self.overlap(withEllipse, 0))
        #                 withEllipse.position = oldpos2
        # #                 print("swell calc4",self._id,withEllipse._id,"tc",tstart,"F0",F0)
        #
        #                 return cpt
        # #         print(self._id,withEllipse._id,F0,F0dot,F0dot2,"going inside tc_finder3")
        #         try:
        # #             if (self._id,withEllipse._id) == (23,26):
        # #                 print("#1 going to tc_finder3",self._id,withEllipse._id,"F0",F0,"tstart",tstart,"tend",tend)
        #
        #             dtswell = self.cell.cellspace.system.system_properties.swelling_rate
        #             aspectratio = self.cell.cellspace.system.system_properties.aspect_ratio
        #
        #             if dtswell is None:
        #                 dtswell = 0.
        #
        #             tcf = tc_finder3(ellipse1=self, ellipse2=withEllipse, tstart=tstart, tend=tend, dtswell = dtswell, aspectratio = aspectratio)
        #
        #         except RuntimeError as runerr:
        #             val = runerr.args
        #             if val is not None:
        #                 withEllipse.position = oldpos2
        #                 raise ConvergenceError(val)
        #             else:
        #                 withEllipse.position = oldpos2
        #                 return None
        #         else:
        #             if tcf is not None:
        #                 (tc, lambda0, series, rte) = tcf
        #
        #                 if isfinite(tc) and tc >= 0 and self.overlap_potential(withEllipse, lambda0, tc) < epsilon:
        # #                     cp = self.overlap_collision_point(withEllipse, lambda0, tc)
        #                     cp1 = self.overlap_collision_point(withEllipse, lambda0, tc)
        #                     cp1_unwrap = self.cell.cellspace.system.boundary.unwrap(cp1)
        #                     cp2 = self.overlap_normal_vector(withEllipse,lambda0,tc)
        #                     cp2_unwrap = self.cell.cellspace.system.boundary.unwrap(cp2)
        # #                     rc = self.position + self.velocity*tc + (1.0-lambda0)*np.matmul(cp1[0],cp2_unwrap)
        # #                     rc_unwrap = self.cell.cellspace.system.boundary.unwrap(rc)
        #                     cp3 = np.append(cp1_unwrap,cp2_unwrap)
        #                     cp = np.append(cp3,lambda0)
        # #                     cp = np.append(cp1_unwrap,cp2_unwrap)
        #                     cpt = np.concatenate(([tc],cp), axis=0)
        # #                     print("swell calc5",self.major,self.minor,self._id,withEllipse._id,"tc",tc,"F0dot",F0dot)#"overlap potential",self.overlap_potential(withEllipse,lambda0,tc))
        #
        #                     withEllipse.position = oldpos2
        #                     return cpt
        #
        #         withEllipse.position = oldpos2
        return None  # FIXME: return $t_c^tilde$ instead

    def circle_approximation(self, withEllipse):
        """
        checked -Aiyin
        Approximate maximum and minimum collision time of two ellipses by calculating incircle and circumcircle overlaps.
        
        This method does not contain any program flow logic, i.e. does not raise any errors if times are negative, etc.
        """
        asum = self.major + withEllipse.major
        bsum = self.minor + withEllipse.minor
        csum = self.minor2 + withEllipse.minor2

        if bsum < csum:
            bsum = csum

        r = -self.cell.cellspace.system.boundary.delta_dir(self.position, withEllipse.position)
        raiyin = self.position - withEllipse.position
        r2shifted = self.position - r
        v = self.velocity - withEllipse.velocity
        b = np.dot(r, v)
        rr = np.dot(r, r)
        rr2 = np.dot(raiyin, raiyin)
        vv = np.dot(v, v)
        argouter = b ** 2 - vv * (rr - asum ** 2)
        arginner = b ** 2 - vv * (rr - bsum ** 2)

        # rel_angvel = np.linalg.norm(np.cross(self.angvel,self.position)-np.cross(withEllipse.angvel,withEllipse.position))

        # Finding t_min
        if sqrt(rr) <= asum:  # circumference overlapping
            if vv == 0:  # center of mas not moving
                tmin = 0
                tmax = 1000  # just a big number greater than the brownian time
            else:  # center of mass moving
                if b >= 0:  # further apart
                    tmin = 0
                    tmax = asum / vv  # time it takes to move outside each others major axes diameter
                else:
                    tmin = 0
                    tmax = (-b + sqrt(argouter)) / vv  # corresponds to the time it takes to move further apart again
        else:  # circumference not overlapping
            if vv == 0:  # not moving
                return None
            else:
                if b >= 0:  # moving further apart
                    return None
                else:  # moving closer together
                    if argouter >= 0:
                        tmin = (-b - sqrt(argouter)) / vv  # smaller root corresponds to the time of overlap
                        tmax = (-b + sqrt(argouter)) / vv
                    else:
                        return None

        return (tmin, tmax)

    """
    Serialization
    """

    def set_summary(self, bmf_data):
        self.major = bmf_data.pop('major')
        self.minor = bmf_data.pop('minor')
        self.minor2 = bmf_data.pop('minor2')
        Orientable.set_summary(self, bmf_data)

    def get_summary(self):
        a = Orientable.get_summary(self)
        a['major'] = self.major
        a['minor'] = self.minor
        a['minor2'] = self.minor2
        return a

    def get_summary_dynamics(self):
        a = Orientable.get_summary_dynamics(self)
        return a

    def get_summary_statics(self):
        a = Orientable.get_summary_statics(self)
        a['major'] = self.major
        a['minor'] = self.minor
        a['minor2'] = self.minor2
        return a

    def __init__(self, **kwargs):
        """
        @param kwargs can be either keywords params:
            @param position,velocity 2D vector
            @param angle,angvel,major,minor scalar value
        or
            @param handle A handle to a dataset to load corresponding data from.
        """
        # check if loading of cython modules worked and set functions accordingly

        if USE_CYTHON:
            #             self.E = self.__E_cython
            #             self.gradE = self.__gradE_cython
            #             self._project_on = self._project_on_cython
            #             self._get_normal_projection_interval = self._get_normal_projection_interval_cython
            self.overlap_potential = self.__f_cython
            self.overlap_potential_dt = self.__fdot_cython
            #             self.overlap_potential_dlambda0 = self.__h_cython
            self.overlap_potential_dlambda0 = self.__fABdlambda_cython
            self.overlap_collision_point = self.__rC_cython
            self.overlap_normal_vector = self.__nC_cython
            self.overlap_velocity_vector = self.__vC_cython
            self.overlap_potential_d2lambda0 = self.__fABd2lambda_cython
        # else:
        #             self.E = self.__E
        #             self.gradE = self.__gradE

        major = kwargs.pop('major')  # pop(..) throws an exception if key does not exist
        minor = kwargs.pop('minor')
        minor2 = kwargs.pop('minor2')
        self.major = major
        self.minor = minor
        self.minor2 = minor2
        #         if self.cell.cellspace.system.system_properties.swelling_rate is not None:
        #             self.swellrate = self.cell.cellspace.system.system_properties.swelling_rate
        #         else:
        #             self.swellrate = 0.
        #         self.axes = major, minor, minor2
        Orientable.__init__(self, **kwargs)

    #     @property
    #     def axes(self):
    #         return (self._major, self._minor, self._minor2)
    #
    #     @axes.setter
    #     def axes(self, axes):
    #         (major, minor, minor2) = axes
    #         self.__axescheck(major, minor)
    #         self._major = major
    #         self._minor = minor
    #         self._minor2 = minor2
    #
    #     @property
    #     def major(self):
    #         return self._major
    #
    #     """
    #     @major.setter
    #     def major(self, value):
    #         self.__axescheck(value, self.minor)
    #         self._major = value
    #     """
    #
    #     @property
    #     def minor(self):
    #         return self._minor
    #
    #     """
    #     @minor.setter
    #     def minor(self, value):
    #         self.__axescheck(self.minor, value)
    #         self._minor = value
    #     """
    #
    #     @property
    #     def minor2(self):
    #         return self._minor2
    #
    #     """
    #     @minor.setter
    #     def minor2(self, value):
    #         self.__axescheck(self.minor2, value)
    #         self._minor2 = value
    #     """
    #
    #     @property
    def volume(self):
        # NOTE: major and minor are here the full axes, not half
        return 4. * np.pi * self.major ** 3 / 3.


#
#     @area.setter
#     def area(self, area):
#         """
#         This method should not be used upon construction and is just implemented
#         for completeness. It scales the axes linearly to match the desired 'area'.
#         @param area Desired area to which the ellipse should be linearly shrinked or expanded to.
#         """
#         gamma = area / self.area
#         self.major = self.major * sqrt(gamma)
#         self.minor = self.minor * sqrt(gamma)
#  
#  
#     def __axescheck(self, major, minor):
#         #if major is not None and minor is not None:
#         nmajor = np.linalg.norm(np.array(major))
#         nminor = np.linalg.norm(np.array(minor))
#         if nmajor > 0 and nminor > 0:
#             if nmajor < nminor:
#                 raise ValueError("Ellipse major < minor! (major semi-minor axis is smaller than semi-minor axis)")
#         else:
#             raise ValueError("Ellipse axes must be > 0!")
# 
#

class SphereOrientable(Orientable):
    def mass(self):
        return np.float64(
            1)  # np.pi * self.major * self.minor # FIXME: only valid for constant density # FIXME: proportionality constant...

    def update(self):

        deltaT = Orientable.update(self)
        return deltaT

    def collision(self, with_circle, collision_point):

        self.lastCollision[with_circle] = self.cell.cellspace.system.system_properties.localtime()
        with_circle.lastCollision[self] = self.lastCollision[with_circle]
        sigma = self.radius + with_circle.radius  # should be the same as norm(delta_r)
        delta_r = self.cell.cellspace.system.boundary.delta_dir(with_circle.position, self.position)  # r1-r2
        delta_v = self.velocity - with_circle.velocity  # v1-v2
        delta_rv = np.float64(np.dot(delta_r, delta_v))

        # self.velocity = self.velocity - (delta_r * delta_rv/sigma**2) * 2/(1/self.mass() + 1/with_circle.mass())
        # with_circle.velocity = with_circle.velocity + (delta_r * delta_rv/sigma**2) * 2/(1/self.mass() + 1/with_circle.mass())
        self.velocity = np.float64(self.velocity - 2 * delta_rv / (
            (1 + self.mass() / with_circle.mass()) * sigma) * delta_r / sigma )
        with_circle.velocity = np.float64( with_circle.velocity + 2 * delta_rv / (
            (1 + with_circle.mass() / self.mass()) * sigma) * delta_r / sigma )
        # TODO put angular velocity updates here -sophia SOPHIA

    def collision_for_swelling(self, withEllipse,
                               collisionPoint):  # FIXME: NOTE withEllipse! Add support for other geometries!
        """
        Calculate and apply the results of a collision of the current 'Ellipse' with the 'withEllipse' at a certain 'collisionPoint'.
        @param withEllipse Reference to another Ellipse as collision target
        @param collisionPoint Coordinates of the collision point
        Now collisionPoint also includes normal velctor at point of collision -A
        Note omegas are in body coordinates
        """
        self.lastCollision[withEllipse] = self.cell.cellspace.system.system_properties.localtime()
        withEllipse.lastCollision[self] = self.lastCollision[withEllipse]

        rc = collisionPoint[0:3]  # point of collision, already unwrapped
        normal = collisionPoint[
                 3:6]  # gives the end point of the normal vector wrt to with the origin at rc, this has to be translated to the origin of the simbox
        lambda0 = collisionPoint[-1]

        normal_min = self.cell.cellspace.system.boundary.delta_dir(np.array([0., 0., 0.]),
                                                                   normal)  # translated vector from rc to the simbox
        e_perpendicular = normal_min / np.linalg.norm(normal_min)

        #         e_perpendicular_body = lambda obj: pytobody(obj.angle,e_perpendicular)

        #         print("e_perpendicular",e_perpendicular,"collision point",rc,"va",np.linalg.norm(self.velocity),np.linalg.norm(withEllipse.velocity))
        rcx = lambda obj: self.cell.cellspace.system.boundary.delta_dir(rc, obj.position)
        #         rcx_body = lambda obj: pytobody(obj.angle,rcx(obj))

        e_parallel = lambda obj: self.cell.cellspace.system.boundary.delta_dir(np.array([0., 0., 0.]),
                                                                               np.cross(rcx(obj), e_perpendicular))
        #         e_parallel_body = lambda obj: self.cell.cellspace.system.boundary.delta_dir(np.array([0.,0.,0.]),np.cross(rcx_body(obj),e_perpendicular_body(obj)))

        #         s1 = self.angle[3] # get last element of the array
        #         s2 = withEllipse.angle[3]
        #         p1 = self.angle[:3] # get first three elements of the array
        #         p2 = withEllipse.angle[:3]
        #
        #         P1 = np.array([[0,-p1[2],p1[1]],[p1[2],0,-p1[0]],[-p1[1],p1[0],0]])
        #         P2 = np.array([[0,-p2[2],p2[1]],[p2[2],0,-p2[0]],[-p2[1],p2[0],0]])
        #
        #         Q1 = 2.* (p1*p1.transpose() - s1*P1 + (s1*s1 - 0.5)*np.identity(3) )
        #         Q2 = 2.* (p2*p2.transpose() - s2*P2 + (s2*s2 - 0.5)*np.identity(3) )

        swellrate = self.cell.cellspace.system.system_properties.swelling_rate

        if swellrate is not None:

            swell_matrix = np.identity(3) * swellrate
            s = lambda obj: obj.angle[3]
            p = lambda obj: obj.angle[:3]

            P = lambda obj: np.array(
                [[0, -p(obj)[2], p(obj)[1]], [p(obj)[2], 0, -p(obj)[0]], [-p(obj)[1], p(obj)[0], 0]])
            Q = lambda obj: 2. * (
                p(obj) * p(obj).transpose() - s(obj) * P(obj) + (s(obj) * s(obj) - 0.5) * np.identity(3))
            Oinv = lambda obj: np.array([[1. / obj.major, 0, 0], [0, 1. / obj.minor, 0], [0, 0, 1. / obj.minor2]])
            Gamma = lambda obj: Q(obj).transpose() * (Oinv(obj) * swell_matrix) * Q(obj)

            v_term = lambda obj: obj.velocity + np.cross(rcx(obj), obj.angvel) + np.matmul(Gamma(obj), -rcx(obj))

        else:

            v_term = lambda obj: obj.velocity + np.cross(rcx(obj), obj.angvel)

            #         Oinv1 = np.array([[1./self.major,0,0],[0,1./self.minor,0],[0,0,1./self.minor2]])
            #         Oinv2 = np.array([[1./withEllipse.major,0,0],[0,1./withEllipse.minor,0],[0,0,1./withEllipse.minor2]])

            #         Gamma1 = Q1.transpose()*(Oinv1*swell_matrix)*Q1
            #         Gamma2 = Q2.transpose()*(Oinv2*swell_matrix)*Q2


            #         v_term = lambda obj: obj.velocity + np.cross(rcx(obj),obj.angvel) + np.matmul(Gamma(obj),-rcx(obj))

        Theta_term = lambda obj: obj.mass() * (
            obj.minor ** 2 + obj.minor ** 2) / 5.  # lambda obj: 2.0 * obj.mass() * obj.major**2 / 5.0 #Ia,Ib

        vn = np.matmul(e_perpendicular.transpose(),
                       v_term(withEllipse) - v_term(self))  # np.dot(e_perpendicular,v_term(withEllipse) - v_term(self))

        denom = (
            1. / self.mass() + 1. / withEllipse.mass() + (np.linalg.norm(e_parallel(self)) ** 2) / Theta_term(self) + (
                np.linalg.norm(e_parallel(withEllipse)) ** 2) / Theta_term(withEllipse))

        del_pab = 2.0 * vn / denom

        # just to apply shortest image distance befor F0 computation
        r = -self.cell.cellspace.system.boundary.delta_dir(self.position, withEllipse.position)
        r2shifted = self.position - r
        oldpos2 = np.copy(withEllipse.position)
        withEllipse.position = r2shifted

        print("COLLISION between", self._id, withEllipse._id, "del_pab", del_pab, "F0",
              self.overlap_potential(withEllipse, lambda0, 0.), "lambda", lambda0)

        #         withEllipse.position = oldpos2

        old_veloc_a = self.velocity
        old_veloc_b = withEllipse.velocity
        old_angvel_a = self.angvel
        old_angvel_b = withEllipse.angvel

        del_pab_vector = del_pab * e_perpendicular

        if del_pab < 0.:

            #             print("original",self.angvel,withEllipse.angvel)
            self.velocity = old_veloc_a + del_pab * e_perpendicular / self.mass()
            withEllipse.velocity = old_veloc_b - del_pab * e_perpendicular / withEllipse.mass()
            self.angvel = old_angvel_a - del_pab * e_parallel(self) / Theta_term(self)
            withEllipse.angvel = old_angvel_b + del_pab * e_parallel(withEllipse) / Theta_term(withEllipse)

            lambda1 = self.major / (withEllipse.major + self.major)
            fz = lambda l: self.overlap_potential_dlambda0(withEllipse, l, 1e-8)
            lambda1 = newton(func=fz, x0=lambda1)

            F0dot = self.overlap_potential(withEllipse, lambda1, 1e-8) - self.overlap_potential(withEllipse, lambda0,
                                                                                                0.)

            print("F0", self.overlap_potential(withEllipse, lambda0, 0.), "F0dot", F0dot)

        else:
            print("delpab > 0")

        withEllipse.position = oldpos2

    def correction(self, withEllipse):
        print("correction move is being applied", self._id, withEllipse._id)
        """
        Calculate and apply a correction move between the current and 'withEllipse' particles.
        @param withEllipse Correction with this ellipse is desired
        """
        # Set currect time step as 'last collision time'
        self.lastCollision[withEllipse] = self.cell.cellspace.system.system_properties.localtime()
        withEllipse.lastCollision[self] = self.lastCollision[withEllipse]

        r = -self.cell.cellspace.system.boundary.delta_dir(self.position, withEllipse.position)

        v1 = self.velocity
        v2 = withEllipse.velocity
        vsum = self.cell.cellspace.system.boundary.unwrap(v1 + v2)
        vdiff = self.cell.cellspace.system.boundary.delta_dir(v1, v2)
        #         unit_vsum = vsum/np.linalg.norm(vsum)

        # let e_parallel = e1, here you just want to construct some 3 random new euclidean coordinates -Aiyin
        e1 = r / np.linalg.norm(r)  # first unit vector parallel to the center of mass separation

        e2 = np.cross(e1, vsum)

        e2 = e2 / np.linalg.norm(e2)  # unit vector 2
        e3 = np.cross(e1, e2)  # unit vector 3

        if not self.pinned and not withEllipse.pinned:  # both not pinned
            m1 = self.mass()
            m2 = withEllipse.mass()
            self.velocity = 1 / (2 * m1) * (
                np.dot(e1, (np.dot(m1 * v1 + m2 * v2, e1) + np.linalg.norm(m1 * v1 - m2 * v2))) + np.dot(e2, np.dot(
                    m1 * v1 + m2 * v2, e2)) + np.dot(e3, np.dot(m1 * v1 + m2 * v2, e3)))
            withEllipse.velocity = 1 / (2 * m2) * (
                np.dot(e1, (np.dot(m1 * v1 + m2 * v2, e1) - np.linalg.norm(m2 * v1 - m2 * v2))) + np.dot(e2, np.dot(
                    m1 * v1 + m2 * v2, e2)) + np.dot(e3, np.dot(m1 * v1 + m2 * v2, e3)))
        elif self.pinned and not withEllipse.pinned:  # only self pinned
            withEllipse.velocity = 0.5 * (
                np.dot(e_parallel, (np.dot(v2, e_parallel) - np.linalg.norm(v2))) + np.dot(e_perpendicular,
                                                                                           np.dot(v2, e_perpendicular)))
        elif not self.pinned and withEllipse.pinned:  # only withEllipse pinned
            self.velocity = 0.5 * (
                np.dot(e_parallel, (np.dot(v1, e_parallel) + np.linalg.norm(v1))) + np.dot(e_perpendicular,
                                                                                           np.dot(v1, e_perpendicular)))
        else:  # both pinned
            pass

    def overlap(self, withCircle, dt):
        lambda0 = self.radius / (self.radius + withCircle.radius)
        fdl0 = self.overlap_potential_dlambda0(withCircle, 0, dt)
        fdl1 = self.overlap_potential_dlambda0(withCircle, 1, dt)
        if fdl0 == 0 or fdl1 == 0:
            psi_term = -1  # overlap #ugly
        else:
            try:
                #                 print("position",self.position,withEllipse.position)
                #                 print("velocity",self.velocity,withEllipse.velocity)
                #                 print("angle",self.angle,withEllipse.angle)
                #                 print("angvel",self.angvel,withEllipse.angvel)
                fdlz = lambda l: self.overlap_potential_dlambda0(withCircle, l, dt)
                b = newton(func=fdlz, x0=lambda0)  # secant method

                if 0 < b < 1:
                    lambda0 = b
                    psi_term = self.overlap_potential(withCircle, lambda0, dt)

            except RuntimeError:  # should not happen as we have a concave function and checked signs above, but just leave it in case of numerical failure
                print("Error")

        return (psi_term < 0.)

    """
    Ellipse overlap functions according to Torquato, Perram & Wertheim
    """

    def overlap_potential(self, withEllipse, lambda0,
                          time):  # set f in constructor to cython version of python version if available
        raise NotImplementedError()

    def __f_cython(self, withCircle, lambda0, time):
        if isfinite(time):
            qnewA = pyquat(self.angle, self.angvel, time)
            qnewB = pyquat(withCircle.angle, withCircle.angvel, time)
            #             if (self._id,withEllipse._id) == (62,126):
            #                 print("time",time,np.linalg.norm(self.position - withEllipse.position))
            #                 print("position",self.position + self.velocity * time,withEllipse.position + withEllipse.velocity*time)
            #                 print("angle new",qnewA,qnewB)
            #                 print("angvel",self.angvel,withEllipse.angvel)
            #                 print("F0",pyfAB(lambda0, self.major, withEllipse.major, self.minor, withEllipse.minor, self.minor2, withEllipse.minor2, qnewA, qnewB, self.position + self.velocity * time, withEllipse.position + withEllipse.velocity * time))
            swellrate = self.cell.cellspace.system.system_properties.swelling_rate
            aspectratio = self.cell.cellspace.system.system_properties.aspect_ratio
            if swellrate is not None:
                majorA = self.major + swellrate * time
                majorB = withCircle.major + swellrate * time
                minorA = majorA / aspectratio
                minorB = majorB / aspectratio
            else:
                majorA = self.major
                majorB = withCircle.major
                minorA = self.minor
                minorB = withCircle.minor

            return pyfAB(lambda0, majorA, majorB, minorA, minorB, minorA, minorB, qnewA, qnewB,
                         self.position + self.velocity * time, withCircle.position + withCircle.velocity * time)
        else:
            return float('inf')  # FIXME

    # derivative of overlap_potential by time
    def overlap_potential_dt(self, withEllipse, lambda0,
                             time):  # set h in constructor to cython version of python version if available
        raise NotImplementedError()

    def __fdot_cython(self, withEllipse, lambda0, time):
        qnewA = pyquat(self.angle, self.angvel, time)
        qnewB = pyquat(withEllipse.angle, withEllipse.angvel, time)
        swellrate = self.cell.cellspace.system.system_properties.swelling_rate
        aspectratio = self.cell.cellspace.system.system_properties.aspect_ratio
        if swellrate is not None:
            majorA = self.major + swellrate * time
            majorB = withEllipse.major + swellrate * time
            minorA = majorA / aspectratio
            minorB = majorB / aspectratio
        else:
            majorA = self.major
            majorB = withEllipse.major
            minorA = self.minor
            minorB = withEllipse.minor

        return pyfdotAB(lambda0, majorA, majorB, minorA, minorB, minorA, minorB, qnewA, qnewB,
                        self.position + self.velocity * time, withEllipse.position + withEllipse.velocity * time,
                        self.velocity, withEllipse.velocity, self.angvel, withEllipse.angvel)

    # derivative of overlap_potential by lambda0
    def overlap_potential_dlambda0(self, withEllipse, lambda0,
                                   time):  # set h in constructor to cython version of python version if available
        raise NotImplementedError()

    #     def __h_cython(self, withEllipse, lambda0, time):
    #         qnewA = pyquat(self.angle,self.angvel,time)
    #         qnewB = pyquat(withEllipse.angle, withEllipse.angvel, time)
    #         return pyhAB(lambda0, self.major, withEllipse.major, self.minor, withEllipse.minor, self.minor2, withEllipse.minor2, qnewA, qnewB, self.position + self.velocity * time, withEllipse.position + withEllipse.velocity * time)

    def __fABdlambda_cython(self, withEllipse, lambda0, time):
        qnewA = pyquat(self.angle, self.angvel, time)
        qnewB = pyquat(withEllipse.angle, withEllipse.angvel, time)
        swellrate = self.cell.cellspace.system.system_properties.swelling_rate
        aspectratio = self.cell.cellspace.system.system_properties.aspect_ratio
        if swellrate is not None:
            majorA = self.major + swellrate * time
            majorB = withEllipse.major + swellrate * time
            minorA = majorA / aspectratio
            minorB = majorB / aspectratio
        else:
            majorA = self.major
            majorB = withEllipse.major
            minorA = self.minor
            minorB = withEllipse.minor

        return pyhAB(lambda0, majorA, majorB, minorA, minorB, minorA, minorB, qnewA, qnewB,
                     self.position + self.velocity * time, withEllipse.position + withEllipse.velocity * time)

    def overlap_potential_d2lambda0(self, withEllipse, lambda0,
                                    time):  # set h in constructor to cython version of python version if available
        raise NotImplementedError()

    def __fABd2lambda_cython(self, withEllipse, lambda0, time):
        qnewA = pyquat(self.angle, self.angvel, time)
        qnewB = pyquat(withEllipse.angle, withEllipse.angvel, time)
        swellrate = self.cell.cellspace.system.system_properties.swelling_rate
        aspectratio = self.cell.cellspace.system.system_properties.aspect_ratio
        if swellrate is not None:
            majorA = self.major + swellrate * time
            majorB = withEllipse.major + swellrate * time
            minorA = majorA / aspectratio
            minorB = majorB / aspectratio
        else:
            majorA = self.major
            majorB = withEllipse.major
            minorA = self.minor
            minorB = withEllipse.minor

        return pydhAB(lambda0, majorA, majorB, minorA, minorB, minorA, minorB, qnewA, qnewB,
                      self.position + self.velocity * time, withEllipse.position + withEllipse.velocity * time)

    def overlap_collision_point(self, withEllipse, lambda0, time):
        raise NotImplementedError()

    def __rC_cython(self, withEllipse, lambda0, time):
        qnewA = pyquat(self.angle, self.angvel, time)
        qnewB = pyquat(withEllipse.angle, withEllipse.angvel, time)
        swellrate = self.cell.cellspace.system.system_properties.swelling_rate
        aspectratio = self.cell.cellspace.system.system_properties.aspect_ratio
        if swellrate is not None:
            majorA = self.major + swellrate * time
            majorB = withEllipse.major + swellrate * time
            minorA = majorA / aspectratio
            minorB = majorB / aspectratio
        else:
            majorA = self.major
            majorB = withEllipse.major
            minorA = self.minor
            minorB = withEllipse.minor
        return pyrC(lambda0, majorA, majorB, minorA, minorB, minorA, minorB, qnewA, qnewB,
                    self.position + self.velocity * time, withEllipse.position + withEllipse.velocity * time)

    def __nC_cython(self, withEllipse, lambda0, time):
        qnewA = pyquat(self.angle, self.angvel, time)
        qnewB = pyquat(withEllipse.angle, withEllipse.angvel, time)
        swellrate = self.cell.cellspace.system.system_properties.swelling_rate
        aspectratio = self.cell.cellspace.system.system_properties.aspect_ratio
        if swellrate is not None:
            majorA = self.major + swellrate * time
            majorB = withEllipse.major + swellrate * time
            minorA = majorA / aspectratio
            minorB = majorB / aspectratio
        else:
            majorA = self.major
            majorB = withEllipse.major
            minorA = self.minor
            minorB = withEllipse.minor

        return pynC(lambda0, majorA, majorB, minorA, minorB, minorA, minorB, qnewA, qnewB,
                    self.position + self.velocity * time, withEllipse.position + withEllipse.velocity * time)

    def __vC_cython(self, withEllipse, lambda0, time):
        qnewA = pyquat(self.angle, self.angvel, time)
        qnewB = pyquat(withEllipse.angle, withEllipse.angvel, time)
        swellrate = self.cell.cellspace.system.system_properties.swelling_rate
        aspectratio = self.cell.cellspace.system.system_properties.aspect_ratio
        if swellrate is not None:
            majorA = self.major + swellrate * time
            majorB = withEllipse.major + swellrate * time
            minorA = majorA / aspectratio
            minorB = majorB / aspectratio
        else:
            majorA = self.major
            majorB = withEllipse.major
            minorA = self.minor
            minorB = withEllipse.minor

        return pyvC(lambda0, majorA, majorB, minorA, minorB, minorA, minorB, qnewA, qnewB,
                    self.position + self.velocity * time, withEllipse.position + withEllipse.velocity * time)

    def tpw_collision_time(self, withCircle, T=None):
        # FIXME sophia: normalization, boundary conditions and   / should be 1 but its 1.5...
        """
        self.E not included
        Collision time algorithm by Torquato, Perram and Wertheim
        @param withCircle check collision with 'withEllipse'
        @param T check for collisions in the time interval [0,T]
        """
        r = -self.cell.cellspace.system.boundary.delta_dir(self.position,
                                                           withCircle.position) # SOPHIA gives minimum distance concidering boundary conditions
        rnorm = self.cell.cellspace.system.boundary.delta_norm(self.position,
                                                               withCircle.position)  # SOPHIA wählt kleinsten Abstand
        r2shifted = self.position - r
        oldpos2 = np.copy(withCircle.position)
        withCircle.position = r2shifted

        if rnorm < self.radius + withCircle.radius:
            raise ParticleOverlap("%s and %s are overlapping." % (self.get_name(), withCircle.get_name()))

        coll_time = self.circle_approximation(withCircle)  # gibt collision time

        if coll_time is None:
            withCircle.position = oldpos2
            return None

        # SOPHIA
        # 1.) update
        x = np.float64(self.position + coll_time * self.velocity)
        xneighbour = np.float64(withCircle.position + coll_time * withCircle.velocity)
        # 2.) normal vector
        normal_vector = xneighbour - x  # Sophia, check if this is minimum image FIXME
        # 3.) collision point
        collision_point = np.float64(x + normal_vector / 2.)
        if np.linalg.norm(normal_vector) != 0:
            normal_vector = np.float64(normal_vector / np.linalg.norm(normal_vector))

        cp1 = np.append(collision_point, normal_vector)
        cp = np.concatenate(([coll_time], cp1), axis=0)

        withCircle.position = oldpos2

        return cp

    def find_swelling_time(self, withEllipse, T=None):
        """
        checked -Aiyin but self.E not included

        @param withEllipse check collision with 'withEllipse'
        @param T check for collisions in the time interval [0,T]
        Returns collision time and collision point (and maybe normal vector is useful too) -A
        """
        # if self.overlap(withEllipse, 0) and self.overlap(withEllipse, 1e-10):
        #    raise ParticleOverlap("%s and %s are overlapping." % (self.get_name(), withEllipse.get_name()))
        ct_approx = self.circle_approximation(withEllipse)

        if ct_approx is not None:
            (tmin, tmax) = ct_approx
        else:
            return None

        # steps (a) and (b)
        if T is not None:
            if tmin > T:
                return None
            else:
                tend = min(T, tmax)
        else:
            tend = tmax
        tstart = tmin
        # Move ellipse to avoid discontinuities,
        # NO NEED FOR THIS
        #         print("before",self._id,withEllipse._id,self.position,withEllipse.position)

        r = -self.cell.cellspace.system.boundary.delta_dir(self.position, withEllipse.position)
        r2shifted = self.position - r
        oldpos2 = np.copy(withEllipse.position)
        #         if not np.array_equal(withEllipse.position, r2shifted):
        #             print("check position shifting",withEllipse.position,r2shifted)
        withEllipse.position = r2shifted
        #         print("after",self._id,withEllipse._id,self.position,withEllipse.position)

        # step (c)
        lambda0 = self.major / (self.major + withEllipse.major)  # guess of the initial lambda

        # """ # FIXME: should be activated...
        # initial lambda0 may be the max. of F at time tstart, not the lambda0 initial guess above
        fdl0 = self.overlap_potential_dlambda0(withEllipse, 0, tstart)
        fdl1 = self.overlap_potential_dlambda0(withEllipse, 1, tstart)

        if np.sign(fdl0) == np.sign(fdl1):
            withEllipse.position = oldpos2
            return None
        try:
            fz = lambda l: self.overlap_potential_dlambda0(withEllipse, l, tstart)
            lambda0 = newton(func=fz, x0=lambda0)  # secant method
        # print(self._id,withEllipse._id,lambda0)
        except RuntimeError:
            pass
            # NOTE: if you change lambda0 here, i.e. if you calculate if for another t you'll have to change the overlap checking time-point as well
            # """
            # DEBUG CODE: check if cpt lies on boundary of both ellipses
            # print("collision point lies on boundary of ellipses: ", self.E(cpt), "\tand ", withEllipse.E(cpt), ";\toverlap: ",self.overlap(withEllipse, 0))

        #         #OVERRIDE EVERYTHING
        #         tstart = 0
        #         tend = 0.231481481481
        F0 = self.overlap_potential(withEllipse, lambda0, tstart)
        F0dot2 = self.overlap_potential_dt(withEllipse, lambda0, tstart)

        # catch circles or exact coll.
        if tstart == tend:  # FIXME: maybe move above, maybe simplify, maybe maybe... :D
            cp1 = self.overlap_collision_point(withEllipse, lambda0, tstart)
            cp1_unwrap = self.cell.cellspace.system.boundary.unwrap(cp1)
            cp2 = self.overlap_normal_vector(withEllipse, lambda0, tstart)
            cp2_unwrap = self.cell.cellspace.system.boundary.unwrap(cp2)
            cp3 = np.append(cp1_unwrap, cp2_unwrap)
            cp = np.append(cp3, lambda0)
            cpt = np.concatenate(([tstart], cp), axis=0)
            withEllipse.position = oldpos2
            #             print("swell calc 1,tstart = tend",self._id,withEllipse._id)
            return cpt

        if tstart < 1e-8 and F0 < 0 and self.overlap_potential(withEllipse, lambda0,
                                                               tstart + 1e-8) < 0:  # in case lambda0 is for t=tstart

            withEllipse.position = oldpos2
            # print("OVERLAP, tmin=",tmin," tmax= ",tmax)
            cp1 = self.overlap_collision_point(withEllipse, lambda0, tstart)
            cp1_unwrap = self.cell.cellspace.system.boundary.unwrap(cp1)
            cp2 = self.overlap_normal_vector(withEllipse, lambda0, tstart)
            cp2_unwrap = self.cell.cellspace.system.boundary.unwrap(cp2)
            cp3 = np.append(cp1_unwrap, cp2_unwrap)
            cp = np.append(cp3, lambda0)
            cpt = np.concatenate(([tstart], cp), axis=0)
            print("swell calc2, OVERLAP!", self._id, withEllipse._id, "major", self.major, self.minor, "tc", tstart,
                  "F0", F0)
            #             sys.exit("Overlap detected")
            raise ParticleOverlap("%s and %s are overlapping." % (self.get_name(), withEllipse.get_name()))
            return cpt
            # """

        #         epsilon = 1e-4 # 1e-3 - 1e-4 according to Torquato
        #         # step (d)
        #         # FIXME: secant method is used atm, because I suspect that the Torquato time derivative is not correct... at least my calculation of the torquato derivative
        #         F0dot = (self.overlap_potential(withEllipse, lambda0, tstart+1e-4) - F0)/1e-4 # I don't trust the fdotAB implementation
        #         #""" # FIXME: should be activated
        #
        #         if abs(F0) < epsilon: # ellipses are overlapping or nearly touching
        #             #FIXME: Fdot not working properly!!
        #             #print("F0=",F0, "; F0dot=",F0dot)
        #             if F0dot >= 0: # ellipse distance increasing
        #                 F0 = epsilon
        # #                 print("swell calc3",self._id,withEllipse._id,"ellipse distance increasing")
        #             else: # ellipses are approaching
        #                 cp1 = self.overlap_collision_point(withEllipse, lambda0, tstart)
        #                 cp1_unwrap = self.cell.cellspace.system.boundary.unwrap(cp1)
        #                 cp2 = self.overlap_normal_vector(withEllipse,lambda0,tstart)
        #                 cp2_unwrap = self.cell.cellspace.system.boundary.unwrap(cp2)
        #                 cp3 = np.append(cp1_unwrap,cp2_unwrap)
        #                 cp = np.append(cp3,lambda0)
        # #                 rc = self.position + self.velocity*tstart + (1.0-lambda0)*np.matmul(cp1[0],cp2_unwrap)
        # #                 rc_unwrap = self.cell.cellspace.system.boundary.unwrap(rc)
        # #                 cp = np.append(rc_unwrap,cp2_unwrap)
        #                 cpt = np.concatenate(([tstart],cp), axis=0)
        #                 #print("EXIT2: cpoint lies on boundary of ellipses: ", self.E(cpt), "\tand ", withEllipse.E(cpt), ";\toverlap: ",self.overlap(withEllipse, 0))
        #                 withEllipse.position = oldpos2
        # #                 print("swell calc4",self._id,withEllipse._id,"tc",tstart,"F0",F0)
        #
        #                 return cpt
        # #         print(self._id,withEllipse._id,F0,F0dot,F0dot2,"going inside tc_finder3")
        #         try:
        # #             if (self._id,withEllipse._id) == (23,26):
        # #                 print("#1 going to tc_finder3",self._id,withEllipse._id,"F0",F0,"tstart",tstart,"tend",tend)
        #
        #             dtswell = self.cell.cellspace.system.system_properties.swelling_rate
        #             aspectratio = self.cell.cellspace.system.system_properties.aspect_ratio
        #
        #             if dtswell is None:
        #                 dtswell = 0.
        #
        #             tcf = tc_finder3(ellipse1=self, ellipse2=withEllipse, tstart=tstart, tend=tend, dtswell = dtswell, aspectratio = aspectratio)
        #
        #         except RuntimeError as runerr:
        #             val = runerr.args
        #             if val is not None:
        #                 withEllipse.position = oldpos2
        #                 raise ConvergenceError(val)
        #             else:
        #                 withEllipse.position = oldpos2
        #                 return None
        #         else:
        #             if tcf is not None:
        #                 (tc, lambda0, series, rte) = tcf
        #
        #                 if isfinite(tc) and tc >= 0 and self.overlap_potential(withEllipse, lambda0, tc) < epsilon:
        # #                     cp = self.overlap_collision_point(withEllipse, lambda0, tc)
        #                     cp1 = self.overlap_collision_point(withEllipse, lambda0, tc)
        #                     cp1_unwrap = self.cell.cellspace.system.boundary.unwrap(cp1)
        #                     cp2 = self.overlap_normal_vector(withEllipse,lambda0,tc)
        #                     cp2_unwrap = self.cell.cellspace.system.boundary.unwrap(cp2)
        # #                     rc = self.position + self.velocity*tc + (1.0-lambda0)*np.matmul(cp1[0],cp2_unwrap)
        # #                     rc_unwrap = self.cell.cellspace.system.boundary.unwrap(rc)
        #                     cp3 = np.append(cp1_unwrap,cp2_unwrap)
        #                     cp = np.append(cp3,lambda0)
        # #                     cp = np.append(cp1_unwrap,cp2_unwrap)
        #                     cpt = np.concatenate(([tc],cp), axis=0)
        # #                     print("swell calc5",self.major,self.minor,self._id,withEllipse._id,"tc",tc,"F0dot",F0dot)#"overlap potential",self.overlap_potential(withEllipse,lambda0,tc))
        #
        #                     withEllipse.position = oldpos2
        #                     return cpt
        #
        #         withEllipse.position = oldpos2
        return None  # FIXME: return $t_c^tilde$ instead

    def circle_approximation(self, withCircle):
        """
        SOPHIA

        Approximate collision time of two particles by calculating overlaps.

        This method does not contain any program flow logic, i.e. does not raise any errors if times are negative, etc.
        """
        sigma = self.radius + withCircle.radius

        r = -self.cell.cellspace.system.boundary.delta_dir(self.position, withCircle.position)
        r2shifted = self.position - r
        v = self.velocity - withCircle.velocity
        b = np.float64(np.dot(r, v))
        rr = np.float64(np.dot(r, r))
        vv = np.float64(np.dot(v, v))
        argouter = np.float64(b ** 2 - vv * (rr - sigma ** 2))
        # Finding t_min
        if sqrt(rr) <= sigma:  # circumference overlapping
            t_coll = 0.
        else:  # circumference not overlapping
            if vv == 0:  # not moving
                return None
            else:
                if b >= 0:  # moving further apart
                    return None
                else:  # moving closer together
                    if argouter >= 0:
                        t_coll = np.float64((-b - sqrt(argouter)) / vv )
                    # tmax = (-b + sqrt(argouter)) / vv
                    else:
                        return None
        return t_coll

    """
    Serialization
    """

    def set_summary(self, bmf_data):
        self.major = bmf_data.pop('major')
        self.minor = bmf_data.pop('minor')
        self.minor2 = bmf_data.pop('minor2')
        Orientable.set_summary(self, bmf_data)

    def get_summary(self):
        a = Orientable.get_summary(self)
        a['major'] = self.major
        a['minor'] = self.minor
        a['minor2'] = self.minor2
        return a

    def get_summary_dynamics(self):
        a = Orientable.get_summary_dynamics(self)
        return a

    def get_summary_statics(self):
        a = Orientable.get_summary_statics(self)
        a['major'] = self.major
        a['minor'] = self.minor
        a['minor2'] = self.minor2
        return a

    def __init__(self, **kwargs):

        """
        @param kwargs can be either keywords params:
            @param position,velocity 2D vector
            @param angle,angvel,major,minor scalar value
        or
            @param handle A handle to a dataset to load corresponding data from.
        """
        # check if loading of cython modules worked and set functions accordingly

        if USE_CYTHON:
            #             self.E = self.__E_cython
            #             self.gradE = self.__gradE_cython
            #             self._project_on = self._project_on_cython
            #             self._get_normal_projection_interval = self._get_normal_projection_interval_cython
            self.overlap_potential = self.__f_cython
            self.overlap_potential_dt = self.__fdot_cython
            #             self.overlap_potential_dlambda0 = self.__h_cython
            self.overlap_potential_dlambda0 = self.__fABdlambda_cython
            self.overlap_collision_point = self.__rC_cython
            self.overlap_normal_vector = self.__nC_cython
            self.overlap_velocity_vector = self.__vC_cython
            self.overlap_potential_d2lambda0 = self.__fABd2lambda_cython
        # else:
        #             self.E = self.__E
        #             self.gradE = self.__gradE


        major = kwargs.pop('major')  # pop(..) throws an exception if key does not exist
        minor = kwargs.pop('minor')
        minor2 = kwargs.pop('minor2')
        radius = kwargs.pop('radius')
        self.major = major
        self.minor = minor
        self.minor2 = minor2
        self.radius = radius
        #         if self.cell.cellspace.system.system_properties.swelling_rate is not None:
        #             self.swellrate = self.cell.cellspace.system.system_properties.swelling_rate
        #         else:
        #             self.swellrate = 0.
        #         self.axes = major, minor, minor2
        Orientable.__init__(self, **kwargs)

    #     @property
    #     def axes(self):
    #         return (self._major, self._minor, self._minor2)
    #
    #     @axes.setter
    #     def axes(self, axes):
    #         (major, minor, minor2) = axes
    #         self.__axescheck(major, minor)
    #         self._major = major
    #         self._minor = minor
    #         self._minor2 = minor2
    #
    #     @property
    #     def major(self):
    #         return self._major
    #
    #     """
    #     @major.setter
    #     def major(self, value):
    #         self.__axescheck(value, self.minor)
    #         self._major = value
    #     """
    #
    #     @property
    #     def minor(self):
    #         return self._minor
    #
    #     """
    #     @minor.setter
    #     def minor(self, value):
    #         self.__axescheck(self.minor, value)
    #         self._minor = value
    #     """
    #
    #     @property
    #     def minor2(self):
    #         return self._minor2
    #
    #     """
    #     @minor.setter
    #     def minor2(self, value):
    #         self.__axescheck(self.minor2, value)
    #         self._minor2 = value
    #     """
    #
    #     @property
    def volume(self):
        # NOTE: major and minor are here the full axes, not half
        return 4. * np.pi * self.major ** 3 / 3.

    #
    #     @area.setter
    #     def area(self, area):
    #         """
    #         This method should not be used upon construction and is just implemented
    #         for completeness. It scales the axes linearly to match the desired 'area'.
    #         @param area Desired area to which the ellipse should be linearly shrinked or expanded to.
    #         """
    #         gamma = area / self.area
    #         self.major = self.major * sqrt(gamma)
    #         self.minor = self.minor * sqrt(gamma)
    #
    #
    #     def __axescheck(self, major, minor):
    #         #if major is not None and minor is not None:
    #         nmajor = np.linalg.norm(np.array(major))
    #         nminor = np.linalg.norm(np.array(minor))
    #         if nmajor > 0 and nminor > 0:
    #             if nmajor < nminor:
    #                 raise ValueError("Ellipse major < minor! (major semi-minor axis is smaller than semi-minor axis)")
    #         else:
    #             raise ValueError("Ellipse axes must be > 0!")
    #
    #
    def __str__(self):
        return "SphereOrientable at {} with Orientation {}".format(self.position, self.angle)
