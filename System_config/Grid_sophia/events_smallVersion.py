"""
Created on 04.02.2015

@author: mahe
"""

from collections import namedtuple

import sys #-A
import simpy
import numpy as np

#from bedsim.particle import ParticleOverlap, ConvergenceError,\
 #   ParticleSwellOverlap



"""
# This is a sample "Event-type" simulation where only the particle trajectories are updated every timestep
"""

class ParticleOverlap2(ValueError):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

class Event(object):
    """
    Abstract Event class.
    """
    def process(self, time):
        raise NotImplementedError()
    
    def particle_check(self, particles):
        """
        This will printout all the particles in your system
        """
        done = {}
        for particle in particles:
                print("Particles ", particle.get_name())

    def __init__(self, event_manager, time):
        """ Save EventManager reference and call abstract process method. """
        self.event_manager = event_manager
        self.proc = self.event_manager.env.process(self.process(time))

class SystemEvent(Event):
    """
    Event concerning the whole system, i.e. all particles. E.g. drawing new velocities
    """
    
    def draw_velocities_and_correct_variances(self):
        particles = []
        
        for particle in self.event_manager.system._particles:
            if not particle.pinned:
                particles.append(particle)

        # 2. give 'particles' random speed with mean 0
        pn = len(particles)
        rn = np.random.normal(loc=0, scale=1, size=(pn,6)) # random xvel, yvel, zvel and angvelx, angvely, angvelz for each particle; loc => mu, scale => sigma
        rn -= np.mean(rn, axis=0) # calculate the per particle and per component means in order to correct velocities to sum=0
        for (p, r) in zip(particles, rn): # p is reference to a particle, r is 3d numpy array # FIXME: convert rn tolist()?
            p.velocity = r[:3] # returns a numpy.array(), first three terms of r
            p.angvel   = r[-3:] # returns a number, last three terms of r

        # 3. correct variance of velocities
        velvar = np.array([np.dot(r[:3], r[:3]) for r in rn])
        velvar = np.sum(velvar) / len(particles)
        angvelvar = np.array([np.dot(r[-3:], r[-3:]) for r in rn])
        angvelvar = np.sum(angvelvar) / len(particles)
        for p in particles:
             
            p.velocity *= np.sqrt(1 / velvar) # FIXME: mass
            p.angvel   *= np.sqrt(5 / 2.* p.major / angvelvar) 

    
    def __init__(self, event_manager, time):
        Event.__init__(self, event_manager, time)

"""
# Implementation of the actual events
"""
class UpdateEvent(Event):

    def process(self, time):
        yield self.event_manager.env.timeout(time)

        """ Update particles trajectories """

        [p.update() for p in self.event_manager.system._particles] 
        
        """ Save results to file """
        self.event_manager.system.save_to_file()
        
        self.event_manager.env.process(self.process(time))


    def __init__(self, event_manager, time):

        Event.__init__(self, event_manager, time)



class QuitEvent(SystemEvent):
    def process(self, time):
        """
        1. Update all particles
        2. Write particle properties to file
        """
        yield self.event_manager.env.timeout(time)
        
        print("for quit event saving")
        self.event_manager.system.save_to_file()


class StartEvent(SystemEvent):
    def process(self, time): # time = 0
        """
        The start event. This event predicts first particle events and schedules system events.
        After that StartEvent terminates without rescheduling.
        """
        yield self.event_manager.env.timeout(time)

        self.draw_velocities_and_correct_variances()
        
        self.particle_check(self.event_manager.system._particles)

        for particle in self.event_manager.system._particles:
            self.event_manager.events_per_particle[particle] = [] # initialize per particle event lists with an empty list

        self.event_manager.system.save_to_file()
        
        # time interval where particles are updated
        timestep = 1e-2
        
        UpdateEvent(self.event_manager, timestep)


"""
# Event Manager
Event Manager should be a special implementation of an abstract base class for general
dynamics simulations.
"""
class EventManager(object):
    def start(self):
        StartEvent(self, 0)
        QuitEvent(self, self.system.system_properties.lifetime - self.system.system_properties.summary_timestep/2) # NOTE: simpy simulates from [tstart, tend)
        self.env.run(until=self.system.system_properties.lifetime)

    def __init__(self, system):
        self.system = system # Reference to the system
        self.env = simpy.Environment() # create simpy environment
        self.events_per_particle = {}  # create per particle event lists
