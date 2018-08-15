"""
Created on 04.02.2015

@author: mahe
"""

from collections import namedtuple

import sys #-A
import simpy
import numpy as np

from bedsim.particle import ParticleOverlap, ConvergenceError



"""
# Abstract Event class.
"""
class Event(object):
    """
    Abstract Event class.
    """
    count_collision_predictions = 0 # FIXME: remove or leave for stats

    def process(self, time):
        raise NotImplementedError()
    
    def overlap_check(self, particles):
        """
        Manual overlap check. Used e.g. for checking the initial configuration.
        """
        done = {}
        for particle in particles:
            pid1 = int(particle._id)
            for neigh_cell in [particle.cell] + particle.cell.neighbours: # check for collisions of particle in its own cell + neighbour cells
                for neigh_particle in neigh_cell.particles:
                    pid2 = neigh_particle._id
                    if pid1 != pid2:
                        if pid1 < pid2:
                            key = (pid1, pid2)
                        else:
                            key = (pid2, pid1)
                        if key not in done: #don't check particles twice
                            done[key] = True
                            if particle.overlap(neigh_particle, 0):
                                print("Particles ", particle.get_name(), " and ", neigh_particle.get_name(), " overlap!")
                              #  FIXME: maybe raise an exception
    

    # FIXME: use until instead of brownian step!
    # FIXME: DOES NOT CONSIDER DIFFERENT UPDATE STATES OF THE PARTICLES!!! This is circumvented atm by updating all particles when the state changes
    def predict_future_events(self, particles, ignore_last_collision_time=False):
        """
        Predict and schedule all future events for each particle of list particles.
        @param particles List of particles to predict and schedule particle events for
        """
        tb = self.event_manager.system.system_properties.brownian_timestep
        if tb is None:
            until = 0.1 # FIXME: depends on packing fraction!
        else:
            now = self.event_manager.system.system_properties.localtime() % tb
            until = 0.01 # FIXME: brownian_timestep - (current_time % brownian_timestep)
            #until = tb - np.fmod(now, tb)
        
        done = {}
        epsilon = 1e-8  # ask Aleena
        for particle in particles:
            """ 1. predict cell crossings """
            cc = particle.cell.crossing_time(particle.position, particle.velocity)
            if cc is not None:
                (crossing_time, to_cell , crossing_point) = cc
                if until >= crossing_time > epsilon: # this is to remove negative times in cell crossing predictions check this
                    CellCrossingEvent(self.event_manager, crossing_time, [particle], to_cell, crossing_point)

            """ 2. predict collisions """
            pid1 = int(particle._id)
            
            for neigh_cell in [particle.cell] + particle.cell.neighbours: # check for collisions of particle in its own cell + neighbour cells
                
                for neigh_particle in neigh_cell.particles:
                    #if self.event_manager.system.boundary.delta_norm(particle.position,neigh_particle.position) <= (particle.major +neigh_particle.major + 0.2): ##VERLET list
                        pid2 = int(neigh_particle._id)
                        if pid1 < pid2:
                            key = (pid1, pid2)
                        else:
                            key = (pid2, pid1)

                        if particle is not neigh_particle and key not in done: # don't check particles twice
                            Event.count_collision_predictions += 1
                            done[key] = True
                            try:
                                ct = particle.tpw_collision_time(neigh_particle, T=until)
                                #FIXME SOPHIA vielleicht läfut das Programm ja dann?

                                if ct is not None:

                                    # if collision time algorithm has no upper limit, also compare to Brownian time step if set and discard event if later than Brownian time step
                                    #print(ct[0], particle._id, neigh_particle._id)
                                    #print('time',ct[0], particle._id, neigh_particle._id)
    #                                CollisionEvent(self.event_manager, ct[0], [particle, neigh_particle], ct[1:7])

                                    """ Now check if last collision of these two particles occurred at least 10**-8 timesteps ago. """
                                    """ First check if a collision of the two particles ever occurred in the past """

                                    if neigh_particle in particle.lastCollision:
                                        absolute_collision_time = self.event_manager.system.system_properties.localtime() + ct[0] # absolute collision time # FIXME: change in rescaled system


                                        if absolute_collision_time - particle.lastCollision[neigh_particle] > epsilon or ignore_last_collision_time:
                                            CollisionEvent(self.event_manager, ct[0], [particle, neigh_particle], ct[1:8])
    #                                         print("EVENT: CollisionEvent1 at", self.event_manager.system.system_properties.localtime() + ct[0], "by particles ", [pid1,pid2])
                                    else:
                                        """ no collision in the past => just collide """
                                        CollisionEvent(self.event_manager, ct[0], [particle, neigh_particle], ct[1:8])
    #                                     print("EVENT: CollisionEvent2 at", self.event_manager.system.system_properties.localtime() + ct[0], "by particles ", [pid1,pid2])


                            except ParticleOverlap: # this is actually a ParticleOverlap error inheriting a ValueError
                                """ catch overlap """
    #                             print("inside predict_future_events catch particle overlap")
                                if neigh_particle in particle.lastCollision:
                                    print('overlap: ') ## SOPHIA Juni
                                    print(particle._id,neigh_particle._id, 'distance', np.linalg.norm(particle.position - neigh_particle.position))
                                    now = self.event_manager.system.system_properties.localtime()
                                    if abs(now - particle.lastCollision[neigh_particle]) > epsilon or ignore_last_collision_time:

                                        CorrectionEvent(self.event_manager, 0.0, [particle, neigh_particle])

                                    #print("ignoring last overlap", particle._id,neigh_particle._id) ## SOPHIA Juni
                                else:
                                    CorrectionEvent(self.event_manager, 0.0, [particle, neigh_particle])

    def __init__(self,event_manager, time):
        """ Save EventManager reference and call abstract process method"""
        self.event_manager = event_manager
        self.proc = self.event_manager.env.process(self.process(time))




"""
# General Event classification
"""
class ParticleEvent(Event):
    """
    Event concerning a certain number of particles.
    """
    def deactivate_related(self):
        """
        Deactivates related events, i.e. deactivate all ParticleEvents which concern the particles in self.particles
        """
        for particle in self.particles:
#             print("inside deactivate_related",particle._id)
            for partevent in self.event_manager.events_per_particle[particle]:
                partevent.active = False
            self.event_manager.events_per_particle[particle] = [] # remove obsolete events from events_per_particle lists to simplify list access
            
    def __init__(self, event_manager, time, particles):
        self.active = True
        self.particles = particles
        # add event to eventManager.events_per_particle
        #######[print("FOOOO: ", event_manager.events_per_particle[particle]) for particle in self.particles]
        [event_manager.events_per_particle[particle].append(self) for particle in self.particles]
        Event.__init__(self, event_manager, time)


class SystemEvent(Event):
    """
    Event concerning the whole system, i.e. all particles.
    """
    def draw_velocities(self):
        # consider pinning appropriately
        particles = []
        for particle in self.event_manager.system._particles:
            if not particle.pinned:
                particles.append(particle)

        # draw velocities
        pn = len(particles)
        rn = np.random.normal(loc=0, scale=1, size=(pn,6)) # random xvel, yvel, angvel for each particle; loc => mu, scale => sigma
        rn -= np.mean(rn, axis=0) # calculate the per particle and per component means in order to correct velocities to sum=0
        for (p, r) in zip(particles, rn): # p is reference to a particle, r is 3d numpy array # FIXME: convert rn tolist()?
            p.velocity = r[:3] # returns a numpy.array()
            p.angvel   = r[-3:] # returns a number
            #p.velocity = r[:2] / p.mass() # returns a numpy.array()
            #p.angvel   = r[-1] /(p.mass() * (p.major**2 + p.minor**2)/4) # returns a number # FIXME: schaue ob Trägheitstensor ok ist!!

    
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
            p.angvel   *= np.sqrt(5 / 2.* p.major / angvelvar) # FIXME: mass


    def draw_translational_velocities(self):
        particles = []
        for particle in self.event_manager.system._particles:
            if not particle.pinned:
                particles.append(particle)

        # 2. give 'particles' random speed with mean 0
        pn = len(particles)
        rn = np.random.normal(loc=0, scale=1, size=(pn,3)) # random xvel, yvel, zvel and angvelx, angvely, angvelz for each particle; loc => mu, scale => sigma
        for (p, r) in zip(particles, rn): # p is reference to a particle, r is 3d numpy array # FIXME: convert rn tolist()?
            p.velocity = r[:3] # returns a numpy.array(), first three terms of r


    def draw_rotational_velocities(self):
        particles = []
        for particle in self.event_manager.system._particles:
            if not particle.pinned:
                particles.append(particle)

        # 2. give 'particles' random speed with mean 0
        pn = len(particles)
        rn = np.random.normal(loc=0, scale=1, size=(pn,3)) # random xvel, yvel, zvel and angvelx, angvely, angvelz for each particle; loc => mu, scale => sigma
        for (p, r) in zip(particles, rn): # p is reference to a particle, r is 3d numpy array # FIXME: convert rn tolist()?
            p.angvel   = r[:3] # returns a number, last three terms of r

        for p in particles:
            p.angvel   *= np.sqrt(5 / 2.* p.major**2) # FIXME: mass

    def set_velocities_to_zero_for_swelling(self):
         
        particles = []
        for particle in self.event_manager.system._particles:
            if not particle.pinned:
                particles.append(particle)
  
        # 2. give 'particles' 0 velocity and angular velocity
        #old_traj = [] # Aiyin, do I need to save old trajectories
        for p in particles: # p is reference to a particle
            #old_traj.append([p.velocity,p.angvel])
            p.velocity = np.zeros(3) 
            p.angvel   = np.zeros(3) 

    
    def __init__(self, event_manager, time):
        Event.__init__(self, event_manager, time)


"""
# Implementation of the actual events
"""
class CollisionEvent(ParticleEvent):
    """
    Event for collision of two particles.
    """

    count = 0 # REMOVE

    def process(self, time):

        yield self.event_manager.env.timeout(time)
        """ the following code is executed after timeout """
        if self.active:
            """ perform collision and deactivate related events """
            [p.update() for p in self.event_manager.system._particles] ## FIXME: update neigh or consider different times of last update in predict_future events!!!
            #print("inside collision, position",self.event_manager.system._particles[2].position)
            #self.event_manager.system.save_to_file() ## Save collsions to file
#            if (self.particles[0]._id,self.particles[1]._id) == (20,24):
#                print("before",self.particles[0].velocity,self.particles[1].velocity,self.particles[0]._id,self.particles[1]._id)

            self.particles[0].collision(self.particles[1], self.collision_point) # perform collision calculations
#            if (self.particles[0]._id, self.particles[1]._id) == (20, 24):
#                print("after", self.particles[0].velocity,self.particles[1].velocity,self.particles[0]._id,self.particles[1]._id)
#                sys.exit()
            self.deactivate_related()
            self.predict_future_events(self.particles, ignore_last_collision_time=False) # is prediction of new events needed? => yes
            CollisionEvent.count += 1 # REMOVE
        else:
            """ deactivate process """
            yield self.event_manager.env.event()

    def __init__(self, event_manager, time, particles, collision_point):
        """
        @param event_manager Reference to the EventManager
        @param time Relative time (timeout) when collision happens
        @param particles Involved particles
        @param collision_point Collision coordinate
        """
        self.collision_point = collision_point
        ParticleEvent.__init__(self, event_manager, time, particles)

# class SwellEvent(ParticleEvent):
#     """
#     Event for collision of two particles.
#     """
#     count = 0 # REMOVE
#     def process(self, time):
#         yield self.event_manager.env.timeout(time)
#         """ the following code is executed after timeout """
#         if self.active:
#             #print("ctime: ", self.event_manager.system.system_properties.localtime(), " and waited for ", time)
#             """ perform collision and deactivate related events """
# 
#             [p.swell_update() for p in self.event_manager.system._particles] ## FIXME: update neigh or consider different times of last update in predict_future events!!!
# #             if (self.particles[0]._id,self.particles[1]._id) == (62,126):
#             self.event_manager.system.save_to_file() # placed by Aiyin
#             
#             print("axes before",self.particles[0]._id,self.particles[1]._id,"major",self.particles[0].major,self.particles[0].minor,self.particles[0].minor2,self.particles[].major,self.particles[1].minor,self.particles[1].minor2)
#             self.particles[0].collision_for_swelling(self.particles[1], self.collision_point) # perform collision calculations even if actual velocity is zero? or maybe do a separate algorithm for swelling collision? - Aiyin
#             print("axes after",self.particles[0]._id,self.particles[1]._id,"major",self.particles[0].major,self.particles[0].minor,self.particles[0].minor2,self.particles[].major,self.particles[1].minor,self.particles[1].minor2)
#             print("TIME NOW",self.event_manager.system.system_properties.localtime())
#             self.deactivate_related()
# #             print("inside process Collision, predict future events again")
#             self.predict_future_events(self.particles, ignore_last_collision_time=False) # is prediction of new events needed? => yes
# #             CollisionEvent.count += 1 # REMOVE
#         else:
#             """ deactivate process """
#             yield self.event_manager.env.event()
# 
#     def __init__(self, event_manager, time, particles, collision_point):
#         """
#         @param event_manager Reference to the EventManager
#         @param time Relative time (timeout) when collision happens
#         @param particles Involved particles
#         @param collision_point Collision coordinate
#         """
#         self.collision_point = collision_point
#         ParticleEvent.__init__(self, event_manager, time, particles)

class PredictionEvent(ParticleEvent):
    """
    Event for manually scheduling a prediction of future events (useful if collision time algorithm does not converge)
    """
    def process(self, time):
        yield self.event_manager.env.timeout(time)
#         print("Prediction event started at ", self.event_manager.env.now)
        """ the following code is executed after timeout """
        if self.active:
            """ update and check for future events """
            #[p.update() for p in self.particles] # coord update
            [p.update() for p in self.event_manager.system._particles] # FIXME: update neigh
            #self.deactivate_related() # FIXME: is this correct
            self.predict_future_events([self.particles[0]]) # only predict for one particle, because the other one should be in neigh-list # FIXME: only predict for those two particles!
            #self.predict_future_events(self.particles)
        else:
            """ deactivate process """
            yield self.event_manager.env.event()

    def __init__(self, event_manager, time, particles):
        """
        @param event_manager Reference to the EventManager
        @param time Relative time (timeout) when event happens
        @param particles Involved particles
        """
        ParticleEvent.__init__(self, event_manager, time, particles)


class CorrectionEvent(ParticleEvent):
    count = 0 # REMOVE
    def process(self, time):
        yield self.event_manager.env.timeout(time)
        """ the following code is executed after timeout """
        if self.active:
            #print("\n>> Performing correction")
            """ perform collision and deactivate related events """
            #if time > 0: # the standard use case is that time=0. This implementation, however, supports arbitrary times
            #    [p.update() for p in self.particles]
            #[p.update() for p in self.particles] # just make sure that update is performed (this is indeed not assured)
            [p.update() for p in self.event_manager.system._particles] # FIXME: update neigh
            self.particles[0].correction(self.particles[1])
            self.deactivate_related()
            self.predict_future_events(self.particles) # is prediction of new events needed? => yes
            CorrectionEvent.count += 1
        else:
            """ deactivate process """
            yield self.event_manager.env.event()

    def __init__(self, event_manager, time, particles):
        """
        @param event_manager Reference to the EventManager
        @param time Relative time (timeout) when collision happens
        @param particles Involved particles
        """
        ParticleEvent.__init__(self, event_manager, time, particles)


class CellCrossingEvent(ParticleEvent):
    count = 0 # REMOVE
    def process(self, time):
        
        yield self.event_manager.env.timeout(time)
        
        if self.active:
            [p.update() for p in self.event_manager.system._particles] # FIXME: update neigh
#             print("inside CellCrossingEvent ------- predict again")
            for p in self.particles: # self.particles is a 1 element list, but this looks nicer
                    p.update() # update is done in cell_cross method => is needed here for cumulative position ### FIXME: Position sollte NUR in update-Methode verändert werden!
                    #if p._id == 14:
                    #    print("Particle 14 crossed from ", p.cell," to cell ", self.to_cell)
                    p.cell_cross(self.to_cell, self.crossing_point)
#             print("for cell crossing saving")
#             self.event_manager.system.save_to_file() # placed by Aiyin
            
            self.predict_future_events(self.particles) # is prediction of new events needed? => yes, for crossing particle
            CellCrossingEvent.count += 1
        else:
            """ deactivate process """
            yield self.event_manager.env.event()

    def __init__(self, event_manager, time, particles, to_cell, crossing_point):
        self.to_cell = to_cell
        self.crossing_point = crossing_point
        ParticleEvent.__init__(self, event_manager, time, particles)


class QuitEvent(SystemEvent):
    def process(self, time):
        """
        1. Update all particles
        2. Write particle properties to file
        """
        yield self.event_manager.env.timeout(time)
        print("inside Quit Event after yielding", self.event_manager.env.now)
        print("\nExecuted correction events: ", CorrectionEvent.count)
        print("Executed collision events: ", CollisionEvent.count)
        print("Particle collision predictions: ", Event.count_collision_predictions)
        print("Particle cell crossing predictions: ", CellCrossingEvent.count)
        print("Volume fraction reached", self.event_manager.system.compute_volume_fraction())
        [p.update() for p in self.event_manager.system._particles]
        self.event_manager.system.save_to_file()
        # FIXME: ensure to quit


class StartEvent(SystemEvent):
    def process(self, time): # time = 0
        """
        The start event. This event predicts first particle events and schedules system events.
        After that StartEvent terminates without rescheduling.
        """
        yield self.event_manager.env.timeout(time)
        # check initial configuration for overlap
        self.overlap_check(self.event_manager.system._particles)
        # initialize event_per_particle lists
        for particle in self.event_manager.system._particles:
            self.event_manager.events_per_particle[particle] = [] # initialize per particle event lists with an empty list

        ## save initial configuration to dynamics table
        ## FIXME: maybe all of this can be done more elegantly
        self.event_manager.system.save_to_file() # 2.
        #self.event_manager.system.system_properties.summary_timestep_number += 1
#         self.write_trajectories(self.event_manager.system._particles)        
        #MaintenanceEvent(self.event_manager, self.event_manager.system.system_properties.summary_timestep*100) # FIXME: timestep! 100 seems good for brownian motion
             
        if self.event_manager.system.system_properties.brownian_timestep is None:
#             MaintenanceEvent(self.event_manager, self.event_manager.system.system_properties.summary_timestep*10) # FIXME: timestep!
            MaintenanceEvent(self.event_manager,self.event_manager.system.system_properties.summary_timestep) # FIXME: timestep! #self.event_manager.system.system_properties.summary_timestep
            NewtonianStartEvent(self.event_manager, 0)
            
        else:

            MaintenanceEvent(self.event_manager, self.event_manager.system.system_properties.summary_timestep*10) # FIXME: timestep!
            BrownianStartEvent(self.event_manager, 0)
            RotationalBrownianStartEvent(self.event_manager,0)

            if self.event_manager.system.system_properties.swelling_rate is not None:
                SwellingStartEvent(self.event_manager,0)

#         print("Going to SummaryEvent with a delay of ", self.event_manager.system.system_properties.summary_timestep)
        # Note: procession scheduling is handled in the constructor
        SummaryEvent(self.event_manager, self.event_manager.system.system_properties.summary_timestep) # turnoff
        



class NewtonianStartEvent(SystemEvent):
    def process(self, time):
        """
        1. only take non pinned particles
        2. distribute velocities and correct mean
        3. correct variance
        """
        yield self.event_manager.env.timeout(time)
        self.draw_velocities_and_correct_variances() #!Uncomment this
        self.predict_future_events(self.event_manager.system._particles) # predict events for all particles
        NewtonianEvent(self.event_manager, 0.001) # FIXME: add proper time difference, this should be zero? : 0.01


# FIXME: NOT NECESSARY IF CT-ALGORITHM WORKS RELIABLY
class NewtonianEvent(SystemEvent):
    """
    Newtonian Event: this event is a helper event and should not be used if the collisison time algorithm works reliably for long times!
    """
    def process(self, time):
        """
        1. Update all particles
        2. Predict events for each particle
        """
        yield self.event_manager.env.timeout(time)
        #print("inside Newtonian Event after yielding", self.event_manager.env.now,self.event_manager.system._particles[2].position)
        #self.event_manager.system.save_to_file() # placed by Aiyin
        [particle.update() for particle in self.event_manager.system._particles] # 1.1

        self.predict_future_events(self.event_manager.system._particles)

        # add new newtonian event with same timeout
        self.event_manager.env.process(self.process(time)) # NOTE: time is relative in simpy


class BrownianStartEvent(SystemEvent):
    def process(self, time):
        yield self.event_manager.env.timeout(time)
#         print("inside BrownianStartEvent", self.event_manager.env.now)
        self.draw_translational_velocities()
        self.event_manager.system.system_properties.brownian_timestep_number += 1 # i think this is just a counter
        
        # 4. repredict multi-particle events for all particles
        self.predict_future_events(self.event_manager.system._particles, ignore_last_collision_time=True)

        BrownianEvent(self.event_manager, self.event_manager.system.system_properties.brownian_timestep)
#         NewtonianEvent(self.event_manager, 0.0001)


class BrownianEvent(SystemEvent):
    def __deactivate_all_particleevents(self): ### FIXME: maybe rebuild the event tree instead
        """
        Deactivate all events of ParticleEvent type
        """
        for particle in self.event_manager.system._particles:
            if self.event_manager.events_per_particle != {}: # catch the case that all particles are not moving and no particle events were predicted
                for partevent in self.event_manager.events_per_particle[particle]:
                    partevent.active = False
            self.event_manager.events_per_particle[particle] = [] # remove obsolete events from events_per_particle lists to simplify list access

        
    def process(self, time):
        """
        1. Deactivate all future ParticleEvents
        2. Update all particles
        3. Redistribute velocities
        4. Predict events for each particle (cell uppper right in BDsim2d)
        """
        yield self.event_manager.env.timeout(time)
        self.__deactivate_all_particleevents() # 1.
        [particle.update() for particle in self.event_manager.system._particles] # 2.

        # 3.
#         print("inside BROWNIANEvent", self.event_manager.env.now)
        self.draw_translational_velocities()
        self.event_manager.system.system_properties.brownian_timestep_number += 1

        # 4. repredict multi-particle events for all particles
        self.predict_future_events(self.event_manager.system._particles, ignore_last_collision_time=True)
#         print("Apply CorrectionEvent after an overlap")
        # 5. add new brownian event with same timeout
        self.event_manager.env.process(self.process(time)) # NOTE: time is relative in simpy


class RotationalBrownianStartEvent(SystemEvent):
    def process(self, time):
        yield self.event_manager.env.timeout(time)
        #         print("inside BrownianStartEvent", self.event_manager.env.now)
        self.draw_rotational_velocities()
        self.event_manager.system.system_properties.brownian_timestep_number += 1  # i think this is just a counter

        # 4. repredict multi-particle events for all particles
        self.predict_future_events(self.event_manager.system._particles, ignore_last_collision_time=True)

        a2 = self.event_manager.system._particles[0].major**2. # note that this is assuming you have a monodisperse system
        RotationalBrownianEvent(self.event_manager,  3.*a2/10.*self.event_manager.system.system_properties.brownian_timestep)


# NewtonianEvent(self.event_manager, 0.0001)


class RotationalBrownianEvent(SystemEvent):
    def __deactivate_all_particleevents(self):  ### FIXME: maybe rebuild the event tree instead
        """
        Deactivate all events of ParticleEvent type
        """
        for particle in self.event_manager.system._particles:
            if self.event_manager.events_per_particle != {}:  # catch the case that all particles are not moving and no particle events were predicted
                for partevent in self.event_manager.events_per_particle[particle]:
                    partevent.active = False
            self.event_manager.events_per_particle[
                particle] = []  # remove obsolete events from events_per_particle lists to simplify list access

    def process(self, time):
        """
        1. Deactivate all future ParticleEvents
        2. Update all particles
        3. Redistribute velocities
        4. Predict events for each particle (cell uppper right in BDsim2d)
        """
        yield self.event_manager.env.timeout(time)
        self.__deactivate_all_particleevents()  # 1.
        [particle.update() for particle in self.event_manager.system._particles]  # 2.

        # 3.
        #         print("inside BROWNIANEvent", self.event_manager.env.now)
        self.draw_rotational_velocities()
        self.event_manager.system.system_properties.brownian_timestep_number += 1

        # 4. repredict multi-particle events for all particles
        self.predict_future_events(self.event_manager.system._particles, ignore_last_collision_time=True)
        #         print("Apply CorrectionEvent after an overlap")
        # 5. add new brownian event with same timeout
        self.event_manager.env.process(self.process(time))  # NOTE: time is relative in simpy
#
# #LINEAR
# class SummaryEvent(SystemEvent):
#     def process(self, time):
#         """
#         1. Update all particles
#         2. Write particle properties to file
#         """
#         yield self.event_manager.env.timeout(time)
#         [particle.update() for particle in self.event_manager.system._particles] # 1.
#         self.event_manager.system.save_to_file() # 2. - I commented this, this should be uncommented -A
#         self.event_manager.system.system_properties.summary_timestep_number += 1
#
#         # add new summary event with same timeout
#
#         self.event_manager.env.process(self.process(time)) # NOTE: time is relative in simpy

#LOG
class SummaryEvent(SystemEvent):
    def process(self, time):
        """
        1. Update all particles
        2. Write particle properties to file
        """
        yield self.event_manager.env.timeout(time)
        [particle.update() for particle in self.event_manager.system._particles] # 1.
        #logarithmic scale with initial value if summary_timestep and base of 1.03
        tscale = self.event_manager.system.system_properties.summary_timestep*1.02**self.event_manager.system.system_properties.summary_timestep_number
        now = self.event_manager.env.now
        time = tscale - now

        self.event_manager.system.save_to_file() # 2.
        self.event_manager.system.system_properties.summary_timestep_number += 1

        # add new summary event with same timeout
        #print(now)
        self.event_manager.env.process(self.process(time)) # NOTE: time is relative in simpy

# ##COLLISION
# class SummaryEvent(SystemEvent):
#     def process(self, time):
#         """
#         1. Update all particles
#         2. Write particle properties to file
#         """
#         yield self.event_manager.env.timeout(time)
#         [particle.update() for particle in self.event_manager.system._particles] # 1.
#         #logarithmic scale with initial value if summary_timestep and base of 1.03
#         tscale = self.event_manager.system.system_properties.summary_timestep*1.02**self.event_manager.system.system_properties.summary_timestep_number
#         now = self.event_manager.env.now
#         time = tscale - now
#
#         self.event_manager.system.save_to_file() # 2.
#         self.event_manager.system.system_properties.summary_timestep_number += 1
#
#         # add new summary event with same timeout
#         #print(now)
#         self.event_manager.env.process(self.process(time)) # NOTE: time is relative in simpy



class MaintenanceEvent(SystemEvent):
    """
    General maintenance event. Could be used to clean the tree or reassign cells to avoid the cumulation of numerical errors.
    """
    def process(self, time):
        yield self.event_manager.env.timeout(time)
        [particle.update() for particle in self.event_manager.system._particles] # 1.
        [particle.cell.cellspace.reassign_particle(particle) for particle in self.event_manager.system._particles]
        
        # add new maintenance event with same timeout
        self.event_manager.env.process(self.process(time)) # NOTE: time is relative in simpy


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
