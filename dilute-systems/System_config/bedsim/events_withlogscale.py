"""
Created on 04.02.2015

@author: mahe
"""

from collections import namedtuple

import sys #-A
import simpy
import numpy as np
import time

from bedsim.particle import ParticleOverlap, ConvergenceError, NotyetOverlapping


"""
# Abstract Event class.
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
    count_collision_predictions = 0 # FIXME: remove or leave for stats
    print("inside Event class checked List and repredicted list initialized")

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
                                
    def overlap_check2(self, particles):
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
                            Fab = particle.overlap2(neigh_particle, 0)
                            if Fab < -1.1e-4:
                                print("Particles ", particle.get_name(), " and ", neigh_particle.get_name(), " overlap!", Fab)


    # FIXME: use until instead of brownian step!
    # FIXME: DOES NOT CONSIDER DIFFERENT UPDATE STATES OF THE PARTICLES!!! This is circumvented atm by updating all particles when the state changes
    def predict_future_events(self, particles, ignore_last_collision_time=False):
        """
        Predict and schedule all future events for each particle of list particles.
        @param particles List of particles to predict and schedule particle events for
        @param ct, either a type float (gives the next prediction time) or type np.ndarray (gives the next collision time with point of collision and normal vector)
        """

#         tb = self.event_manager.system.system_properties.brownian_timestep
#         if tb is None:
#             until = 0.01 # FIXME: depends on packing fraction!
#         else:
#             now = self.event_manager.system.system_properties.localtime() % tb
#             until = 0.01 # FIXME: brownian_timestep - (current_time % brownian_timestep)

        until = 1e-2
        
        done = {}
        epsilon = 1e-8
        
        for particle in particles:
            """ 1. predict cell crossings """
            cc = particle.cell.crossing_time(particle.position, particle.velocity)
            if cc is not None:
                (crossing_time, to_cell , crossing_point) = cc
                if until >= crossing_time > epsilon: # this is to remove negative times in cell crossing predictions check thi
                    CellCrossingEvent(self.event_manager, crossing_time, [particle], to_cell, crossing_point)

            """ 2. predict collisions """
            pid1 = int(particle._id)
            # check the 26 neighbouring cells surrounding  and the cell, where the particle itself belong (27 total)
            for neigh_cell in [particle.cell] + particle.cell.neighbours: # check for collisions of particle in its own cell + neighbour cells
                
                # for each particle in each of the neighbouring cells, this can reach to 20 (ave) for dense ellipsoidal systems
                for neigh_particle in neigh_cell.particles:
                    # minimize the neighbor cells further
                    if self.event_manager.system.boundary.delta_norm(particle.position,neigh_particle.position) <= (particle.major + neigh_particle.major + 0.1):

                        pid2 = int(neigh_particle._id)
                        if pid1 < pid2:
                            key = (pid1, pid2)
                        else:
                            key = (pid2, pid1)
                            
                        if particle is not neigh_particle and key not in done: # don't check particles twice
                            Event.count_collision_predictions += 1
                            done[key] = True

                                
                            try:
                                ct = particle.tpw_collision_time3(neigh_particle, T=until)

                                if ct is not None: # if collision time algorithm has no upper limit, also compare to Brownian time step if set and discard event if later than Brownian time step

                                    if isinstance(ct, float):
#                                         print("here 1",key,self.event_manager.system.system_properties.localtime())
                                        if neigh_particle in particle.lastCollision: # check if the two particles just collided, these particles are highly likely to collide again
                                            PredictionEvent(self.event_manager,ct, [particle, neigh_particle]) #1e-4 is the one used inside the tpw_collion time so it has to be bigger than this

                                    else:
                                        CollisionEvent(self.event_manager, ct[0], [particle, neigh_particle], ct[1:8])
                                        
#                                 except NotyetOverlapping as ne:
#                                     """
#                                     Take a second look at these particles
#                                     make sure that future overlaps are repredicted in the future i.e. after 2e-4
#                                     val2 : obtained from tpw_collision_time3, time needed to call for PredictionEvent for that pair
#                                     Event.checkedList: list of particle pairs and the last time they were predicted
#                                     """
#                                     (val2,) = ne.args
#                                               
#                                     if neigh_particle in particle.lastCollision: # check if the two particles just collided, these particles are highly likely to collide again
#                                         if key in Event.checkedList:
#                                             predictingtime = self.event_manager.system.system_properties.localtime() - Event.checkedList[key] - val2
#             
#                                             if predictingtime > 0.01: # if the last prediction time was too long ago (in this case it means > 0.01), predict again, in a volume fraction of 0.5 I saw that some particles needed reevaluation after 0.15                
#             
#                                                 Event.checkedList[key] = self.event_manager.system.system_properties.localtime()
#                                                 PredictionEvent(self.event_manager,val2, [particle, neigh_particle]) #2e-4 is the one used inside the tpw_collion time so it has to be bigger than this
#                                         else:
#             
#                                             Event.checkedList[key] = self.event_manager.system.system_properties.localtime()
#                                             PredictionEvent(self.event_manager,val2, [particle, neigh_particle])
#             
                            except ConvergenceError as ce:
                                ((val,),) = ce.args 
                                if val is not None and val > 1e-12:
                                    PredictionEvent(self.event_manager, val, [particle, neigh_particle])
                                        
                                        
    def predict_pair(self, particles, ignore_last_collision_time=False):
        """
        For particle pairs that are likely to collide, check if they are colliding
        """

#         tb = self.event_manager.system.system_properties.brownian_timestep
#         if tb is None:
#             until = 0.01 # FIXME: depends on packing fraction!
#         else:
#             now = self.event_manager.system.system_properties.localtime() % tb
#             until = 0.01 # FIXME: brownian_timestep - (current_time % brownian_timestep)

        until = 1e-2
        epsilon = 1e-8
        
        [particle, neigh_particle] = particles
        
        
        key = (particle._id, neigh_particle._id)
        Event.count_collision_predictions += 1
            
        try:
            ct = particle.tpw_collision_time3(neigh_particle, T=until)
            
            if ct is not None:
                if isinstance(ct, float):
#                     print("here 2",key)
                    PredictionEvent(self.event_manager,ct, [particle, neigh_particle])
                else:
                    CollisionEvent(self.event_manager, ct[0], [particle, neigh_particle], ct[1:8])

        except ConvergenceError as ce:
            ((val,),) = ce.args 
            if val is not None and val > 1e-12:
                PredictionEvent(self.event_manager, val, [particle, neigh_particle])

                                
    def goto_max_volume(self, particles, ignore_last_collision_time=False):
        """
        Predict the maximum volume the system can reach at this event
        @param particles List of particles to predict and schedule particle events for
        """
 
         
        swellrate = self.event_manager.system.system_properties.swelling_rate
 
        until = 0.1 # FIXME: just find a range of time where the swelling can be predicted
         
        done = {}
        epsilon = 1e-8
         
        ct0 = [ ]
        n = len(particles)
        for particle in particles:
 
            """ 2. predict collisions """            
            for neigh_cell in [particle.cell] + particle.cell.neighbours: # check for collisions of particle in its own cell + neighbour cells
                 
                for neigh_particle in neigh_cell.particles:
                    
                    try:
                        ct = particle.find_swelling_time(neigh_particle, T=until)            
                        if ct is not None:
                            ct0.append(ct[0])
                            print(particle._id,neigh_particle._id,ct[0],'inside goto_max volume')
   
                    except ParticleSwellOverlap: # this is actually a ParticleOverlap error inheriting a ValueError
                        """ catch overlap """
                        print("inside predict_future_events catch particle overlap, maybe contract the particles a little",particle._id,neigh_particle._id)
                        ct0.append(0.)
                        CorrectionEvent(self.event_manager, 0.0, [particle, neigh_particle])
                        sys.exit("Swell problem")
 
                    except ConvergenceError as ce:
                        print("inside predict_max_volume,ConvergenceError")
                        ((val,),) = ce.args # FIXME: proper Exception class
                        if val is not None and val > 1e-12:
                            #print("\n>>> VLA=",val)
                            PredictionEvent(self.event_manager, val, [particle, neigh_particle]) 
        
        

        system_volume = self.event_manager.system.compute_volume_fraction()
        required_volume = self.event_manager.system.system_properties.packing_fraction
        k = self.event_manager.system.system_properties.aspect_ratio
         
        if system_volume < required_volume:
            major_temp = particles[0].major
         
            for particle in particles:
                 
                if ct0 != []:

                    st = min(ct0) - 1e-3
                    particle.major += swellrate * st
                    particle.minor = particle.major/k
                    particle.minor2 = particle.minor
             
            current_volume = self.event_manager.system.compute_volume_fraction()
             
            # if the volume exceeds the required volume fraction, reduce it to the required volume
            if current_volume > required_volume:
                change_of_volume = current_volume - required_volume
                print("change of volume", change_of_volume)
                box_volume = self.event_manager.system.box_volume()
                del_major = (3.* change_of_volume * k * k * box_volume / (4. * n * np.pi))**(1./3.)
                print("del_major",del_major)
                 
                for particle in particles:
                    particle.major = (particle.major**3 - del_major**3)**(1./3.)
                    particle.minor = particle.major/k
                    particle.minor2 = particle.minor
                     
                current_volume = self.event_manager.system.compute_volume_fraction()
                self.volume_reached = True
                print("current volume after reduction",current_volume)
                 
            print("SYSTEM VOLUME",current_volume)                    
                                
    

#     def goto_max_volume(self, particles, ignore_last_collision_time=False):
#         """
#         Predict the maximum volume the system can reach at this event
#         @param particles List of particles to predict and schedule particle events for
#         """
#  
#          
#         swellrate = self.event_manager.system.system_properties.swelling_rate
#  
#         until = 0.1 # FIXME: just find a range of time where the swelling can be predicted
#          
#         done = {}
#         epsilon = 1e-8
#          
#         ct0 = [ ]
#         n = len(particles)
#         for particle in particles:
#  
#             """ 2. predict collisions """            
#             for neigh_cell in [particle.cell] + particle.cell.neighbours: # check for collisions of particle in its own cell + neighbour cells
#                  
#                 for neigh_particle in neigh_cell.particles:
#                     try:
#                         ct = particle.find_swelling_time(neigh_particle, T=until)
#                                                      
#                         if ct is not None:
#                             ct0.append(ct[0])
#    
#                     except ParticleOverlap: # this is actually a ParticleOverlap error inheriting a ValueError
#                         """ catch overlap """
#                         print("inside predict_future_events catch particle overlap, maybe contract the particles a little",particle._id,neigh_particle._id)
#                         ct0.append(0.)
#  
#                     except ConvergenceError as ce:
#                         print("inside predict_max_volume,ConvergenceError")
#                         ((val,),) = ce.args # FIXME: proper Exception class
#                         if val is not None and val > 1e-12:
#                             #print("\n>>> VLA=",val)
#                             PredictionEvent(self.event_manager, val, [particle, neigh_particle]) 
#          
#         system_volume = self.event_manager.system.compute_volume_fraction()
#         required_volume = self.event_manager.system.system_properties.packing_fraction
#         k = self.event_manager.system.system_properties.aspect_ratio
#          
#          
#         if system_volume < required_volume:
#             major_temp = particles[0].major
#          
#             for particle in particles:
#                  
#                 if ct0 != []:
#                     st = min(ct0) - 1e-6
#                     particle.major += swellrate * st
#                     particle.minor = particle.major/k
#                     particle.minor2 = particle.minor
#              
#             current_volume = self.event_manager.system.compute_volume_fraction()
#              
#             # if the volume exceeds the required volume fraction, reduce it to the required volume
#             if current_volume > required_volume:
#                 change_of_volume = current_volume - required_volume
#                 print("change of volume", change_of_volume)
#                 box_volume = self.event_manager.system.box_volume()
#                 print(n)
#                 del_major = (3.* change_of_volume * k * k * box_volume / (4. * n * np.pi))**(1./3.)
#                 print("del_major",del_major)
#                  
#                 for particle in particles:
#                     particle.major = (particle.major**3 - del_major**3)**(1./3.)
#                     particle.minor = particle.major/k
#                     particle.minor2 = particle.minor
#                      
#                 current_volume = self.event_manager.system.compute_volume_fraction()
#                 self.volume_reached = True
#                 print("current volume after reduction",current_volume)
#                  
#             print("SYSTEM VOLUME",current_volume)                    


    def __init__(self, event_manager, time):
        """ Save EventManager reference and call abstract process method. """
        self.event_manager = event_manager
        self.proc = self.event_manager.env.process(self.process(time))
        self.volume_reached = False


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
            #p.velocity = r[:2] / p.mass() # returns a numpy.array()
            #p.angvel   = r[-1] /(p.mass() * (p.major**2 + p.minor**2)/4) # returns a number # FIXME: schaue ob Trägheitstensor ok ist!!

        # 3. correct variance of velocities
        velvar = np.array([np.dot(r[:3], r[:3]) for r in rn])
        velvar = np.sum(velvar) / len(particles)
        angvelvar = np.array([np.dot(r[-3:], r[-3:]) for r in rn])
        angvelvar = np.sum(angvelvar) / len(particles)
        for p in particles:
            
            p.velocity *= np.sqrt(1 / velvar) # FIXME: mass
            p.angvel   *= np.sqrt(5 / 2.* p.minor / angvelvar) # FIXME: mass
            
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
#             self.event_manager.system.save_to_file() # placed by Aiyin
            
            self.particles[0].collision(self.particles[1], self.collision_point) # perform collision calculations
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
    Event for manually scheduling a prediction of future events 
    """
    def process(self, time):
        yield self.event_manager.env.timeout(time)
#         print("Prediction event started at ", self.event_manager.env.now)
        """ the following code is executed after timeout """
        if self.active:
            """ update and check for future events """
            [p.update() for p in self.event_manager.system._particles]
            self.predict_pair(self.particles)
            
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
            [p.update() for p in self.event_manager.system._particles]
            for p in self.particles: # self.particles is a 1 element list, but this looks nicer
                    p.update() # update is done in cell_cross method => is needed here for cumulative position ### FIXME: Position sollte NUR in update-Methode verändert werden!

                    p.cell_cross(self.to_cell, self.crossing_point)
            
            self.predict_future_events(self.particles)
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
        print("for quit event saving")
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

        self.event_manager.system.save_to_file() # 2.

        if self.event_manager.system.system_properties.brownian_timestep is None:
            print("doing Newtonian Dynamics")
            print("going to MaintenanceEvent with a delay of", self.event_manager.system.system_properties.summary_timestep)
            MaintenanceEvent(self.event_manager,1e-2) # FIXME: timestep! #self.event_manager.system.system_properties.summary_timestep
            print("going to NewtonianStartEvent event with no delay")
            NewtonianStartEvent(self.event_manager, 0)
            
        else:
            print("going to MaintenanceEvent for BD with a delay of",self.event_manager.system.system_properties.summary_timestep*10)
            MaintenanceEvent(self.event_manager, 1e-2)  #self.event_manager.system.system_properties.brownian_timestep/10.
            print("going to BrownianStartEvent with no delay")
#             SystemPredictionEvent(self.event_manager,0.01)
            BrownianStartEvent(self.event_manager, 0)
            print("going to SwellingEvent with no delay")

        print("Going to SummaryEvent with a delay of ", self.event_manager.system.system_properties.summary_timestep)
        SummaryEvent(self.event_manager, self.event_manager.system.system_properties.summary_timestep)
        



class NewtonianStartEvent(SystemEvent):
    def process(self, time):
        """
        1. only take non pinned particles
        2. distribute velocities and correct mean
        3. correct variance
        """
        yield self.event_manager.env.timeout(time)
        self.draw_velocities_and_correct_variances()
        self.predict_future_events(self.event_manager.system._particles) # predict events for all particles
        NewtonianEvent(self.event_manager, 0.1) # FIXME: add proper time difference, this should be zero? : 0.01


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
        #print("inside Newtonian Event after yielding", self.event_manager.env.now)
#         self.event_manager.system.save_to_file() # placed by Aiyin
        [particle.update() for particle in self.event_manager.system._particles] # 1.1
        
        self.predict_future_events(self.event_manager.system._particles)
#         self.event_manager.system.save_to_file()
        # add new newtonian event with same timeout
        self.event_manager.env.process(self.process(time)) # NOTE: time is relative in simpy


class BrownianStartEvent(SystemEvent):
    def process(self, time):
        yield self.event_manager.env.timeout(time)
#         print("inside BrownianStartEvent", self.event_manager.env.now)
        self.draw_velocities_and_correct_variances()
        self.event_manager.system.system_properties.brownian_timestep_number += 1 # i think this is just a counter
        
        # 4. repredict multi-particle events for all particles
        self.predict_future_events(self.event_manager.system._particles, ignore_last_collision_time=True)

        BrownianEvent(self.event_manager, self.event_manager.system.system_properties.brownian_timestep)
#         self.event_manager.system.save_to_file()
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
        Event.checkedList = {} # clear checkedList too
        
        [particle.update() for particle in self.event_manager.system._particles] # 2.
        
        self.draw_velocities_and_correct_variances()
        self.event_manager.system.system_properties.brownian_timestep_number += 1
        # 4. repredict multi-particle events for all particles
        self.predict_future_events(self.event_manager.system._particles, ignore_last_collision_time=True)
        
        self.event_manager.env.process(self.process(time)) # NOTE: time is relative in simpy
        
class SwellingStartEvent(SystemEvent):
    def process(self, time):
        yield self.event_manager.env.timeout(time)
        self.set_velocities_to_zero_for_swelling()

        self.goto_max_volume(self.event_manager.system._particles, ignore_last_collision_time=True)
        print("save first swell event to file")
        self.event_manager.system.save_to_file() # 2. - I commented this, this should be uncommented -A
#         self.event_manager.system.save_to_file() # 2. - I commented this, this should be uncommented -A
        self.draw_velocities_and_correct_variances()        
        self.predict_future_events(self.event_manager.system._particles, ignore_last_collision_time=True)
        SwellingEvent(self.event_manager, self.event_manager.system.system_properties.brownian_timestep*2.)


class SwellingEvent(SystemEvent):
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
        3. Set velocities and angular velocities to zero
        4. Predict events for each particle 
        5. Assign new velocities and angular velocities
        """
        yield self.event_manager.env.timeout(time)
        self.__deactivate_all_particleevents() 
        
        [particle.update() for particle in self.event_manager.system._particles] # 2.
        
        if self.overlap_check(self.event_manager.system._particles) is None:
        
            self.set_velocities_to_zero_for_swelling()
            self.goto_max_volume(self.event_manager.system._particles, ignore_last_collision_time=True)
            print("save swelling event to file")
            self.event_manager.system.save_to_file() # 2. - I commented this, this should be uncommented -A
            self.draw_velocities_and_correct_variances()
            self.predict_future_events(self.event_manager.system._particles, ignore_last_collision_time=True)
            self.event_manager.env.process(self.process(time)) # NOTE: time is relative in simpy
        else:
            print("THERE'S SWELLING OVERLAP")
            self.event_manager.env.process(self.process(time+1e-4))


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
        self.event_manager.env.process(self.process(time)) # NOTE: time is relative in simpy


class MaintenanceEvent(SystemEvent):
    """
    General maintenance event. Could be used to clean the tree or reassign cells to avoid the cumulation of numerical errors.
    """
    def process(self, time):
        yield self.event_manager.env.timeout(time)

        [particle.update() for particle in self.event_manager.system._particles] # 1.
#         [particle.cell.cellspace.reassign_particle(particle) for particle in self.event_manager.system._particles]

        self.predict_future_events(self.event_manager.system._particles)
#         # add new maintenance event with same timeout
        self.event_manager.env.process(self.process(time)) # NOTE: time is relative in simpy
        
class SystemPredictionEvent(SystemEvent):
    """
    General maintenance event. Could be used to clean the tree or reassign cells to avoid the cumulation of numerical errors.
    """
    def process(self, time):
        yield self.event_manager.env.timeout(time)
        # add the next two lines or not, probably an overkill
#         [particle.update() for particle in self.event_manager.system._particles] # 1.
#         [particle.cell.cellspace.reassign_particle(particle) for particle in self.event_manager.system._particles]
#         self.overlap_check2(self.event_manager.system._particles)
        self.predict_future_events(self.event_manager.system._particles)
        
        self.event_manager.env.process(self.process(time)) # NOTE: time is relative in simpy


"""
# Event Manager
Event Manager should be a special implementation of an abstract base class for general
dynamics simulations.
"""
class EventManager(object):
    def start(self):
        print("StartEvent",self.env.now)
        StartEvent(self, 0)
        print("after StartEvent",self.system.system_properties.lifetime - self.system.system_properties.summary_timestep/2)
        print(self.system.system_properties.lifetime,self.system.system_properties.summary_timestep)
        QuitEvent(self, self.system.system_properties.lifetime - self.system.system_properties.summary_timestep/2) # NOTE: simpy simulates from [tstart, tend)
        print("after QuitEvent",self.env.now)
        self.env.run(until=self.system.system_properties.lifetime)

    def __init__(self, system):
        self.system = system # Reference to the system
        self.env = simpy.Environment() # create simpy environment
        self.events_per_particle = {}  # create per particle event lists
