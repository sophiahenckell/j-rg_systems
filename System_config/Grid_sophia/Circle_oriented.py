#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 14 13:36:06 2018
Eine Klasse, die der Ellipse(Orientable) ähnlich sein soll, nur für Kugeln
@author: sophia
"""

class Circle_Sophia(Orientable):
    
    def mass(self):
        return np.float64(1) #np.pi * self.major * self.minor 
    # FIXME: only valid for constant density 
    # FIXME: proportionality constant...
    
    def update(self): #????
        deltaT = Orientable.update(self)
        
        return deltaT
    
    def collision(self, with_circle, collision_point): ##?? ellipse

        self.lastCollision[with_circle] = self.cell.cellspace.system.system_properties.localtime()
        with_circle.lastCollision[self] = self.lastCollision[with_circle] ##???

        sigma = self.radius + with_circle.radius # should be the same as norm(delta_r)
        delta_r = self.cell.cellspace.system.boundary.delta_dir(with_circle.position, self.position) # r1-r2
        delta_v = self.velocity - with_circle.velocity # v1-v2
        delta_rv = np.dot(delta_r, delta_v)
        
        #self.velocity = self.velocity - (delta_r * delta_rv/sigma**2) * 2/(1/self.mass() + 1/with_circle.mass())
        #with_circle.velocity = with_circle.velocity + (delta_r * delta_rv/sigma**2) * 2/(1/self.mass() + 1/with_circle.mass())
        
        self.velocity = self.velocity - 2 * delta_rv /((1+self.mass()/with_circle.mass()) * sigma) * delta_r/sigma
        with_circle.velocity = with_circle.velocity + 2 * delta_rv /((1+with_circle.mass()/self.mass()) * sigma) * delta_r/sigma
    
    def collision_for_swelling(self, withEllipse, collisionPoint): ##??? brauch ich nicht oder?
        return None
    
    def correction(self, withEllipse): #TODO 
        return None
    
    def overlap(self, with_circle):
        delta = self.cell.cellspace.system.boundary.delta_dir(self.position, with_circle.position)
        #return np.linalg.norm(delta) < self.radius + with_circle.radius
        return abs(np.linalg.norm(delta) - self.radius + with_circle.radius) < 10**-8
    
    