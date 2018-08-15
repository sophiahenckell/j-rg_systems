#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  9 15:35:31 2018

@author: sophia
"""

import simpy
def clock(env, name, tick):
    while True:
        print(name, env.now)
        yield env.timeout(tick)
        
env = simpy.Environment()
env.process(clock(env, 'fast', 0.05))
env.process(clock(env,'slow', 1))
env.run(until = 3)

