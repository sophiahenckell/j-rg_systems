#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 25 16:22:35 2018

@author: sophia
"""
import numpy as np

def norm(x,y,z):
    k = np.array([x,y,z])
    l = np.linalg.norm(k)
    return l
    return 0*l
    

test = norm(1,2,3)