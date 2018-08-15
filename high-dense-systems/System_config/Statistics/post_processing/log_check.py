#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 15 13:38:59 2018

@author: sophia
"""

import numpy as np
import matplotlib.pyplot as plt

k = np.array([0, 1, 1.03, 1.0609, 1.1255, 57.37, 77.103])
i = np.array( [0,1,2,3,4,5,6])

fig, ax = plt.subplots()
#ax.set_xscale('log', basex = 1.03)
ax.set_yscale('log', basex = 1.03)
ax.plot(k)
plt.show()


