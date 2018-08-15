#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 12 12:02:48 2018
Dieses Programm nimmt die ausgewerteten h5 dateien
(ausgewertet über: python3 hdf2stat.py --filename data/Ellipse....h5)
und stellt MSD und MSR grafisch dar. 
@author: sophia
"""

import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('/data/scc/sophia/MSDCalculations1-n256-t4000-phi0.4-k1.2-tb0.1-ts10/Ellipse-Benchmark-n256-k1.2-phi0.4-id1-equilibrate.h5-msd.dat')
vector = np.linspace(0,4000, len(data))
plt.plot(vector,data)


#data2 = np.loadtxt('data/Ellipse-Benchmark-n1000-k1.7-phi0.5-id1.h5-msr.dat')
#plt.loglog(vector, data2)

plt.title('Ellipse-Benchmark-n100-k0-phi0.2-id1.h5-msd.eps')
plt.show()

