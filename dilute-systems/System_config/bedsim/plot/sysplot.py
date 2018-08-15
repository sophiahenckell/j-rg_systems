import math

import numpy as np
import matplotlib
#matplotlib.use("Agg")

from pylab import figure, show
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
#import matplotlib.animation as manimation


# plot.draw, set_data
class Plot(object):

    def plot_cellspace(self):
        # Plot cells
        self.ax.triplot(self._system.cellspace._grid_points[:,0], self._system.cellspace._grid_points[:,1], self._system.cellspace.triang.simplices.copy())
        #for xy in list(self._system.cellspace._grid_points):
        #    xy2 = tuple(list(xy))
        #    self.ax.annotate('(%s, %s)' % xy2, xy=xy2, textcoords='offset points')


    def plot_ellipses(self): # static initialization
        # Plot particles
        self.ells = [Ellipse(xy=e.position, width=2*e.major, height=2*e.minor, angle=math.degrees(e.angle)) for e in self._system._particles]
        for e in self.ells:
            self.ax.add_artist(e)
            #e.set_clip_box(self.ax.bbox)
            #e.set_alpha(0.1)

            k = e.height/e.width
            if k>0.5:
                r=1
                b=0
            else:
                r=0
                b=1
            e.set_facecolor([r,0,b]) # FIXME: encode angle or particle type

    def update_system(self):
        # adjust ellipses
        for (ep, es) in zip(self.ells, self._system._particles): # FIXME: self.ells = [] after add_artist(e) for e in ells
            ep.center = es.position
            ep.angle = math.degrees(es.angle)


    def plot_system(self, outfile):
        self.update_system()

        #ax.set_ylim(-1,3)
        #ax.set_xlim(-1,3)

        #self.fig.canvas.draw() # apply ellipse changes
        plt.savefig(outfile)
        #plt.close()

    
    
    def plot_system_fmantrag(self, outfile):
        # Plot particles
        ells = [Ellipse(xy=e.position, width=2*e.major, height=2*e.minor, angle=math.degrees(e.angle)) for e in self._system._particles]
        for e in ells:
            self.ax.add_artist(e)
            e.set_clip_box(self.ax.bbox)
            #e.set_alpha(0.1)

            k = e.height/e.width
            if k>0.5:
                g=0.54902
                b=0
            else:
                g=0
                b=1
            e.set_facecolor([0,g,b]) # FIXME: encode angle or particle type
         
        # add circles to check area calculation for packing fraction
        """   
        c1 = plt.Circle((1,1), 0.5, color='r', fill=False)
        c2 = plt.Circle((1,1), 0.8, color='b', fill=False)
        ax.add_artist(c1)
        ax.add_artist(c2)
        """

        self.ax.set_ylim(-0.1,2.1)
        self.ax.set_xlim(-0.1,2.1)
        #plt.grid()
        plt.axis('off')
        plt.savefig(outfile)


    def system_sync(self):
        """
        Initial system plot.
        Call this function only when the system has fully loaded!
        """
        self.plot_cellspace()
        self.plot_ellipses()
        

    def __init__(self, system):
        self._system = system
        
        self.fig = plt.figure()
        plt.grid()
        
        self.ax = self.fig.add_subplot(111)
        self.ax.set_aspect(1)
        
        self.ells = []