from math import degrees
import numpy as np
#from pylab import figure, show
#from matplotlib.patches import Ellipse
import matplotlib.patches 
import matplotlib.pyplot as plt
from matplotlib.collections import EllipseCollection
#from matplotlib.pyplot import text
#import matplotlib.animation as animation

try:
    from moviepy.editor import VideoClip
    from moviepy.video.io.bindings import mplfig_to_npimage
    USE_MOVIEPY = True
except ImportError:
    print("moviepy has not been found. No video will be generated, only vec frames according to the settings.")
    USE_MOVIEPY = False

from bedsim.files import BedfileReader

from bedsim.cell import DelaunayCellspace


FPS = 30
USE_BOUNDARY_REDRAW = True
DRAW_CELLS = False
DRAW_PARTICLES = True
DRAW_PARTICLE_IDS = False
DRAW_TIME = True
DRAW_CELL_NEIGHBOURS = False
DRAW_FIG_AXES = True
FIGSIZE = (20,20)
SAVE_VEC_FRAMES = True # save a vector graphics frame after each "VEC_FRAME_STEP"
VEC_FRAME_STEP = 10


class arb(object):
    pass


class Animate(object):
    
    """
    Moviepy related functions (based on matplotlib routines below)
    """
    def moviepy_update(self, i):
        self.update(i)
        if SAVE_VEC_FRAMES:
            if i*FPS % VEC_FRAME_STEP == 0:
                self.save_vec_frame(i)
        return mplfig_to_npimage(self.fig) if USE_MOVIEPY else None
    
    def moviepy_animate(self):
        self.prepare()
        if USE_MOVIEPY:
            self.ani = VideoClip(self.moviepy_update, duration=self.bf.particles.get_number_of_frames()/FPS) # FIXME: duration in seconds
        else:
            for i in range(self.bf.particles.get_number_of_frames()):
                print(i)
                self.moviepy_update(i)
                
    def moviepy_save_to_file(self, filename):
        if USE_MOVIEPY:
            self.ani.write_videofile(filename, fps=FPS)
    
    
    
    def save_vec_frame(self, i):
        self.fig.savefig("%s-frame-ts%d.pdf" % (self.output_filename, i*30))

    
    """
    Matplotlib related functions
    """
    def update(self, i):
        """
        Update dynamic particle properties.
        """
        if DRAW_PARTICLES:
            particles = next(self.particle_data)
            if DRAW_TIME:
                # FIXME: evtl. in Simulationszeit konvertieren, also *saving_timestep
                #self.ax.set_title(r'$t_s$=%3.2f' % (i*FPS))
                #self.ax.set_title(r'$t_s$=%03d' % (i*FPS))
                self.ax.set_title(r'$t_s$=%d' % (i*FPS))
                #self.ax.set_title(r'$t_s$=%d' % (i*FPS), fontsize=80)
            
            if not USE_BOUNDARY_REDRAW:
                for (ep, es) in zip(self.ells, particles): # FIXME: self.ells = [] after add_artist(e) for e in ells
                    ep.center = es['position']
                    if not isinstance(ep, matplotlib.patches.Circle): # FIXME: das ist noch nicht so wirklich schön
                        ep.angle = degrees(es['angle'])
                    if DRAW_PARTICLE_IDS:
                        self.labels[es['id']].set_position(xy=es['position'])
    
            ### VERSION MIT WIEDERHOLDUNGEN AN RÄNDERN
            else:
                [box_width, box_height] = self.__box_size
                translations = [np.array([dx, dy]) for dx in [box_width, -box_width, 0] for dy in [box_height, -box_height, 0]]
                for (ep, es) in zip(self.ells2, particles): # FIXME: self.ells = [] after add_artist(e) for e in ells
                    ep.set_offsets( es['position'] + translations )
                    if (ep._heights != ep._widths).all():
                        ep._angles = np.repeat((es['angle']), 9)
                
                    if DRAW_PARTICLE_IDS:
                        self.labels[es['id']].set_position(xy=es['position'])
                        
            ### ENDE DER 2. Version

        


    def prepare(self):
        """
        Initialization of the particle patches.
        """
        self.plot_cellspace()
        #particles = next(self.particle_data)
        particles = next(self.particle_static_data)
        
        
        if not DRAW_FIG_AXES:
            self.ax.axis('off')
            
        
        if DRAW_PARTICLES:
            if not USE_BOUNDARY_REDRAW:
                #self.ells = [Ellipse(xy=e['position'], width=2*e['major'], height=2*e['minor'], angle=degrees(e['angle'])) for e in particles]
                #self.ells = [getattr(matplotlib.patches, e['type'])(xy=e['position'], width=2*e['major'], height=2*e['minor'], angle=degrees(e['angle'])) for e in particles]
                self.ells = []
                self.labels = {}
                for e in particles: # FIXME: das ist noch nicht so wirklich schön
                    if e['type'] == "Ellipse":
                        self.ells.append(matplotlib.patches.Ellipse(xy=e['position'], width=2*e['major'], height=2*e['minor'], angle=degrees(e['angle'])))
                    elif e['type'] == "Circle":
                        self.ells.append(matplotlib.patches.Circle(xy=e['position'], radius=e['radius']))
                    if DRAW_PARTICLE_IDS:
                        self.labels[e['id']] = self.ax.text(e['position'][0], e['position'][1], ("%d" % e['id']), fontsize=12, color='white')
                
                for e in self.ells:
                    self.ax.add_artist(e)
        
                    k = e.height/e.width
                    
                    # Antrag (R/G)
                    """"(r,g,b)=(0,0,1)
                    if k>0.5:
                        (g,b) = (0.54902, 0) 
                    """
                    # R/B
                    (r,g,b) = (0,1,0)
                    if k>0.5:
                        (r,b)=1,0
                    
                    e.set_facecolor([r,g,b]) # FIXME: encode angle or particle type
                
            else:
                ### VERSION MIT WIEDERHOLDUNGEN AN RÄNDERN
                self.ells2 = []
                [box_width, box_height] = self.__box_size
                translations = [np.array([dx, dy]) for dx in [box_width, -box_width, 0] for dy in [box_height, -box_height, 0]]
                for e in particles:
                    if e['type'] == "Ellipse":
                        # FIXME: this should be done in boundary module!
                        widths = np.repeat(2*e['major'], 9)
                        heights = np.repeat(2*e['minor'], 9)
                        angles = np.repeat(degrees(e['angle']), 9)
                        (r,g,b)=(0,0,1)
                        if e['minor']/e['major']>0.5:
                            (r,g,b)=(1,0,0) # probably not too useful, because pinned particles have same color
#                             (g,b)=(0.54902,0)
                    elif e['type'] == "Circle":
                        widths = np.repeat(2*e['radius'], 9)
                        heights = np.repeat(2*e['radius'], 9)
                        angles = np.repeat(0, 9)
                        (r,g,b)=(1,0,0)
                    if e['pinned'] == True:
                        (r,g,b)=(1,0,0)
                    if DRAW_PARTICLE_IDS:
                        self.labels[e['id']] = self.ax.text(e['position'][0], e['position'][1], ("%d" % e['id']), fontsize=12, color='white')
                        
                    XY = e['position'] + translations
                    ec = EllipseCollection(widths, heights, angles, units='x', offsets=XY, transOffset=self.ax.transData)
                    ec.set_facecolor([r,g,b])
                    self.ells2.append(ec)
                    self.ax.add_collection(ec)

        

    def plot_cellspace(self):
        """
        Plots the Delaunay triangulation of the cellspace.
        """
        system = arb()
        self.cellspace = DelaunayCellspace(system=system, grid_points=self.grid_data)
        
        if DRAW_CELLS:
            self.ax.triplot(self.cellspace._grid_points[:,0], self.cellspace._grid_points[:,1], self.cellspace.triang.simplices.copy())
        self.__grid_to_box_corners()
        (xmin, xmax, ymin, ymax) = self.__box_corners
        self.ax.set_xlim(xmin, xmax)
        self.ax.set_ylim(ymin, ymax)
        
        if DRAW_CELL_NEIGHBOURS:
            ### load boundaries as well to show appropriate neighbours
            system.cellspace = self.cellspace
            from bedsim.boundary import PeriodicBox
            self.boundary = PeriodicBox(system=system)
            system.boundary = self.boundary
            ###
            np.random.seed(3)
            for cell in self.cellspace.cells:
                #print("cs = ", cell._simplex)
                #if (cell._simplex == np.array([15, 21, 16])).all():
                #if (cell._simplex == np.array([15, 11, 10])).all():
                #if (cell._simplex == np.array([21, 17, 16])).all():
                if (cell._simplex == np.array([5, 1, 0])).all():
                    cell_points = self.cellspace._grid_points[cell._simplex]
                    mpoint = np.sum(cell_points, axis=0)/len(cell_points)
                    
                    for neigh in cell.neighbours:
                        neigh_cell_points = self.cellspace._grid_points[neigh._simplex]
                        neigh_mpoint = np.sum(neigh_cell_points, axis=0)/len(neigh_cell_points)
                        
                        self.ax.arrow(mpoint[0], mpoint[1], neigh_mpoint[0]-mpoint[0]+np.random.normal(0, 0.1), neigh_mpoint[1]-mpoint[1]+np.random.normal(0, 0.1), head_width=0.05, head_length=0.1, fc='r', ec='r')
                #self.ax.text(mpoint[0], mpoint[1], "I'm a cell")
                

    """ # OBSOLETE:
    def animate(self):
        self.ani = animation.FuncAnimation(self.fig, self.update, frames=self.bf.particles.get_number_of_frames(), interval=10, blit=True, init_func=self.prepare) # FIXME: frames

    def save_to_file(self, filename):
        self.ani.save(filename, fps=25, extra_args=['-vcodec', 'libx264'])
    """

    def load_from_file(self, filename):
        self.bf = BedfileReader(filename)
        self.particle_data = self.bf.particles.load_dynamics() # FIXME: im Moment ist noch alles dynamisch. Das wird sich aber bald ändern!
        self.particle_static_data = self.bf.particles.load_statics()
        self.grid_data = self.bf.system.grid
        
    def __init__(self, input_filename, output_filename):
        self.ells = []
        self.load_from_file(input_filename)
        self.fig, self.ax = plt.subplots(figsize=FIGSIZE)
        self.ax.set_aspect(1)
        self.ax.set_ylim(0,10) # FIXME
        self.ax.set_xlim(0,10) # FIXME
        self.labels = {}
        
        self.__box_size = np.array([10, 10])
        self.__box_center = np.array([5,5])
        self.__box_corners = (0, 10, 0, 10)
        
        # matplotlib
#         self.animate()
        #self.save_to_file(output_filename)
        self.output_filename = output_filename
        
        # moviepy
        self.moviepy_animate()
        self.moviepy_save_to_file(output_filename)
        
    
    def __grid_to_box_corners(self): ## FIXME: generalize this
        """
        Calculate the corners of the simulation box from the cell grid.
        """
        [x,y] = self.cellspace._grid_points.transpose()
        [xmin, xmax, ymin, ymax] = [np.amin(x), np.amax(x), np.amin(y), np.amax(y)]
        [width, height] = [xmax-xmin, ymax-ymin]
        self.__box_size = np.array([width, height])
        self.__box_center = np.array([(xmin+xmax)/2, (ymin+ymax)/2])
        self.__box_corners = (xmin, xmax, ymin, ymax)
        self.extent = np.fabs(np.array([xmax, ymax]) - np.array([xmin, ymin]))  
