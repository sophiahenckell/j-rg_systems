import numpy as np

from bedsim.files import BedfileReader


class arb(object):
    pass

class Statprop(object):
    
    def _uglyplot_msd(self, outfilename):
        import matplotlib.pyplot as plt
        plt.figure(figsize=(4,3))
         
        # FIXME: replace 10 with actual time, 0.1/index 5 with Brownian step/summary step
        x = np.linspace(start=0, stop=10000, num=len(self.msd)+1)[1:]
        print(len(self.msd)+1)
        y = self.msd
        #D0 = self.msd[5]/0.1 # FIXME: generalize
        #y = self.msd/(x*D0) # FIXME: generalize lower limit as brownian step from h5 file
         
        #y = 5*self.msd/x # FIXME: generalize lower limit as brownian step from h5 file
         
        #y = []
        #for xx,yy in zip(x, self.msd):
        #    y.append( yy/(xx * D0) )
 
        plt.xlabel("t")
        plt.ylabel("MSD")
         
        #plt.grid(color='#cccccc', linestyle='--', linewidth=1)
         
        plt.yscale('log')
        plt.xscale('log')
         
        plt.xlim([0.1, 10000])
        #plt.xlim([0,50])
        #plt.plot(x, y, 'ro')
        plt.plot(x, y)
        plt.show()
        plt.savefig(outfilename)
 
    def _uglyplot_msr(self, outfilename):
        import matplotlib.pyplot as plt
        plt.figure(figsize=(4,3))
        # FIXME: replace 10 with actual time, 0.1/index 5 with Brownian step/summary step
        x = np.linspace(start=0, stop=10000, num=len(self.msr)+1)[1:]
        y = self.msr
        plt.xlabel("t")
        plt.ylabel("MSR")
        plt.yscale('log')
        plt.xscale('log')
        plt.xlim([0.1, 10000])
        plt.plot(x, y)
        #plt.show()
        plt.savefig(outfilename)
        
    def _uglyplot_skstatic(self, outfilename):
        import matplotlib.pyplot as plt
        plt.figure(figsize=(5,4))
        # FIXME: replace 10 with actual time, 0.1/index 5 with Brownian step/summary step
        
#         x = self.sk_static[0]
#         y = self.sk_static[1]
        
        x = np.linspace(start=0, stop=40, num=len(self.sk_static)+1)[1:]
        y = self.sk_static
        plt.xlabel("q")
        plt.ylabel("sk")
#         plt.yscale('log')
#         plt.xscale('log')
        plt.xlim([0, 40])
        plt.plot(x, y)
        #plt.show()
        plt.savefig(outfilename)
 
 
#     def _uglyplot_energy(self, outfilename):
#         import matplotlib.pyplot as plt
#         plt.figure(figsize=(4,3))
#          
#         # FIXME: replace 10 with actual time, 0.1/index 5 with Brownian step/summary step
#         y = self.energy
#         x = np.linspace(start=0, stop=1, num=len(y)+1)[1:]
#          
#         plt.xlabel("t")
#         plt.ylabel("E")
#      
#         #plt.grid(color='#cccccc', linestyle='--', linewidth=1)        
#         #plt.yscale('log')
#         #plt.xscale('log')
#                  
#         #plt.plot(x, y, 'ro')
#         plt.plot(x, y)
#         #plt.show()
#         plt.savefig(outfilename)
# 
# 
#     def _plot_packing_fraction(self):
#         import matplotlib.pyplot as plt
#         plt.figure(figsize=(4,3))
#         
#         # FIXME: replace 10 with actual time, 0.1/index 5 with Brownian step/summary step
#         x = np.linspace(start=0, stop=50, num=len(self.pf)+1)[1:]
#         y = self.pf
# 
#         plt.xlabel("t")
#         plt.ylabel("phi")
#         
#         #plt.grid(color='#cccccc', linestyle='--', linewidth=1)
#         
#         #plt.xlim([0.1,50])
#         #plt.xlim([0,50])
#         #plt.plot(x, y, 'ro')
#         plt.plot(x, y)
#         plt.show()
#         #plt.savefig(outfilename)
# 
# 
#     #def msd(self, force_recalc=False): # FIXME: implement
# 
# 
#     def ftplot(self, ftdata, outfilename):
#         z = ftdata
#         import matplotlib.pyplot as plt
#         #plt.figure(figsize=(4,3),dpi=600)
#         fig, ax = plt.subplots(figsize=(4,3),dpi=300)
#         #plt.gca().tight_layout()
#         plt.gcf().subplots_adjust(bottom=0.15)
#         
#         plt.xlabel(r'$q_1$')
#         plt.ylabel(r'$q_2$')
#         
#         """
#         #plt.yscale('log')
#         #plt.xscale('log')
#         x = np.transpose(z)[0]
#         y = np.transpose(z)[1]
#         n = len(x)
#         #plt.plot(x[1:],y[1:], 'ro')
#         plt.hist2d(x[2:n/2], y[2:n/2], bins=200, norm=LogNorm())
#         
#         #plt.hist(z.ravel(), bins=100)
#         #plt.plot(z.ravel(), 'ro')
#         """
#         
#         from matplotlib.colors import LogNorm
# 
#         #(xmin, xmax, ymin, ymax) = (-5,5,-5,5)
#         #(xres, yres) = (100,100)
#         #(xmin, xmax, ymin, ymax) = (-2.5,2.5,-2.5,2.5)
#         # MW system
#         '''
#         #s=4.22/2
#         #s=2.965 # MW phi=0.7 (kmax sc)
#         #s=2.74587 # MW phi=0.6 (kmax sc)
#         s=4.21 # kmax (2*b)
#         
#         # kmin
#         s=4 * 0.0858086 # 4*kmin for phi=0.6 in sc crystal
#         plt.xticks([-0.2,0,0.2])
#         '''
# 
#         # P3 system kmax
#         phi = 0.697
#         kmaxagreen = 3.73131/np.sqrt(phi)
#         kmaxbgreen = 5.13571/np.sqrt(phi)
#         kmaxablue  = 3.17404/np.sqrt(phi)
#         kmaxbblue  = 9.7687/np.sqrt(phi) # biggest
#         s = kmaxbblue # no crystal 
#         
#         # P3 system
#         #s=11.7 # P3 phi=0.74 W!
#         #s=11.55 # P3 phi=0.75 W!
#         #s=11.4 # P3 phi=0.76 W!
#         #s=11.25 # P3 phi=0.77 W!
#         #s=11.21 # 0.76
#         #s=11.13 # 0.77
# 
# 
#         #s=7*0.146 # kmin
#         
#         (xmin, xmax, ymin, ymax) = (-s,s,-s,s)
#         
#         
#         #(xres, yres) = (200,200)
#         (xres, yres) = (500,500)
#         extent = (xmin, xmax, ymin, ymax)
#         x = np.arange(xmin, xmax, (xmax-xmin)/xres)
#         y = np.arange(ymin, ymax, (ymax-ymin)/yres)
#         z = np.array([np.absolute(z(xx,yy))**2 for yy in y for xx in x])
#         N = int(len(z)**.5)
#         z = z.reshape(N, N)
# 
#         plt.imshow(z, extent = extent, norm = LogNorm(), origin='lower')
#         plt.colorbar()
#         plt.savefig(outfilename)
# 
#         # DRAW CIRCLES P3-2
#         circle1 = plt.Circle((0,0),kmaxagreen,color='green',fill=False, alpha=0.6, lw=2) # color='white'; alpha=0.6
#         circle2 = plt.Circle((0,0),kmaxbgreen,color='green',fill=False, alpha=0.6, lw=2)
#         circle3 = plt.Circle((0,0),kmaxablue,color='blue',fill=False, alpha=0.6, lw=2)
#         circle4 = plt.Circle((0,0),kmaxbblue,color='blue',fill=False, alpha=0.6, lw=2)
#         #circle1 = plt.Circle((0,0),kmaxa,color='red',fill=False)
#         ax.add_patch(circle1)
#         ax.add_patch(circle2)
#         ax.add_patch(circle3)
#         ax.add_patch(circle4)
#         # END OF CIRCLES P3-2
#         
#         #plt.show()
#         plt.savefig(outfilename + ".circle.pdf")
#         
#         
#         #plt.show()
#         #plt.savefig(outfilename)
#         plt.close()
#         return z
        

    
    def __init__(self, input_filename):
        print("inside Statprop")
        self.bf = BedfileReader(input_filename)
        
        # FIXME: save to h5, not to csv with numpy
        
        # FIXME: temporary deact

        self.msd = np.array(list(self.bf.msd()))
        np.savetxt(input_filename + '-msd.dat', self.msd, delimiter=',') # save msd to file (quick & dirty method.. just for testing)
        self._uglyplot_msd(input_filename + "-msd.eps")
   
        self.msr = np.array(list(self.bf.msr()))
        np.savetxt(input_filename + '-msr.dat', self.msr, delimiter=',') # save msr to file (quick & dirty method.. just for testing)
        self._uglyplot_msr(input_filename + "-msr.eps")
        
        self.sk_static = np.array(list(self.bf.sk_static()))
        np.savetxt(input_filename + '-sk_static.dat', self.sk_static, delimiter=',') # save msr to file (quick & dirty method.. just for testing)
        self._uglyplot_skstatic(input_filename + "-sk_static.eps")
        
        
        """
        self.energy = np.array(list(self.bf.energy_conservation()))
        np.savetxt(input_filename + '-energy.dat', self.energy, delimiter=',')
        self._uglyplot_energy(input_filename + "-energy.eps")
        """
        
        #self.pf = self.bf.calc_packing_fraction()
        #np.savetxt(input_filename + '-packing_fraction2.dat', self.pf, delimiter=',') # save msd to file (quick & dirty method.. just for testing)
        #self._plot_packing_fraction()
        
        #self.ftplot(self.bf.ft_static(), input_filename + '-ft.pdf')
        #zlast = self.ftplot(self.bf.ft_last(), input_filename + '-ft-last.pdf')
        #np.savetxt(input_filename + '-ft-last.dat', zlast, delimiter=',')
        
        """ #calc ft for each saved time 
        fts = self.bf.ft()
        i=0
        for ft in fts:
            i+=1
            zdata = self.ftplot(ft, input_filename + ('-ft-%d.svg' % i))
            np.savetxt(input_filename + ('-ft-%d.dat' % i), zdata, delimiter=',')
        """
