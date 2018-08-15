#TODO werte checken und kleines Zeitintervall ergänzen
#TODO skalieren und achsenbeschriftung

import numpy as np
import matplotlib.pyplot as plt
# def create_image(directory):
#     data1 = np.loadtxt(directory+'msd.dat', dtype=float, usecols= 0)
#     data2 = np.loadtxt(directory + 'msd.dat', dtype=float, usecols=1)
#     print(data1, data2)
#     plt.plot(data1, data2)
#     plt.show()
#     return None
#
# create_image('/data/scc/sophia/Second_MSDCalculations2-n256-t4000-phi0.494-k1.2-tb0.01-ts10/')

'''dieses Programm schreibt ganz grob die MSDs für Herrn Prof. Dr. Fuchs zusammen'''
data1_x = np.loadtxt('/home/newton/sophia/Desktop/MASTER/spheres/spheres_ugly/data/MSD_Data_6June/msd_n256t4000phi0494tb0.01ts10_msd.dat', dtype=float, usecols= 0)
data1_y = np.loadtxt('/home/newton/sophia/Desktop/MASTER/spheres/spheres_ugly/data/MSD_Data_6June/msd_n256t4000phi0494tb0.01ts10_msd.dat', dtype=float, usecols= 1)

data2_x = np.loadtxt('/home/newton/sophia/Desktop/MASTER/spheres/spheres_ugly/data/MSD_Data_6June/msd_n500t4000phi03tb0.01ts10_msd.dat', dtype=float, usecols= 0)
data2_y = np.loadtxt('/home/newton/sophia/Desktop/MASTER/spheres/spheres_ugly/data/MSD_Data_6June/msd_n500t4000phi03tb0.01ts10_msd.dat', dtype=float, usecols= 1)

data3_x = np.loadtxt('/data/scc/sophia/Second_MSDCalculations2-n500-t4000-phi0.45-k1.2-tb0.01-ts10/msd.dat', dtype=float, usecols= 0)
data3_y = np.loadtxt('/data/scc/sophia/Second_MSDCalculations2-n500-t4000-phi0.45-k1.2-tb0.01-ts10/msd.dat', dtype=float, usecols= 1)

data4_x = np.loadtxt('/data/scc/sophia/Second_MSDCalculations2-n500-t4000-phi0.494-k1.2-tb0.01-ts10/msd.dat', dtype=float, usecols= 0)
data4_y = np.loadtxt('/data/scc/sophia/Second_MSDCalculations2-n500-t4000-phi0.494-k1.2-tb0.01-ts10/msd.dat', dtype=float, usecols= 1)

##ERWEITERUNG UM KLEINE ZEITSCHRITTE
data5_x = np.loadtxt('/data/scc/sophia/Second_MSDCalculations2-n500-t4000-phi0.494-k1.2-tb0.01-ts5/msd.dat', dtype = float, usecols= 0)
data5_y = np.loadtxt('/data/scc/sophia/Second_MSDCalculations2-n500-t4000-phi0.494-k1.2-tb0.01-ts5/msd.dat', dtype = float, usecols= 1)

data6_x = np.loadtxt('/data/scc/sophia/Second_MSDCalculations2-n500-t10-phi0.494-k1.2-tb0.01-ts0.1/msd.dat', dtype = float, usecols= 0)
data6_y = np.loadtxt('/data/scc/sophia/Second_MSDCalculations2-n500-t10-phi0.494-k1.2-tb0.01-ts0.1/msd.dat', dtype = float, usecols= 1)

plt.loglog(data1_x, data1_y, label= 'n256_phi0.494_tb0.01')
plt.loglog(data2_x, data2_y, label = 'n500_phi0.3_tb0.01')
plt.loglog(data3_x, data3_y, label = 'n500_phi0.45_tb0.01')
plt.loglog(data4_x, data4_y, label = 'n500_phi0.494_tb0.01')
plt.loglog(data5_x, data5_y, label = 'n500-phi0.494_tb0.01 saving 5')
plt.loglog(data6_x, data6_y, label = 'n500-phi0.494_tb0.01 saving 0.1')
plt.title('MSD')
plt.xlabel(r'$ \tau $', fontsize = 18)
plt.ylabel(r'$ \langle \delta r ^2 (\tau)\rangle $ ', fontsize = 18)
plt.legend(loc ='upper left', fontsize = 12)
plt.show()
