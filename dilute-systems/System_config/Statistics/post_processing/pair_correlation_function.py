import numpy as np
import matplotlib.pyplot as plt
#TODO skalieren und achsenbeschriftung :)
data1_x = np.loadtxt('/home/newton/sophia/Desktop/MASTER/spheres/spheres_ugly/data/paircorrelation_radialdistributionfunction/n256phi0.494tb0.01_frame30von51.txt', usecols= 1)
data1_y = np.loadtxt('/home/newton/sophia/Desktop/MASTER/spheres/spheres_ugly/data/paircorrelation_radialdistributionfunction/n256phi0.494tb0.01_frame30von51.txt', usecols= 2)

data2_x = np.loadtxt('/home/newton/sophia/Desktop/MASTER/spheres/spheres_ugly/data/paircorrelation_radialdistributionfunction/n500_phi0.3_tb0.01.txt', usecols= 1)
data2_y = np.loadtxt('/home/newton/sophia/Desktop/MASTER/spheres/spheres_ugly/data/paircorrelation_radialdistributionfunction/n500_phi0.3_tb0.01.txt', usecols= 2)

data3_x = np.loadtxt('/home/newton/sophia/Desktop/MASTER/spheres/spheres_ugly/data/paircorrelation_radialdistributionfunction/n500phi0.45tb0.01.txt', usecols= 1)
data3_y = np.loadtxt('/home/newton/sophia/Desktop/MASTER/spheres/spheres_ugly/data/paircorrelation_radialdistributionfunction/n500phi0.45tb0.01.txt', usecols= 2)

data4_x = np.loadtxt('/home/newton/sophia/Desktop/MASTER/spheres/spheres_ugly/data/paircorrelation_radialdistributionfunction/n500phi0.494tb0.01.txt', usecols= 1)
data4_y = np.loadtxt('/home/newton/sophia/Desktop/MASTER/spheres/spheres_ugly/data/paircorrelation_radialdistributionfunction/n500phi0.494tb0.01.txt', usecols= 2)

plt.plot(data1_x/2., data1_y, label = 'n256phi0.494tb0.01') #skalieren da r = 1 und sigma 2r ist deswegen /2.
plt.plot(data2_x/2., data2_y, label = 'n500_phi0.3_tb0.01')
plt.plot(data3_x/2., data3_y, label = 'n500phi0.45tb0.01')
plt.plot(data4_x/2., data4_y, label = 'n500phi0.494tb0.01')
plt.title('Pair correlation function for hard spheres')
plt.xlabel(r'$ r / \sigma $', fontsize = 18)
plt.ylabel('radial distribution function $g(r)$', fontsize = 18)
plt.legend()
plt.show()