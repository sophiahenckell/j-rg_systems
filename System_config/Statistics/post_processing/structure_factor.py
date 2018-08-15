import numpy as np
import matplotlib.pyplot as plt

x1_static = np.loadtxt('/data/scc/sophia/Second_MSDCalculations2-n500-t4000-phi0.494-k1.2-tb0.01-ts10/Ellipse-Benchmark-n500-k1.2-phi0.494-id1-equilibrate.h5.dump-static-sk.dat', usecols= 0)
y1_static = np.loadtxt('/data/scc/sophia/Second_MSDCalculations2-n500-t4000-phi0.494-k1.2-tb0.01-ts10/Ellipse-Benchmark-n500-k1.2-phi0.494-id1-equilibrate.h5.dump-static-sk.dat', usecols= 1)
x2_static = np.loadtxt('/data/scc/sophia/Second_MSDCalculations2-n500-t4000-phi0.494-k1.2-tb0.01-ts10/Ellipse-Benchmark-n500-k1.2-phi0.494-id2-equilibrate.h5.dump-static-sk.dat', usecols= 0)
y2_static = np.loadtxt('/data/scc/sophia/Second_MSDCalculations2-n500-t4000-phi0.494-k1.2-tb0.01-ts10/Ellipse-Benchmark-n500-k1.2-phi0.494-id2-equilibrate.h5.dump-static-sk.dat', usecols= 1)
x3_static = np.loadtxt('/data/scc/sophia/Second_MSDCalculations2-n500-t4000-phi0.494-k1.2-tb0.01-ts10/Ellipse-Benchmark-n500-k1.2-phi0.494-id3-equilibrate.h5.dump-static-sk.dat', usecols= 0)
y3_static = np.loadtxt('/data/scc/sophia/Second_MSDCalculations2-n500-t4000-phi0.494-k1.2-tb0.01-ts10/Ellipse-Benchmark-n500-k1.2-phi0.494-id3-equilibrate.h5.dump-static-sk.dat', usecols= 1)
x4_static = np.loadtxt('/data/scc/sophia/Second_MSDCalculations2-n500-t4000-phi0.494-k1.2-tb0.01-ts10/Ellipse-Benchmark-n500-k1.2-phi0.494-id4-equilibrate.h5.dump-static-sk.dat', usecols= 0)
y4_static = np.loadtxt('/data/scc/sophia/Second_MSDCalculations2-n500-t4000-phi0.494-k1.2-tb0.01-ts10/Ellipse-Benchmark-n500-k1.2-phi0.494-id4-equilibrate.h5.dump-static-sk.dat', usecols= 1)
x5_static = np.loadtxt('/data/scc/sophia/Second_MSDCalculations2-n500-t4000-phi0.494-k1.2-tb0.01-ts10/Ellipse-Benchmark-n500-k1.2-phi0.494-id5-equilibrate.h5.dump-static-sk.dat', usecols= 0)
y5_static = np.loadtxt('/data/scc/sophia/Second_MSDCalculations2-n500-t4000-phi0.494-k1.2-tb0.01-ts10/Ellipse-Benchmark-n500-k1.2-phi0.494-id5-equilibrate.h5.dump-static-sk.dat', usecols= 1)

data = []
for i in range(len(x1_static)):
    data = np.append(data, np.mean((y1_static[i], y2_static[i], y3_static[i], y4_static[i], y5_static[i])))

plt.plot(x1_static, data, label = 'averaged over 5 Tasks')
plt.plot(x3_static,y3_static, label = 'example curve (id: 3)')
plt.title('Static structure factor packing fraction 0.494')
plt.xlabel(r'$ \mid q \mid $', fontsize = 18)
plt.ylabel('static structure factor')
plt.legend()
plt.show()





