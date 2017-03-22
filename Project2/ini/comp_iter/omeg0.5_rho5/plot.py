#!/usr/bin/env python2

import numpy as np
import matplotlib.pyplot as plt

thedat = np.loadtxt("times.txt",skiprows=0)
grid = thedat[:,0]
time_cyc = thedat[:,1]
time_cla = thedat[:,2]
iter_cyc = thedat[:,3]
iter_cla = thedat[:,4]

fig, ax1 = plt.subplots()
ax1.semilogy(grid, iter_cyc, label="Cyclical",linewidth=2.0,color='orangered')
ax1.semilogy(grid, iter_cla, label="Classical",linewidth=2.0,color='midnightblue')
ax1.set_xlabel('Grid Size')
ax1.set_ylabel('log(# iterations)')
ax1.set_title('Iterations until Convergence')
ax1.set_xlim(50,500)
ax1.legend()

plt.savefig('iterations.eps')
plt.show()
