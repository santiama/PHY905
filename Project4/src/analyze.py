import numpy as np
import matplotlib.pyplot as plt

array = np.loadtxt('array.dat') 
T = array[:,0]
K = array[:,1]
U = array[:,2]
P = array[:,3]
disp = array[:,4]

data = np.loadtxt('data.dat')
rho = data[0]
T_IC = data[1]
eq_time = data[2]
time_end = data[3]
dt = data[4]
r_cut = data[5]
r_verlet = data[6]
N_bin = data[7]
N = data[8]

entries = np.size(T)
data_block_size = 50
def data_block(x):
 correct_length = entries-entries%data_block_size 
 xx = x[:correct_length].reshape(-1,data_block_size)
 avgs = np.mean(xx, axis = 0)
 var = np.std(xx,axis = 0)
 avg = np.mean(avgs, axis = 0)
 error = np.std(avgs)/np.sqrt(np.size(avgs))
 return avgs, var, avg, error

T_avgs, T_var, T_avg, T_error = data_block(T)
K_avgs, K_var, K_avg, K_error = data_block(K)
U_avgs, U_var, U_avg, U_error = data_block(U)
n = np.size(disp)
x = np.arange(n)
A = np.vstack([x, np.ones(n)]).T
model, resid = np.linalg.lstsq(A, disp)[:2]
r2 = 1 - resid / (np.size(disp) * np.var(disp))
D = model[0]/(6*dt) #diffusion coefficient
D_error = model[0]*(1-r2[0])*100/(6*dt)

time = np.linspace(0,entries,entries)

def plot(x, y):
 plt.figure()
 plt.ylabel(r'y', fontdict=font)
 plt.xlabel(r'x', fontdict=font)
 plt.plot(x, y, 'b', label=r"line")
 plt.legend()
 plt.show()

def plot_energy():
 plt.figure()
 plt.ylabel(r'Energy $(\epsilon)$', fontdict=font)
 plt.xlabel(r'Time ($\Delta t$)', fontdict=font)
 plt.plot(U, 'r', label='Potential energy')
 plt.plot(K, 'g', label='Kinetic energy')
 plt.plot(K+U, 'b', label='Total energy')
 plt.legend(loc=7)
 plt.show()


def plot_pcf():
 histogram = np.loadtxt('histogram.dat')
 bin_size = r_verlet / N_bin
 pcf = np.ones(int(N_bin)-1)
 pcf_x = np.ones(int(N_bin)-1)
 for i in range(1, len(pcf)+1):
  pcf_x[i-1] = bin_size*i
  pcf[i-1] = 2.0*histogram[i-1] / (4*np.pi*bin_size*N*rho*(i*bin_size)**2)
 plt.figure()
 plt.xlim(0,r_verlet)
 plt.ylabel(r'$g(r)$', fontdict=font)
 plt.xlabel(r'$r/\sigma$', fontdict=font)
 plt.hlines(1, 0, r_verlet, linestyles='dashed')
 plt.plot(pcf_x ,pcf, 'r', label=r'Pair correlation function $g(r)$')
 plt.legend()
 plt.show()

def estimated_autocorrelation(x):
 n = len(x)
 variance = x.var()
 x = x-x.mean()
 r = np.correlate(x, x, mode = 'full')[-n:]
 assert np.allclose(r, np.array([(x[:n-k]*x[-(n-k):]).sum() for k in range(n)]))
 result = r/(variance*(np.arange(n, 0, -1)))
 return result

def find_nearest(array,value):
 idx = (np.abs(array-value)).argmin()
 return idx, array[idx]

def plot_temperature():
 plt.xlabel(r'Time ($\Delta t$)', fontdict=font)
 plt.ylabel(r'Temperature ($\epsilon / k$)', fontdict=font)
 plt.plot(T, 'r', label=r'Temperature')
 plt.legend()
 plt.show()

print 'T ', T_avg, T_error
print 'D ', D, D_error

