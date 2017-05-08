import numpy as np
import MDlib_thermal
Cut = 3.0
np.random.seed = 795344

def InitPositions(N, L):
    Pos = np.zeros((N,3), float)
    NLat = int(N**(1./3.) + 1.)
    r = L * (np.arange(NLat, dtype=float)/NLat - 0.5)
    i = 0
    for x in r:
        for y in r:
            for z in r:
                Pos[i] = np.array([x,y,z], float)
                i += 1
                if i >= N:
                    return Pos
    return Pos


def RescaleVel(Vel, T):
    """
    for use in NVT solver
    """
    Vel = Vel - Vel.mean(axis=0)
    KE = 0.5 * np.sum(Vel * Vel)
    VScale = np.sqrt(1.5 * len(Vel) * T / KE)
    Vel = Vel * VScale
    return Vel  


def InitVel(N, T):
    """
	returns randome velocities
	"""
    Vel = np.random.rand(N, 3)
    Vel = RescaleVelocities(Vel, T)
    return Vel


def InitAccel(Pos, L):
    """
	initialize acceleration
	"""
    Accel = np.zeros_like(Pos)
    PEnergy, Accel = MDlib.Force(Pos, L, Cut, Accel)
    return Accel  

def RunTest():
    N = 108
    rho = 0.05
    L = (N / rho)**(1./3.)
    Temp = 1.0
    dt = 0.001
	mavg = 10000
	nequ = 100000
    RescaleSteps = 1000
    WriteSteps = 100
    MaxSteps = WriteSteps*10000
 
    Pos = InitPositions(N, L)
    Vel = InitVel(N, Temp)
    Accel = InitAccel(Pos, L)
	Jxi = []
	Jyi = []
	Jzi = []

    i = 0
    while i < MaxSteps or MaxSteps <= 0:
	# FORTRAN LIB CALL
        Pos, Vel, Accel, Jx, Jy, Jz, KEnergy, PEnergy = MDlib.IntegrateVerletThermal(Pos, Vel, Accel, L, Cut, dt)
		Jxi.append(Jx)
		Jyi.append(Jy)
		Jzi.append(Jz)

        i += 1
        if i % RescaleSteps == 0:
            Vel = RescaleVel(Vel, Temp)
	
	JiJj=0.0
	heatcorr = np.zeros(mavg)
	# just getting thermal cond where its statistically dominant < ma
	for i in range(0, mavg):
		heatcorr(i)=0.0D00
		for j in range(nequ, MaxSteps):
			heatcorr(i)=heatcorr(i)+Jxi(i+j)*Jxi(j)+Jyi(i+j)*Jyi(j)+Jzi(i+j)*Jzi(j)
		heatcorr(i)=heatcorr(i)/(MaxSteps-i-nequ)
		JiJj=JiJj+heatcorr(i)
	lambdacond=dt/Temp**2*JiJj
	
if __name__ == '__main__':
    RunTest()

    
