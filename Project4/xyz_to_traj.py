from ase.io import read, write

tmp = read('pos.xyz',index=':')
write('pos.traj',tmp)
