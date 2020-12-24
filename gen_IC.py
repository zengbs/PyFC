import pyFC
import matplotlib.pyplot as plt
import matplotlib.cm as cm

Nx = 128 
Ny = 128 
Nz = 128 

fc = pyFC.LogNormalFractalCube(ni=Nx, nj=Ny, nk=Nz, kmin=10, mean=1)

fc.gen_cube()

pyFC.write_cube(fc=fc, fname='UM_IC', app=True, prec='single')
