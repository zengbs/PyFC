import numpy as np
from pri2con import Pri2Con
from sphere import SphericalSphere

def _DumpFile(out):


    # Output file name
    filename = "UM_IC"

    # Number of cells along x/y/z
    Nx = out.shape[0] 
    Ny = out.shape[1] 
    Nz = out.shape[2] 

    # Box size
    Lx = 100
    Ly = 100
    Lz = 100

    # Radius of sphere
    Radius = 40

    # Peak gas density
    Rho0_g = 1

    # Square of sound speed inside the isotheraml sphere
    CsSqr = 0.2
  
    #
    Rho0_DM = 1

    # Scale radius of NFW profile
    ScaleRadius = 1

    # Number of cells along each dimension of the input grid.
    N = np.array([Nx, Ny, Nz], dtype='int')
    L = np.array([Lx, Ly, Lz], dtype='float64') 

    FluidInBox = SphericalSphere( L, N, Radius, Rho0_g, CsSqr, Rho0_DM, ScaleRadius )

    FluidInBox.tofile(filename)

out=np.zeros((1024,1024,1024))
_DumpFile(out)
