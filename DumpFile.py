import numpy as np
from pri2con import *

def _DumpFile(out):

    # Number of cells along x/y/z
    Nx = out.shape[0] 
    Ny = out.shape[1] 
    Nz = out.shape[2] 

    # Number of cells along each dimension of the input grid.
    N = np.array([Nx, Ny, Nz], dtype='int')
    
    # Left edge and right edge coordinates of the desired
    # simulation domain which will be used in GAMER.
    le = np.array([ 0,  0,  0])
    re = np.array([Lx, Ly, Lz])
    
    # Cell spacing
    delta = (re-le)/N  # [1/128 1/128 1/128]
    
    # Construct the grid cell edge coordinates
    x = np.linspace(le[0], re[0], num=N[0]+1) 
    y = np.linspace(le[1], re[1], num=N[1]+1)
    z = np.linspace(le[2], re[2], num=N[2]+1)
    
    # Find the grid cell midpoints
    x = 0.5*(x[1:]+x[:-1])
    y = 0.5*(y[1:]+y[:-1])
    z = 0.5*(z[1:]+z[:-1])
    
    # Set the center of box as origin
    x -= 0.5*Lx
    y -= 0.5*Ly
    z -= 0.5*Lz

    # Use the 1-D coordinate arrays to consruct 3D coordinate arrays
    #xx, yy, zz = np.meshgrid(x, y, z, sparse=False, indexing='xy')
    Fluid = np.zeros((5, x.size, y.size, z.size), dtype=np.float64)

    Rho  = out
    Ux   = np.full_like( Rho, 0.0 )
    Uy   = np.full_like( Rho, 0.0 )
    Uz   = np.full_like( Rho, 0.0 )
    Temp = np.full_like( Rho, 1   )
    Pres = Temp*Rho

    Dens, MomX, MomY, MomZ, Engy = Pri2Con( Rho, Ux, Uy, Uz, Pres )
  
    filename = "UM_IC"

    FluidVar = np.zeros((5, Nx, Ny, Nz), dtype=np.float32)

    FluidVar[0] = Rho
    FluidVar[1] = Ux
    FluidVar[2] = Uy
    FluidVar[3] = Uz
    FluidVar[4] = Temp

    FluidVar.tofile(filename)
