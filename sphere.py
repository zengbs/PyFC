import numpy as np
from coordinate import *
from DensityProfile import *


def SphericalSymmetricSphere( Fluid, L, N, Radius, InsideSphere, OutsideSphere ):

    Nx = N[0]
    Ny = N[1]
    Nz = N[2]

    Lx = L[0]
    Ly = L[1]
    Lz = L[2]

    InRho   = InsideSphere[0]
    InUX    = InsideSphere[1]
    InUy    = InsideSphere[2]
    InUz    = InsideSphere[3]
    InPres  = InsideSphere[4]

    OutRho  = OutsideSphere[0]
    OutUX   = OutsideSphere[1]
    OutUy   = OutsideSphere[2]
    OutUz   = OutsideSphere[3]
    OutPres = OutsideSphere[4]


    # Number of cells along each dimension of the input grid.
    ddims = np.array([Nx, Ny, Nz], dtype='int')
    
    # Left edge and right edge coordinates of the desired
    # simulation domain which will be used in GAMER.
    le = np.array([ 0,  0,  0])
    re = np.array([Lx, Ly, Lz])
    
    # Cell spacing
    delta = (re-le)/ddims  # [1/128 1/128 1/128]
    
    # Construct the grid cell edge coordinates
    x = np.linspace(le[0], re[0], num=ddims[0]+1) 
    y = np.linspace(le[1], re[1], num=ddims[1]+1)
    z = np.linspace(le[2], re[2], num=ddims[2]+1)
    
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


    for i in range(Nx):
        for j in range(Ny):
            for k in range(Nz):

                x = i*ddims[0]
                y = j*ddims[1]
                z = k*ddims[2]

                r = np.sqrt(x**2 + y**2 + z**2)

                if ( r < Radius ):
                     Fluid[0][i][j][k] = InsideSphere[0]
                     Fluid[1][i][j][k] = InsideSphere[1]
                     Fluid[2][i][j][k] = InsideSphere[2]
                     Fluid[3][i][j][k] = InsideSphere[3]
                     Fluid[4][i][j][k] = InsideSphere[4]
                elif ( r >= Radius ):
                     Fluid[0][i][j][k] = OutsideSphere[0]
                     Fluid[1][i][j][k] = OutsideSphere[1]
                     Fluid[2][i][j][k] = OutsideSphere[2]
                     Fluid[3][i][j][k] = OutsideSphere[3]
                     Fluid[4][i][j][k] = OutsideSphere[4]
                else:
                     print("NaN encountered !!")
                
    return Fluid
