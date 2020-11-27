import numpy as np
from density_profile import *
from pri2con import Pri2Con


def SphericalSphere( L, N, ParaPhy, ParaNum ):

    Nx, Ny, Nz = N
    Lx, Ly, Lz = L

    # Center of sphere
    Center = np.array([Lx,Ly,Lz])*0.5
 
    # Unbundle numerical parameters
    BPoint, CoarseDr = ParaNum   

    # Coarse grid
    Coarse_r = np.arange(BPoint, 0.5*Lx, CoarseDr)

    InRho   = NumericalDensity( Coarse_r, ParaPhy )
    InUX    = 0 
    InUy    = 0
    InUz    = 0
    InPres  = CsSqr*InRho

    OutRho  = 1e-3*np.interp( Radius, Coarse_r, InRho )
    OutUX   = 0 
    OutUy   = 0
    OutUz   = 0
    OutPres = 1e3*CsSqr*OutRho

    InDens,   InMomX,  InMomY,  InMomZ,  InEngy = Pri2Con(  InRho,  InUX,  InUy,  InUz, InPres  )
    OutDens, OutMomX, OutMomY, OutMomZ, OutEngy = Pri2Con( OutRho, OutUX, OutUy, OutUz, OutPres )


    # Left edge and right edge coordinates of the desired
    # simulation domain which will be used in GAMER.
    le = np.array([ 0,  0,  0])
    re = np.array([Lx, Ly, Lz])
    
    # Cell spacing
    delta = (re-le)/N  # [1/128 1/128 1/128]
    
    FluidInBox = np.zeros((6, Nx, Ny, Nz), dtype=np.float32)

    for i in range(Nx):
        for j in range(Ny):
            for k in range(Nz):

                x = (i+0.5)*delta[0] 
                y = (j+0.5)*delta[1]
                z = (k+0.5)*delta[2]


                r = np.sqrt((x-Center[0])**2 + (y-Center[1])**2 + (z-Center[2])**2)

                if ( r < Radius ):
                     FluidInBox[0][i][j][k] = np.interp( r, Coarse_r, InDens )
                     FluidInBox[1][i][j][k] = np.interp( r, Coarse_r, InMomX ) 
                     FluidInBox[2][i][j][k] = np.interp( r, Coarse_r, InMomY ) 
                     FluidInBox[3][i][j][k] = np.interp( r, Coarse_r, InMomZ ) 
                     FluidInBox[4][i][j][k] = CsSqr*FluidInBox[0][i][j][k]
                     FluidInBox[5][i][j][k] = CsSqr
                else:
                     FluidInBox[0][i][j][k] = OutDens
                     FluidInBox[1][i][j][k] = OutMomX
                     FluidInBox[2][i][j][k] = OutMomY
                     FluidInBox[3][i][j][k] = OutMomZ
                     FluidInBox[4][i][j][k] = OutEngy
                     FluidInBox[5][i][j][k] = CsSqr
                 
                
    return FluidInBox
