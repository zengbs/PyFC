import numpy as np
from density_profile import *
from pri2con import Pri2Con


def SphericalSphere( L, N, ParaPhy, ParaNum, Precision ):

    Nx, Ny, Nz = N
    Lx, Ly, Lz = L

    # Center of sphere
    Center = np.array([Lx,Ly,Lz])*0.5
 
    # Unbundle numerical parameters
    BPoint, CoarseDr, PlotRadius = ParaNum   

    # Unbundle physical parameters
    Radius_g, Rho0_g, Sigma_g, Lambda, Kappa, Temp_g, Constant, Phi0, DevPhi0 = ParaPhy

    # Coarse grid
    Coarse_r = np.arange(BPoint, 0.5*Lx, CoarseDr)


    # Left edge and right edge coordinates of the desired
    # simulation domain which will be used in GAMER.
    le = np.array([ 0,  0,  0])
    re = np.array([Lx, Ly, Lz])
    
    # Cell spacing
    delta = (re-le)/N  # [1/128 1/128 1/128]

    ####################################
    ############   Fluid  ##############
    ####################################

    # Fluid inside box
    InRho, Nothing   = NumericalDensity( Coarse_r, ParaPhy )
    InUX    = 0 
    InUy    = 0
    InUz    = 0
    InPres  = Sigma_g*InRho

    # Fluid outside box
    OutRho  = 1e-3*np.interp( PlotRadius, Coarse_r, InRho )
    OutUX   = 0 
    OutUy   = 0
    OutUz   = 0
    OutPres = 1e3*Sigma_g*OutRho

    # Conversion
    InDens,   InMomX,  InMomY,  InMomZ,  InEngy = Pri2Con(  InRho,  InUX,  InUy,  InUz, InPres  )
    OutDens, OutMomX, OutMomY, OutMomZ, OutEngy = Pri2Con( OutRho, OutUX, OutUy, OutUz, OutPres )

    
    FluidInBox = np.zeros((5, Nx, Ny, Nz), dtype=Precision)


    ####################################
    ##########   Potential  ############
    ####################################
    PotInBox   = np.zeros((Nx, Ny, Nz), dtype=Precision)

    Rho0_D, Sigma_D, Radius_D  = Free2DerivedPara( Rho0_g, Sigma_g, Radius_g, Lambda, Kappa )
    Psi0       = Phi0      / Sigma_D / Sigma_D
    DevPsi0    = DevPhi0   / Sigma_D / Sigma_D

    Potential = NumericalTotalPotential( Coarse_r, Kappa, Lambda, Constant, Psi0, DevPsi0 )
    Potential *= Sigma_D*Sigma_D



    for i in range(Nx):
        print("i=%d" % i)
        for j in range(Ny):
            for k in range(Nz):

                x = (i+0.5)*delta[0] 
                y = (j+0.5)*delta[1]
                z = (k+0.5)*delta[2]
                r = np.sqrt((x-Center[0])**2 + (y-Center[1])**2 + (z-Center[2])**2)

                if ( r < PlotRadius ):
                     FluidInBox[0][i][j][k] = np.interp( r, Coarse_r, InDens )
                     FluidInBox[1][i][j][k] = np.interp( r, Coarse_r, InMomX ) 
                     FluidInBox[2][i][j][k] = np.interp( r, Coarse_r, InMomY ) 
                     FluidInBox[3][i][j][k] = np.interp( r, Coarse_r, InMomZ ) 
                     FluidInBox[4][i][j][k] = Sigma_g*FluidInBox[0][i][j][k]
                else:
                     FluidInBox[0][i][j][k] = OutDens
                     FluidInBox[1][i][j][k] = OutMomX
                     FluidInBox[2][i][j][k] = OutMomY
                     FluidInBox[3][i][j][k] = OutMomZ
                     FluidInBox[4][i][j][k] = OutEngy

                PotInBox [i][j][k] = np.interp( r, Coarse_r, Potential )
                 

    return FluidInBox, PotInBox
