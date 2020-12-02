import numpy as np
from density_profile import *
from pri2con import Pri2Con


def SphericalSphere( L, N, ParaPhy, ParaNum, Precision ):

    # Speed of light (km/s)
    C = 299792.458

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
    InPres  = InRho*(Sigma_g/C)**2

    # Fluid outside box
    OutRho  = 1e-3*np.interp( PlotRadius, Coarse_r, InRho )
    OutUX   = 0 
    OutUy   = 0
    OutUz   = 0
    OutPres = 1e3*OutRho*(Sigma_g/C)**2

    # Conversion
    InDens,   InMomX,  InMomY,  InMomZ,  InEngy = Pri2Con(  InRho,  InUX,  InUy,  InUz, InPres  )
    OutDens, OutMomX, OutMomY, OutMomZ, OutEngy = Pri2Con( OutRho, OutUX, OutUy, OutUz, OutPres )

    
    FluidInBox = np.zeros((5, Nx, Ny, Nz), dtype=Precision)


    ####################################
    ##########   Potential  ############
    ####################################
    GRA_GHOST_SIZE = 2
    PotInBox   = np.zeros((Nx+2*GRA_GHOST_SIZE, Ny+2*GRA_GHOST_SIZE, Nz+2*GRA_GHOST_SIZE), dtype=Precision)

    Rho0_D, Sigma_D, Radius_D  = Free2DerivedPara( Rho0_g, Sigma_g, Radius_g, Lambda, Kappa )
    Psi0       = Phi0      / Sigma_D / Sigma_D
    DevPsi0    = DevPhi0   / Sigma_D / Sigma_D

    Potential = NumericalTotalPotential( Coarse_r, Kappa, Lambda, Constant, Psi0, DevPsi0 )
    Potential *= (Sigma_D/C)**2



    for i in range(Nx+2*GRA_GHOST_SIZE):
        for j in range(Ny+2*GRA_GHOST_SIZE):
            for k in range(Nz+2*GRA_GHOST_SIZE):

                ii = i - GRA_GHOST_SIZE                     
                jj = j - GRA_GHOST_SIZE                     
                kk = k - GRA_GHOST_SIZE                     

                x = (i+0.5)*delta[0] 
                y = (j+0.5)*delta[1]
                z = (k+0.5)*delta[2]

                r = np.sqrt((x-Center[0])**2 + (y-Center[1])**2 + (z-Center[2])**2)

                if ( 0<=ii<Nx and 0<=jj<Ny and 0<=kk<Nz ):
                     if ( r < PlotRadius ):
                          FluidInBox[0][ii][jj][kk] = np.interp( r, Coarse_r, InDens )
                          FluidInBox[1][ii][jj][kk] = np.interp( r, Coarse_r, InMomX ) 
                          FluidInBox[2][ii][jj][kk] = np.interp( r, Coarse_r, InMomY ) 
                          FluidInBox[3][ii][jj][kk] = np.interp( r, Coarse_r, InMomZ ) 
                          FluidInBox[4][ii][jj][kk] = np.interp( r, Coarse_r, InEngy )
                     else:
                          FluidInBox[0][ii][jj][kk] = OutDens
                          FluidInBox[1][ii][jj][kk] = OutMomX
                          FluidInBox[2][ii][jj][kk] = OutMomY
                          FluidInBox[3][ii][jj][kk] = OutMomZ
                          FluidInBox[4][ii][jj][kk] = OutEngy

                PotInBox [i][j][k] = np.interp( r, Coarse_r, Potential )
                 

    return FluidInBox, PotInBox
