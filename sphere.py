import numpy as np
from density_profile import *
from pri2con import Pri2Con
import parameters as par


def SphericalSphere( ):

    par.Parameters() 

    # Center of sphere
    Center = np.array([par.Lx,par.Ly,par.Lz])*0.5/par.Radius_D
 
    # Coarse grid (Normalized by Radius_D)
    Coarse_r = np.arange(par.BPoint, 0.5*par.Lx/par.Radius_D, par.CoarseDr)

    # Left edge and right edge coordinates of the desired
    # simulation domain which will be used in GAMER. (Normalized by Radius_D)
    le = np.array([ 0,  0,  0])
    re = np.array([par.Lx, par.Ly, par.Lz]) / par.Radius_D
    
    # Cell spacing (Normalized by Radius_D)
    N = [par.Nx, par.Ny, par.Nz ]
    delta = (re-le)/N

    ####################################
    ############   Fluid  ##############
    ####################################

    # Fluid inside box
    InRho, Nothing   = NumericalDensity( Coarse_r )
    InUX    = 0 
    InUy    = 0
    InUz    = 0
    InPres  = InRho*(par.Sigma_g/par.C)**2

    # Fluid outside box
    OutRho  = 1e-3*np.interp( par.Radius/par.Radius_D, Coarse_r, InRho )
    OutUX   = 0 
    OutUy   = 0
    OutUz   = 0
    OutPres = 1e3*OutRho*(par.Sigma_g/par.C)**2

    # Conversion
    InDens,   InMomX,  InMomY,  InMomZ,  InEngy = Pri2Con(  InRho,  InUX,  InUy,  InUz, InPres  )
    OutDens, OutMomX, OutMomY, OutMomZ, OutEngy = Pri2Con( OutRho, OutUX, OutUy, OutUz, OutPres )

    
    FluidInBox = np.zeros((5, par.Nx, par.Ny, par.Nz), dtype=par.Precision)


    ####################################
    ##########   Potential  ############
    ####################################
    GRA_GHOST_SIZE = 2
    PotInBox   = np.zeros((par.Nx+2*GRA_GHOST_SIZE, par.Ny+2*GRA_GHOST_SIZE, par.Nz+2*GRA_GHOST_SIZE), dtype=par.Precision)

    Psi0       = par.Phi0     / par.Sigma_D / par.Sigma_D
    DevPsi0    = par.DevPhi0  / par.Sigma_D / par.Sigma_D

    Potential = NumericalTotalPotential( Coarse_r, par.Kappa, par.Lambda, par.Constant, Psi0, DevPsi0 )
    Potential *= (par.Sigma_D/par.C)**2

    EnclosedMass = 0.0


    for i in range(par.Nx+2*GRA_GHOST_SIZE):
        for j in range(par.Ny+2*GRA_GHOST_SIZE):
            for k in range(par.Nz+2*GRA_GHOST_SIZE):
                ii = i - GRA_GHOST_SIZE                     
                jj = j - GRA_GHOST_SIZE                     
                kk = k - GRA_GHOST_SIZE                     

                x = (i+0.5)*delta[0] 
                y = (j+0.5)*delta[1]
                z = (k+0.5)*delta[2]

                r = np.sqrt((x-Center[0])**2 + (y-Center[1])**2 + (z-Center[2])**2)

                if ( 0<=ii<par.Nx and 0<=jj<par.Ny and 0<=kk<par.Nz ):
                     if ( r < par.Radius/par.Radius_D ):
                          FluidInBox[0][ii][jj][kk] = np.interp( r, Coarse_r, InDens )
                          FluidInBox[1][ii][jj][kk] = np.interp( r, Coarse_r, InMomX ) 
                          FluidInBox[2][ii][jj][kk] = np.interp( r, Coarse_r, InMomY ) 
                          FluidInBox[3][ii][jj][kk] = np.interp( r, Coarse_r, InMomZ ) 
                          FluidInBox[4][ii][jj][kk] = np.interp( r, Coarse_r, InEngy )
                          EnclosedMass += np.interp( r, Coarse_r, InRho )
                     else:
                          FluidInBox[0][ii][jj][kk] = OutDens
                          FluidInBox[1][ii][jj][kk] = OutMomX
                          FluidInBox[2][ii][jj][kk] = OutMomY
                          FluidInBox[3][ii][jj][kk] = OutMomZ
                          FluidInBox[4][ii][jj][kk] = OutEngy

                PotInBox [i][j][k] = np.interp( r, Coarse_r, Potential )


    MeanRho         = 3*EnclosedMass/( 4*np.pi*par.Radius**3 )
    FreeFallingTime = 0.5427/( par.NEWTON_G * MeanRho )**0.5

    print("Encloed mass = %e\n" %(EnclosedMass))
    print("Free falling time = %e\n" % (FreeFallingTime))

    return FluidInBox, PotInBox
