import numpy as np
from density_profile import *
from pri2con import Pri2Con
import parameters as par


def SphericalSphere( ):

    par.Parameters() 

    # Center of sphere (normalized by `par.CoreRadius_D`)
    Center = np.array([par.Lz,par.Ly,par.Lx])*0.5/par.CoreRadius_D
 
    # Coarse grid (Normalized by CoreRadius_D)
    Coarse_r = np.arange(par.BPoint, np.sqrt(2)*par.Lx/par.CoreRadius_D, par.CoarseDr)

    # Left edge and right edge coordinates of the desired
    # simulation domain which will be used in GAMER. (Normalized by CoreRadius_D)
    le = np.array([ 0,  0,  0])
    re = np.array([par.Lz, par.Ly, par.Lx]) / par.CoreRadius_D
    
    # Cell spacing (normalized by CoreRadius_D)
    N = [par.Nz, par.Ny, par.Nx ]
    delta = (re-le)/N

    ####################################
    ############   Fluid  ##############
    ####################################

    # Fluid inside box
    InRho, Nothing   = NumericalDensity( Coarse_r )
    InUX    = 0 
    InUy    = 0
    InUz    = 0
    # unnormalized by `par.Sigma_g` but normalized by `par.Const_C`
    # --> since the speed of light is hard-coded to 1 in GAMER
    InPres  = InRho * par.Const_kB*par.Temp_g / (par.Const_MeanMolecularWeight*par.Const_AtomMass*par.Const_Erg2eV) / (par.Const_C*1e5)**2
    InPres  = InRho*(par.Sigma_g/par.Const_C)**2


    # Fluid outside box
    OutRho  = 1e-3*np.interp( par.SphereRadius/par.CoreRadius_D, Coarse_r, InRho )
    OutUX   = 0 
    OutUy   = 0
    OutUz   = 0
    # unnormalized by `par.Sigma_g` but normalized by `par.Const_C`
    # --> since the speed of light is hard-coded to 1 in GAMER
    OutPres = 1e3*OutRho*(par.Sigma_g/par.Const_C)**2

    # Conversion
    InDens,   InMomX,  InMomY,  InMomZ,  InEngy = Pri2Con(  InRho,  InUX,  InUy,  InUz, InPres  )
    OutDens, OutMomX, OutMomY, OutMomZ, OutEngy = Pri2Con( OutRho, OutUX, OutUy, OutUz, OutPres )


    

    ####################################
    ##########   Potential  ############
    ####################################
    GRA_GHOST_SIZE = par.GRA_GHOST_SIZE

    PotInBox   = np.zeros((par.Nz+2*GRA_GHOST_SIZE, par.Ny+2*GRA_GHOST_SIZE, par.Nx+2*GRA_GHOST_SIZE), dtype=par.Precision)

    Psi0       = par.Phi0    / par.Sigma_D**2 
    DevPsi0    = par.DevPhi0 / par.Sigma_D**2

    Potential = NumericalTotalPotential( Coarse_r, Psi0, DevPsi0 )

    # unnormalized by `par.Sigma_D` but normalized by `par.Const_C`
    # --> since the speed of light is hard-coded to 1 in GAMER
    Potential *= (par.Sigma_D/par.Const_C)**2


    

    ######################################
    ## Filling box with fluid variables ##
    ######################################

    EnclosedMass = 0.0
    FluidInBox = np.zeros((5, par.Nz, par.Ny, par.Nx), dtype=par.Precision)
    PresInBox  = np.zeros((par.Nz, par.Ny, par.Nx),    dtype=par.Precision)


    for k in range(par.Nz+2*GRA_GHOST_SIZE):
        for j in range(par.Ny+2*GRA_GHOST_SIZE):
            for i in range(par.Nx+2*GRA_GHOST_SIZE):
                # indices for non-ghost zone
                ii = i - GRA_GHOST_SIZE                     
                jj = j - GRA_GHOST_SIZE                     
                kk = k - GRA_GHOST_SIZE                     

                x = (ii+0.5)*delta[0] 
                y = (jj+0.5)*delta[1]
                z = (kk+0.5)*delta[2]

                r = np.sqrt((x-Center[0])**2 + (y-Center[1])**2 + (z-Center[2])**2)

                #filling `FluidInBox` with fluid variables
                if ( 0<=ii<par.Nx and 0<=jj<par.Ny and 0<=kk<par.Nz ):
                     if ( r < par.SphereRadius/par.CoreRadius_D ):
                          FluidInBox[0][kk][jj][ii] = np.interp( r, Coarse_r, InDens )
                          FluidInBox[1][kk][jj][ii] = np.interp( r, Coarse_r, InMomX ) 
                          FluidInBox[2][kk][jj][ii] = np.interp( r, Coarse_r, InMomY ) 
                          FluidInBox[3][kk][jj][ii] = np.interp( r, Coarse_r, InMomZ ) 
                          FluidInBox[4][kk][jj][ii] = np.interp( r, Coarse_r, InEngy )
                          PresInBox    [kk][jj][ii] = np.interp( r, Coarse_r, InPres )
                          EnclosedMass += np.interp( r, Coarse_r, InRho )
                     else:
                          FluidInBox[0][kk][jj][ii] = OutDens
                          FluidInBox[1][kk][jj][ii] = OutMomX
                          FluidInBox[2][kk][jj][ii] = OutMomY
                          FluidInBox[3][kk][jj][ii] = OutMomZ
                          FluidInBox[4][kk][jj][ii] = OutEngy
                          PresInBox    [kk][jj][ii] = OutPres
 
      

                # filling `PotInBox` with Potential
                PotInBox [k][j][i] = np.interp( r, Coarse_r, Potential )
    

    ####################################
    ##########      ISM     ############
    ####################################
    #PotInBoxCopy    = np.copy(PotInBox)
    #FluidInBoxCopy  = np.copy(FluidInBox)
    #ISM             = NumericalISM( PotInBoxCopy, FluidInBoxCopy, PresInBox, delta, Center )
    #ISM_Temp        = ISM[4]/ISM[0]

    #ISM[0], ISM[1], ISM[2], ISM[3], ISM[4] = Pri2Con(ISM[0], ISM[1], ISM[2], ISM[3], ISM[4])

    ####################################
    #######     
    ####################################
    #Idx = np.indices((par.Nz, par.Ny, par.Nx))[2]
    #Jdx = np.indices((par.Nz, par.Ny, par.Nx))[1]
    #Kdx = np.indices((par.Nz, par.Ny, par.Nx))[0]
    #    
    #X   = (Idx+0.5)*delta[2]-Center[2]
    #Y   = (Jdx+0.5)*delta[1]-Center[1]
    #Z   = (Kdx+0.5)*delta[0]-Center[0] 
    #R   = np.sqrt(X**2+Y**2+Z**2) 

    #KT_mcSqr  = par.CriticalTemp*par.Const_kB          # K* (erg/K)
    #KT_mcSqr /= par.Const_MeanMolecularWeight*par.Const_AtomMass*(par.Const_C*1e5)**2  # erg
    #KT_mcSqr /= par.Const_Erg2eV

    #ISM[0] = np.where( ISM_Temp < KT_mcSqr, ISM[0], FluidInBox[0] )
    #ISM[1] = np.where( ISM_Temp < KT_mcSqr, ISM[1], FluidInBox[1] )
    #ISM[2] = np.where( ISM_Temp < KT_mcSqr, ISM[2], FluidInBox[2] )
    #ISM[3] = np.where( ISM_Temp < KT_mcSqr, ISM[3], FluidInBox[3] )
    #ISM[4] = np.where( ISM_Temp < KT_mcSqr, ISM[4], FluidInBox[4] )


    #FluidInBox[0] = np.where( R<par.SphereRadius/par.CoreRadius_D, ISM[0], FluidInBox[0] ) 
    #FluidInBox[1] = np.where( R<par.SphereRadius/par.CoreRadius_D, ISM[1], FluidInBox[1] ) 
    #FluidInBox[2] = np.where( R<par.SphereRadius/par.CoreRadius_D, ISM[2], FluidInBox[2] ) 
    #FluidInBox[3] = np.where( R<par.SphereRadius/par.CoreRadius_D, ISM[3], FluidInBox[3] ) 
    #FluidInBox[4] = np.where( R<par.SphereRadius/par.CoreRadius_D, ISM[4], FluidInBox[4] ) 

    ####################################
    #######     
    ####################################

    MeanRho         = 3*EnclosedMass/( 4*np.pi*par.SphereRadius**3 )
    FreeFallingTime = 0.5427/( par.NEWTON_G * MeanRho )**0.5

    print("Encloed mass = %e\n"      % (EnclosedMass)   )
    print("Free falling time = %e\n" % (FreeFallingTime))
    return FluidInBox, PotInBox
