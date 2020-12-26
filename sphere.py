from functools import partial
import multiprocessing
import numpy as np
from density_profile import *
from pri2con import Pri2Con
import parameters as par
import sys

par.Parameters() 

def Splitting3DFluid(CoreIdx, delta, Center, Coarse_r, Potential, Inside, Outside):

    NCore = multiprocessing.cpu_count()
    TotalLoading = par.Nz

    if TotalLoading % NCore is not 0:
       print("TotalLoading % NCore is not 0!!")
       exit(0)

    LoadingPerCore  = TotalLoading / NCore

    FluidInBox_Rank = np.zeros((5, LoadingPerCore, par.Ny, par.Nx), dtype=par.Precision)
    PresInBox_Rank  = np.zeros((   LoadingPerCore, par.Ny, par.Nx), dtype=par.Precision)
    PotInBox_Rank   = np.zeros((   LoadingPerCore, par.Ny, par.Nx), dtype=par.Precision)


    InDens,   InMomX,  InMomY,  InMomZ,  InEngy,  InPres = Inside
    OutDens, OutMomX, OutMomY, OutMomZ, OutEngy, OutPres = Outside

    for k in range(LoadingPerCore):
        for j in range(par.Ny):
            for i in range(par.Nx):

                x = (i+0.5)*delta[0] 
                y = (j+0.5)*delta[1]
                z = (CoreIdx*LoadingPerCore+k+0.5)*delta[2]

                r = np.sqrt((x-Center[0])**2 + (y-Center[1])**2 + (z-Center[2])**2)

                #filling `FluidInBox` with fluid variables
                if ( 0<=i<par.Nx and 0<=j<par.Ny and 0<=k<par.Nz ):
                     if ( r < par.SphereRadius/par.CoreRadius_D ):
                          FluidInBox_Rank[0][k][j][i] = np.interp( r, Coarse_r, InDens )
                          FluidInBox_Rank[1][k][j][i] = np.interp( r, Coarse_r, InMomX ) 
                          FluidInBox_Rank[2][k][j][i] = np.interp( r, Coarse_r, InMomY ) 
                          FluidInBox_Rank[3][k][j][i] = np.interp( r, Coarse_r, InMomZ ) 
                          FluidInBox_Rank[4][k][j][i] = np.interp( r, Coarse_r, InEngy )
                          PresInBox_Rank    [k][j][i] = np.interp( r, Coarse_r, InPres )
                     else:
                          FluidInBox_Rank[0][k][j][i] = OutDens
                          FluidInBox_Rank[1][k][j][i] = OutMomX
                          FluidInBox_Rank[2][k][j][i] = OutMomY
                          FluidInBox_Rank[3][k][j][i] = OutMomZ
                          FluidInBox_Rank[4][k][j][i] = OutEngy
                          PresInBox_Rank    [k][j][i] = OutPres

                PotInBox_Rank [k][j][i] = np.interp( r, Coarse_r, Potential )
      
    return [CoreIdx, FluidInBox_Rank, PresInBox_Rank, PotInBox_Rank]



def SphericalSphere( Fractal ):


    # Center of sphere (normalized by `par.CoreRadius_D`)
    Center = np.array([par.Lz,par.Ly,par.Lx])*0.5/par.CoreRadius_D
 
    # Coarse grid (Normalized by CoreRadius_D)
    Coarse_r = np.arange(par.BPoint, np.sqrt(2)*par.Lx/par.CoreRadius_D, par.CoarseDr)

    # Left edge and right edge coordinates of the desired
    # simulation domain which will be used in GAMER. (Normalized by CoreRadius_D)
    le = np.array([ 0.0,  0.0,  0.0])
    re = np.array([par.Lz, par.Ly, par.Lx]) / float(par.CoreRadius_D)
    
    # Cell spacing (normalized by CoreRadius_D)
    N = np.array([par.Nz, par.Ny, par.Nx ], dtype=par.Precision)
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
    if par.Case == "Mukherjee":
       SoubdSpeedSqr  = par.Const_kB*par.Temp_g/(par.Const_MeanMolecularWeight*par.Const_AtomMass*par.Const_Erg2eV)
       InPres  = InRho * SoubdSpeedSqr / par.Const_C**2
    if par.Case == "Standard":
       # `InPres` is nunormalized by `Sigma_g` but normalized by `par.Const_C`
       # --> since the speed of light is hard-coded to 1 in GAMER
       InPres  = InRho*(par.Sigma_g/par.Const_C)**2


    # Fluid outside box
    OutRho  = 1e-3*np.interp( par.SphereRadius/par.CoreRadius_D, Coarse_r, InRho )
    OutUX   = 0 
    OutUy   = 0
    OutUz   = 0
    # unnormalized by `par.Sigma_g` but normalized by `par.Const_C`
    # --> since the speed of light is hard-coded to 1 in GAMER
    if par.Case == "Mukherjee":
       OutPres = 1e3*OutRho*SoubdSpeedSqr / par.Const_C**2
    if par.Case == "Standard":
       OutPres = 1e3*OutRho*(par.Sigma_g/par.Const_C)**2

    # Conversion
    InDens,   InMomX,  InMomY,  InMomZ,  InEngy = Pri2Con(  InRho,  InUX,  InUy,  InUz, InPres  )
    OutDens, OutMomX, OutMomY, OutMomZ, OutEngy = Pri2Con( OutRho, OutUX, OutUy, OutUz, OutPres )


    

    ####################################
    ##########   Potential  ############
    ####################################
    GRA_GHOST_SIZE = par.GRA_GHOST_SIZE


    Psi0       = par.Phi0    / par.Sigma_D**2 
    DevPsi0    = par.DevPhi0 / par.Sigma_D**2

    Potential = NumericalTotalPotential( Coarse_r, Psi0, DevPsi0 )

    # unnormalized by `par.Sigma_D` but normalized by `par.Const_C`
    # --> since the speed of light is hard-coded to 1 in GAMER
    Potential *= (par.Sigma_D/par.Const_C)**2


    

    ######################################
    ## Filling box with fluid variables ##
    ######################################

    #EnclosedMass = 0.0
   
   
    Inside  = [ InDens,   InMomX,  InMomY,  InMomZ,  InEngy,  InPres ]
    Outside = [ OutDens, OutMomX, OutMomY, OutMomZ, OutEngy, OutPres ]

    # multi-processing
    NCore = multiprocessing.cpu_count()

    p = multiprocessing.Pool(NCore)

    FunPartial = partial(Splitting3DFluid, delta=delta, Center=Center,
                         Coarse_r=Coarse_r, Potential=Potential,
                         Inside=Inside, Outside=Outside)

    Pack = p.map(FunPartial, range(int(NCore)))

    p.close()
    p.join()

    FluidInBox = np.concatenate( [ Pack[i][1] for i in range(NCore) ], axis=1 )
    PresInBox  = np.concatenate( [ Pack[i][2] for i in range(NCore) ], axis=0 )
    PotInBox   = np.concatenate( [ Pack[i][3] for i in range(NCore) ], axis=0 )
    


    ####################################
    ##########      ISM     ############
    ####################################
    PotInBoxCopy    = np.copy(PotInBox)
    FluidInBoxCopy  = np.copy(FluidInBox)
    ISM             = NumericalISM( PotInBoxCopy, FluidInBoxCopy, PresInBox, delta, Center, Fractal )
    ISM_Temp        = ISM[4]/ISM[0]

    ISM[0], ISM[1], ISM[2], ISM[3], ISM[4] = Pri2Con(ISM[0], ISM[1], ISM[2], ISM[3], ISM[4])

    ###################################
    ######     
    ###################################
    Idx = np.indices((par.Nz, par.Ny, par.Nx))[2]
    Jdx = np.indices((par.Nz, par.Ny, par.Nx))[1]
    Kdx = np.indices((par.Nz, par.Ny, par.Nx))[0]
        
    X   = (Idx+0.5)*delta[2]-Center[2]
    Y   = (Jdx+0.5)*delta[1]-Center[1]
    Z   = (Kdx+0.5)*delta[0]-Center[0] 
    R   = np.sqrt(X**2+Y**2+Z**2) 

    KT_mcSqr   = par.CriticalTemp*par.Const_kB                                      # K* (erg/K)
    KT_mcSqr  /= par.Const_MeanMolecularWeight*par.Const_AtomMass*(par.Const_C)**2  # erg
    KT_mcSqr  /= par.Const_Erg2eV

    ISM        = np.where( ISM_Temp < KT_mcSqr, ISM, FluidInBox )

    FluidInBox = np.where( R<par.SphereRadius/par.CoreRadius_D, ISM, FluidInBox ) 

    #####################################################
    ### Calculate enclosed mass and free-falling time ###
    #####################################################

    #MeanRho         = 3*EnclosedMass/( 4*np.pi*par.SphereRadius**3 )
    #FreeFallingTime = 0.5427/( par.NEWTON_G * MeanRho )**0.5

    #print("Encloed mass = %e\n"      % (EnclosedMass)   )
    #print("Free falling time = %e\n" % (FreeFallingTime))
    return FluidInBox, PotInBox
