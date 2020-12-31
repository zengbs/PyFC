from functools import partial
import multiprocessing
import numpy as np
from density_profile import *
from pri2con import Pri2Con
import parameters as par
import sys
import time

par.Parameters() 

def SphericalSphere( FractalDensity,  FractalUxyz ):

    ####################################
    ############   Fluid  ##############
    ####################################

    GasRho, TotPot  = TotPotGasDensity()
    GasPres = GasRho*par.Cs**2
    GasVelX = GasVelY = GasVelZ = np.zeros(GasRho.shape, dtype=par.Precision)

    # convert primitive to conservative variables
    GasDens,   GasMomX,  GasMomY,  GasMomZ,  GasEngy = Pri2Con( GasRho, GasVelX, GasVelY, GasVelZ, GasPres  )


    # unnormalized by `par.Sigma_D` but normalized by `par.Const_C`
    # --> since the speed of light is hard-coded to 1 in GAMER
    Potential *= (par.Sigma_D/par.Const_C)**2


    

    ######################################
    ## Filling box with fluid variables ##
    ######################################

    #EnclosedMass = 0.0
   
   
    Inside  = [ InDens,   InMomX,  InMomY,  InMomZ,  InEngy,  InPres ]
    Outside = [ OutDens, OutMomX, OutMomY, OutMomZ, OutEngy, OutPres ]

    # multi-processing generation of fluid 
    start = time.time()

    NCore = multiprocessing.cpu_count()

    p1 = multiprocessing.Pool(NCore)

    FunPartial = partial(Splitting3DFluid, delta=delta, Center=Center,
                         Coarse_r=Coarse_r, Inside=Inside, Outside=Outside)

    Pack = p1.map(FunPartial, range(int(NCore)))

    p1.close()
    p1.join()

    FluidInBox = np.concatenate( [ Pack[i][1] for i in range(NCore) ], axis=1 )
    PresInBox  = np.concatenate( [ Pack[i][2] for i in range(NCore) ], axis=0 )
    
    p2 = multiprocessing.Pool(NCore)

    FunPartial = partial(Splitting3DPot, delta=delta, Center=Center,
                         Coarse_r=Coarse_r, Potential=Potential)

    Pack = p2.map(FunPartial, range(int(NCore)))

    p2.close()
    p2.join()

    PotInBox  = np.concatenate( [ Pack[i][1] for i in range(NCore) ], axis=0 )

    end = time.time()

    print( "Elapsed time = %e with %d CPU cores"  % (end - start, NCore))
    sys.stdout.flush()
    ####################################
    ##########      ISM     ############
    ####################################
    PotInBoxCopy    = np.copy(PotInBox)
    FluidInBoxCopy  = np.copy(FluidInBox)
    Center = np.array([par.Lz,par.Ly,par.Lx])*0.5/par.CoreRadius_D
    ISM             = NumericalISM( PotInBoxCopy, FluidInBoxCopy, PresInBox, delta, Center, FractalDensity,  FractalUxyz )
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
    KT_mcSqr  /= par.Const_MolecularWeight*par.Const_AtomicMassUnit*(par.Const_C)**2  # erg
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
