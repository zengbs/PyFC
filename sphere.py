import multiprocessing
import numpy as np
from density_profile import *
from pri2con import Pri2Con
import parameters as par
import sys
import time

par.Parameters() 

def SetIC( ):

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
    #Potential *= (par.Sigma_D/par.Const_C**2

    return GasDens,   GasMomX,  GasMomY,  GasMomZ,  GasEngy, Potential

    ####################################
    ##########      ISM     ############
    ####################################
    #PotInBoxCopy    = np.copy(PotInBox)
    #FluidInBoxCopy  = np.copy(FluidInBox)
    #Center = np.array([par.Lz,par.Ly,par.Lx])*0.5/par.CoreRadius_D
    #ISM             = NumericalISM( PotInBoxCopy, FluidInBoxCopy, PresInBox, delta, Center, FractalDensity,  FractalUxyz )
    #ISM_Temp        = ISM[4]/ISM[0]

    #ISM[0], ISM[1], ISM[2], ISM[3], ISM[4] = Pri2Con(ISM[0], ISM[1], ISM[2], ISM[3], ISM[4])
