import multiprocessing
import numpy as np
from density_profile import *
from fluid import Pri2Con
import parameters as par
import sys
import time


def SetIC( FractalDensity, FractalUxyz ):

    ####################################
    ############   Fluid  ##############
    ####################################

    GasRho,     TotPot     = TotPotGasDensity()
    TrunGasRho, TrunTotPot = Truncate(GasRho, TotPot)

    if par.dens_FractalOn:
       FractalTrunGasRho = np.where(GasRho == TrunGasRho, TrunGasRho*FractalDensity, TrunGasRho )
    else:
       FractalTrunGasRho = TrunGasRho 

    TrunGasPres        = TrunGasRho*par.Eta

    if par.Uxyz_FractalOn:
       print("Not support yet!!")
       exit(0)
    else:
       GasVelX = GasVelY = GasVelZ = np.zeros(GasRho.shape, dtype=par.Precision)

    # convert primitive to conservative variables
    GasDens,   GasMomX,  GasMomY,  GasMomZ,  GasEngy = Pri2Con( FractalTrunGasRho, GasVelX, GasVelY, GasVelZ, TrunGasPres  )

    return GasDens,   GasMomX,  GasMomY,  GasMomZ,  GasEngy, TotPot
