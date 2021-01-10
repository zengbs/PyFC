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

    FractalTrunGasRho = np.where(GasRho == TrunGasRho, TrunGasRho*FractalDensity, TrunGasRho )

    TrunGasPres        = TrunGasRho*par.Eta

    GasVelX = GasVelY = GasVelZ = np.zeros(GasRho.shape, dtype=par.Precision)

    # convert primitive to conservative variables
    GasDens,   GasMomX,  GasMomY,  GasMomZ,  GasEngy = Pri2Con( FractalTrunGasRho, GasVelX, GasVelY, GasVelZ, TrunGasPres  )

    return GasDens,   GasMomX,  GasMomY,  GasMomZ,  GasEngy, TotPot
