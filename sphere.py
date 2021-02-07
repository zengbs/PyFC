import multiprocessing
import numpy as np
from density_profile import *
from fluid import Pri2Con
import parameters as par
import sys
import time


def SetIC( FractalDensity, FractalUx, FractalUy, FractalUz ):

    ####################################
    ############   Fluid  ##############
    ####################################

    GasRho,     TotPot     = TotPotGasDensity()
    TrunGasRho, TrunTotPot, TrunFractalUx, TrunFractalUy, TrunFractalUz = Truncate(GasRho, TotPot, FractalUx, FractalUy, FractalUz )

    FractalTrunGasRho = np.where(GasRho == TrunGasRho, TrunGasRho*FractalDensity, TrunGasRho )

    TrunGasPres        = TrunGasRho*par.Eta

    # convert primitive to conservative variables
    GasDens,   GasMomX,  GasMomY,  GasMomZ,  GasEngy = Pri2Con( FractalTrunGasRho, TrunFractalUx, TrunFractalUy, TrunFractalUz, TrunGasPres  )

    return GasDens,   GasMomX,  GasMomY,  GasMomZ,  GasEngy, TotPot
