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

    InsideFracTrunGasRho = np.where( par.PeakGasMassDensity / TrunGasRho > par.FracDensityRatio, False, True)

    FracTrunGasRho       = np.where( InsideFracTrunGasRho, TrunGasRho*FractalDensity, TrunGasRho )
 
    InsideFracTrunGasRho = np.logical_and( InsideFracTrunGasRho, FracTrunGasRho > 1e-28/par.UNIT_D )

    FracTrunGasRho       = np.where( InsideFracTrunGasRho, FracTrunGasRho, GasRho )

    TrunGasPres          = GasRho*par.Eta

    # convert primitive to conservative variables
    GasDens,   GasMomX,  GasMomY,  GasMomZ,  GasEngy = Pri2Con( FracTrunGasRho, TrunFractalUx, TrunFractalUy, TrunFractalUz, TrunGasPres  )

    return GasDens,   GasMomX,  GasMomY,  GasMomZ,  GasEngy, TotPot
