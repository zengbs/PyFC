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

    GasRho     = TotGasDensity()

    GasPres    = GasRho*par.Eta

    # convert primitive to conservative variables
    GasDens,   GasMomX,  GasMomY,  GasMomZ,  GasEngy = Pri2Con( GasRho, 0.0, 0.0, 0.0, GasPres  )

    return GasDens,   GasMomX,  GasMomY,  GasMomZ,  GasEngy
