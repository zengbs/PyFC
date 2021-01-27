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

    GasRho  = GasDensity(FractalDensity)

    GasPres = par.AmbientDens*par.Eta

    GasVelX = GasVelY = GasVelZ = 0

    # convert primitive to conservative variables
    GasDens,   GasMomX,  GasMomY,  GasMomZ,  GasEngy = Pri2Con( GasRho, GasVelX, GasVelY, GasVelZ, GasPres  )

    return GasDens,   GasMomX,  GasMomY,  GasMomZ,  GasEngy
