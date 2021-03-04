import multiprocessing
import numpy as np
from density_profile import TotGasDensity
from potential import TotPotential
from fluid import Pri2Con
import parameters as par
import sys
import time


def SetIC( FractalDensity, FractalUx, FractalUy, FractalUz ):

    ####################################
    ############   Fluid  ##############
    ####################################
    print("Computing gas density ...")
    sys.stdout.flush()
    GasRho     = TotGasDensity()
    print("Computing gas density ... done !!")
    sys.stdout.flush()

    print("Computing potential ...")
    sys.stdout.flush()
    Pot        = TotPotential()
    print("Computing potential ... done !!")
    sys.stdout.flush()

    GasPres    = GasRho*par.Eta

    # convert primitive to conservative variables
    GasDens,   GasMomX,  GasMomY,  GasMomZ,  GasEngy = Pri2Con( GasRho, 0.0, 0.0, 0.0, GasPres  )

    return GasDens,   GasMomX,  GasMomY,  GasMomZ,  GasEngy, Pot
