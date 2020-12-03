import numpy as np
from pri2con import Pri2Con
from sphere import SphericalSphere
from density_profile import Free2DerivedPara
import parameters as par

def _DumpFile():
    par.Parameters()

    ############################
    ###     Dump density    ####
    ############################

    # Number of cells along each dimension of the input grid.
    N = np.array([par.Nx, par.Ny, par.Nz], dtype='int')
    L = np.array([par.Lx, par.Ly, par.Lz], dtype=par.Precision) 


    FluidInBox, PotInBox = SphericalSphere()

    FluidInBox.tofile("UM_IC")

    ############################
    ###     Dump potential   ###
    ############################

    PotInBox.tofile("ExtPotTable")

par.Parameters()
out=np.zeros((par.Nx, par.Ny, par.Nz))
_DumpFile()
