import sys
import numpy as np
import parameters as par


# **********************************************************************
# description: create four 3D array storing 
#              (1) x-coordinate
#              (2) y-coordinate
#              (3) z-coordinate
#              (4) the distance between cell itself and center
# input      : None
# output     : 3D array stroing x, y, z, and radial distance
# **********************************************************************
def Create3DCoordinateArray(Nx, Ny, Nz):
    Idx       = np.indices((Nx, Ny, Nz), dtype=par.Precision)[1]
    Jdx       = np.indices((Nx, Ny, Nz), dtype=par.Precision)[2]
    Kdx       = np.indices((Nx, Ny, Nz), dtype=par.Precision)[0]
              
    X         = Idx+0.5-par.Nx*0.5
    Y         = Jdx+0.5-par.Ny*0.5
    Z         = Kdx+0.5-par.Nz*0.5
    X        *= par.delta[0]
    Y        *= par.delta[1]
    Z        *= par.delta[2]
    R         = np.sqrt(X**2+Y**2+Z**2)
    return X, Y, Z, R


  

# **********************************************************************
# description: get isothermal gas density
# input      : None
# output     : gas density 3D array but excluding ghost zone (g/cm**3)
# **********************************************************************
def GasDensity(FractalDensity):
    X, Y, Z, R = Create3DCoordinateArray(par.Nx, par.Ny, par.Nz)

    AmbientDens = np.full(FractalDensity.shape, par.AmbientDens, dtype=par.Precision)

    FractalDensity = FractalDensity*AmbientDens

    InsideFractal = np.logical_and(-0.15*par.Lz<Z, Z<0.15*par.Lz)
    GasDensity = np.where(InsideFractal, FractalDensity , AmbientDens)


    #GasDensity = np.where( GasDensity > AmbientDens*par.dens_mean, GasDensity, AmbientDens )

    GasPres = par.AmbientDens*par.Eta

    FractalTemperature = GasPres / FractalDensity

    FractalTemperature *= par.Const_AtomicMassUnit*par.Const_MolecularWeight*par.Const_C**2/par.Const_kB # (Kelvin)

#    GasDensity = np.where( FractalTemperature < par.Temperature, GasDensity, AmbientDens )

    return GasDensity
