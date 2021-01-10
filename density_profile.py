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
              
    X         = Idx-par.GRA_GHOST_SIZE+0.5-par.Nx*0.5
    Y         = Jdx-par.GRA_GHOST_SIZE+0.5-par.Ny*0.5
    Z         = Kdx-par.GRA_GHOST_SIZE+0.5-par.Nz*0.5
    X        *= par.delta[0]
    Y        *= par.delta[1]
    Z        *= par.delta[2]
    R         = np.sqrt(X**2+Y**2+Z**2)
    return X, Y, Z, R


# **********************************************************************
# description: halo potential
# input      : None
# output     : potential 3D array including ghost zone
# **********************************************************************
def HaloPotential(origin):
    if origin is False:
       Nx         = par.Nx+2*par.GRA_GHOST_SIZE
       Ny         = par.Ny+2*par.GRA_GHOST_SIZE
       Nz         = par.Nz+2*par.GRA_GHOST_SIZE

       X, Y, Z, R = Create3DCoordinateArray(Nx, Ny, Nz)
    else:
       R          = 0.
    Pot        = par.V_halo**2
    Pot       *= np.log(R**2+par.d_halo**2)
    return Pot


# **********************************************************************
# description: disk potential (Miyamoto-Nagai disk)
# input      : None
# output     : potential 3D array including ghost zone
# **********************************************************************
def DiskPotential(origin):
    if origin is False:
       Nx         = par.Nx+2*par.GRA_GHOST_SIZE
       Ny         = par.Ny+2*par.GRA_GHOST_SIZE
       Nz         = par.Nz+2*par.GRA_GHOST_SIZE
   
       X, Y, Z, R = Create3DCoordinateArray(Nx, Ny, Nz)
    else:
       X = Y = Z = R = 0.

    Var1       = np.sqrt(Z**2+par.b**2)
    Var2       = par.a + Var1
    Var3       = X**2 + Y**2 + Var2**2

    Pot        = -par.NEWTON_G*par.DiskMass
    Pot       /= np.sqrt(Var3)
    return Pot


# **********************************************************************
# description: bulge potential
# input      : None
# output     : potential 3D array including ghost zone
# **********************************************************************
def BulgePotential(origin):
    if origin is False:
       Nx         = par.Nx+2*par.GRA_GHOST_SIZE
       Ny         = par.Ny+2*par.GRA_GHOST_SIZE
       Nz         = par.Nz+2*par.GRA_GHOST_SIZE

       X, Y, Z, R = Create3DCoordinateArray(Nx, Ny, Nz)
    else:
       R = 0.

    Pot        = -par.NEWTON_G*par.BulgeMass
    Pot       /= R+par.d_bulge
    return Pot
  

# **********************************************************************
# description: get isothermal gas density
# input      : None
# output     : gas density 3D array but excluding ghost zone (g/cm**3)
# **********************************************************************
def TotPotGasDensity():
    TotPot           = HaloPotential (False)
    TotPot          += DiskPotential (False)
    TotPot          += BulgePotential(False)

    TotPotCenter     = HaloPotential (True)
    TotPotCenter    += DiskPotential (True)
    TotPotCenter    += BulgePotential(True)

    Nx         = par.Nx+2*par.GRA_GHOST_SIZE
    Ny         = par.Ny+2*par.GRA_GHOST_SIZE
    Nz         = par.Nz+2*par.GRA_GHOST_SIZE

    X, Y, Z, R = Create3DCoordinateArray(Nx, Ny, Nz)

    GasDensity   = np.exp(-(TotPot-TotPotCenter)/par.Eta)
    GasDensity  *= par.PeakGasMassDensity

    PotNx = TotPot.shape[0]
    PotNy = TotPot.shape[1]
    PotNz = TotPot.shape[2]

    # remove ghost zone inside `GasDensity`
    GasDensity   = GasDensity[par.GRA_GHOST_SIZE:PotNx-par.GRA_GHOST_SIZE, :, :]
    GasDensity   = GasDensity[:, par.GRA_GHOST_SIZE:PotNy-par.GRA_GHOST_SIZE, :]
    GasDensity   = GasDensity[:, :, par.GRA_GHOST_SIZE:PotNz-par.GRA_GHOST_SIZE]

    return GasDensity, TotPot


# **********************************************************************
# description: truncation
# input      : None
# output     : 
# **********************************************************************
def Truncate(GasDensity, Pot):

    DensOutSide  = np.full(GasDensity.shape, par.PeakGasMassDensity/par.DensityRatio, dtype=par.Precision)

    GasDensity   = np.where( par.PeakGasMassDensity / GasDensity > par.DensityRatio, DensOutSide, GasDensity)

    # truncate potential
    TotPotCenter     = HaloPotential (True)
    TotPotCenter    += DiskPotential (True)
    TotPotCenter    += BulgePotential(True)

    PotOutside  = TotPotCenter + par.Cs**2*np.log(par.DensityRatio)*np.full(Pot.shape, 1., dtype=par.Precision)

    Pot         = np.where( Pot - TotPotCenter > par.Cs**2*np.log(par.DensityRatio), PotOutside, Pot)

    return GasDensity, Pot
