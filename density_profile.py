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
    Idx       = np.indices((Nz, Ny, Nx), dtype=par.Precision)[2]
    Jdx       = np.indices((Nz, Ny, Nx), dtype=par.Precision)[1]
    Kdx       = np.indices((Nz, Ny, Nx), dtype=par.Precision)[0]
              
    X         = (Idx+0.5)*par.delta[2]-par.Center[2]
    Y         = (Jdx+0.5)*par.delta[1]-par.Center[1]
    Z         = (Kdx+0.5)*par.delta[0]-par.Center[0] 
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


    GasDensity   = np.exp(-(TotPot-TotPotCenter)/par.Cs**2)
    GasDensity  *= par.PeakGasMassDensity
 
    PotNx = TotPot.shape[0]
    PotNy = TotPot.shape[1]
    PotNz = TotPot.shape[2]

    # remove ghost zone inside `GasDensity`
    GasDensity   = GasDensity[par.GRA_GHOST_SIZE:PotNx-par.GRA_GHOST_SIZE, :, :]
    GasDensity   = GasDensity[:, par.GRA_GHOST_SIZE:PotNy-par.GRA_GHOST_SIZE, :]
    GasDensity   = GasDensity[:, :, par.GRA_GHOST_SIZE:PotNz-par.GRA_GHOST_SIZE]

    return GasDensity, TotPot


def Truncate(GasDensity, Pot):
    # truncate density
    DensityRatio = 500.

    DensOutSide  = np.full(GasDensity.shape, par.PeakGasMassDensity/DensityRatio, dtype=par.Precision)

    GasDensity   = np.where( par.PeakGasMassDensity / GasDensity > DensityRatio, DensOutSide, GasDensity)

    # truncate potential
    TotPotCenter     = HaloPotential (True)
    TotPotCenter    += DiskPotential (True)
    TotPotCenter    += BulgePotential(True)

    PotOutside  = TotPotCenter + par.Cs**2*np.log(DensityRatio)*np.full(Pot.shape, 1., dtype=par.Precision)

    Pot         = np.where( Pot - TotPotCenter > par.Cs**2*np.log(DensityRatio), PotOutside, Pot)

    return GasDensity, Pot

# **********************************************************************
# description: 
# input      :
# output     : 
# **********************************************************************
def NumericalISM( PotInBox, FluidInBox, PresInBox, FractalDensity, FractalUxyz ):

    # create an array stored ISM
    ISM       = np.zeros((5, par.Nz, par.Ny, par.Nx), dtype=par.Precision)

    # unnormalized by `par.Const_C**2`
    PotInBox *= par.Const_C**2

    # the box storing potential includes ghost zone
    PotNz     = par.Nz + 2*par.GRA_GHOST_SIZE
    PotNy     = par.Ny + 2*par.GRA_GHOST_SIZE
    PotNx     = par.Nx + 2*par.GRA_GHOST_SIZE

    # remove ghost zone inside `PotInBox`
    PotInBox  = PotInBox[par.GRA_GHOST_SIZE:PotNx-par.GRA_GHOST_SIZE, :, :]
    PotInBox  = PotInBox[:, par.GRA_GHOST_SIZE:PotNy-par.GRA_GHOST_SIZE, :]
    PotInBox  = PotInBox[:, :, par.GRA_GHOST_SIZE:PotNz-par.GRA_GHOST_SIZE]

    # extract the potential values on the equator
    if par.Nz%2 == 0:
       PotOnEquator = 0.5*( PotInBox[:,:,int(par.Nz/2-1)] + PotInBox[:,:,int(par.Nz/2)] )
    else:
       PotOnEquator = PotInBox[:,:,int((par.Nz-1)*0.5)]

    # Construct an array by repeating `PotOnEquator` along z-direction
    PotInBoxExtendZ = np.tile(PotOnEquator,(par.Nz,1,1))

    # Construct indices array
    Idx = np.indices((par.Nz, par.Ny, par.Nx))[2]
    Jdx = np.indices((par.Nz, par.Ny, par.Nx))[1]
    Kdx = np.indices((par.Nz, par.Ny, par.Nx))[0]

    X   = (Idx+0.5)*par.delta[2]-par.Center[2]
    Y   = (Jdx+0.5)*par.delta[1]-par.Center[1]
    Z   = (Kdx+0.5)*par.delta[0]-par.Center[0] 
    R   = np.sqrt(X**2+Y**2)


    # expression below have assumed that the potential at the center of sphere is zero
    if par.Case == "Mukherjee":
       ISM[0] = par.ISM0   * np.exp( -( PotInBox -  PotInBoxExtendZ*par.Epsilon**2 )/par.Sigma_t**2 )
    if par.Case == "Standard":
       a = par.a0 * np.exp(-np.abs(Z)/par.z0)
       ISM[0] = par.PeakNumberDensity * np.exp( -( PotInBox -  PotInBoxExtendZ*a**2           )/par.Sigma_t**2 )

#####################
#    CosTheta = X/R
#    SinTheta = Y/R
#    fig, ax = plt.subplots(1,1)
#    #ax[0].imshow(np.flipud(X[int(par.Nz/2),:,:]), interpolation="None")
#    #ax[1].imshow(np.flipud(Y[int(par.Nz/2),:,:]), interpolation="None")
#    #ax[2].imshow(np.flipud(CosTheta[int(par.Nz/2),:,:]), interpolation="None")
#    cbr=ax.imshow(PresInBox[int(par.Nz/2),:,:], interpolation="None")
#    #ax[1].imshow(Y[int(par.Nz/2),:,:], interpolation="None")
#    #ax[2].imshow(CosTheta[int(par.Nz/2),:,:], interpolation="None")
#    #cbr=ax[3].imshow(np.flipud(SinTheta[int(par.Nz/2),:,:]), interpolation="None")
#    fig.colorbar(cbr)
#    plt.show() 
#####################



    Diff_Phi_R = np.abs( np.gradient(PotInBox,axis=2) * X/R + np.gradient(PotInBox,axis=1) * Y/R )
    if par.Case == "Mukherjee":
       VelocityPhi      = par.Epsilon * np.sqrt( R * Diff_Phi_R )
    if par.Case == "Standard":
       VelocityPhi      = par.a0 * np.sqrt( R * Diff_Phi_R )
       VelocityPhi_ExpZ = VelocityPhi * np.exp(-np.abs(Z)/par.z0) 

    CosTheta = X/R
    SinTheta = Y/R

    # `PotInBox` have been unnormalized by `Const_C**2`
    if par.Case == "Mukherjee":
       ISM[1] = VelocityPhi * SinTheta / par.Const_C
       ISM[2] = VelocityPhi * CosTheta / par.Const_C
    if par.Case == "Standard":
       ISM[1] = VelocityPhi_ExpZ * SinTheta / par.Const_C
       ISM[2] = VelocityPhi_ExpZ * CosTheta / par.Const_C

    ISM[3] = 0
    ISM[4] = PresInBox

    # Fractalize mass density
    ISM[0] *= FractalDensity
    ISM[1] *= FractalUxyz 
    ISM[2] *= FractalUxyz
    ISM[3] *= FractalUxyz

    return ISM
