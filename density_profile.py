import sys
import numpy as np
from scipy.integrate import odeint
import parameters as par
import matplotlib.pyplot as plt


np.set_printoptions(threshold=sys.maxsize)
 
def f(y, rPrime, params):
    Psi, DevPsi = y      # unpack current values of y
    Kappa, Lambda, Constant, Psi0, DevPsi0 = params  # unpack parameters

    derivs = [ DevPsi, -2*DevPsi/rPrime + Constant*( np.exp(-Psi) + Lambda*Lambda/Kappa/Kappa*np.exp(-Kappa*Kappa*Psi) ) ]
    return derivs

def FreePara2DerivedPara( ):
    # Peak density of DM     
    Rho0_D = par.Rho0_g * par.Kappa * par.Kappa / par.Lambda / par.Lambda

    # Core radius of DM
    CoreRadius_D = par.Lambda*par.CoreRadius_g

    # Velocity dispersion of DM
    Sigma_D = par.Kappa * par.Sigma_g
    return Rho0_D, Sigma_D, CoreRadius_D


## Exat NFW density
#def ExactNFWDensity( r, ScaleRadius, Rho0_DM ):
#    a = r / ScaleRadius   
#    NFWDensityProfile1D = Rho0_DM / ( a * np.square( 1 + a ) )
#    return NFWDensityProfile1D
#
## Exat NFW potential
#def ExactNFWPotential(r, Rho0_DM, ScaleRadius):
#    Exact = -4*np.pi*NEWTON_G*Rho0_DM*np.power(ScaleRadius,3) * np.log(1+r/ScaleRadius) / r - -4*np.pi*NEWTON_G*Rho0_DM*ScaleRadius**2
#    return Exact

# Density of isothermal gas sphere as a function of total potential
def IsothermalGasDensity( Phi ):
    SoubdSpeedSqr  = par.kB*par.Temp_g/(par.Mu*par.Matom*par.Const_Erg2eV)
    SoubdSpeedSqr /= 1e10 # (km/cm)**2
    GasDensityProfile = par.Rho0_g * np.exp(  -( Phi - par.Phi0 ) / SoubdSpeedSqr )
    return GasDensityProfile

# Density of isothermal DM sphere as a function of total potential
def IsothermalDMDensity( Phi ):
    DMDensityProfile = par.Rho0_D * np.exp(  -( Phi - par.Phi0 ) / par.Sigma_D**2 )
    return DMDensityProfile

def NumericalDensity( rPrime ):
    Psi0           = par.Phi0    / par.Sigma_D**2
    DevPsi0        = par.DevPhi0 / par.Sigma_D**2
   
    Psi            = NumericalTotalPotential( rPrime, Psi0, DevPsi0 )
    Phi            = Psi         * par.Sigma_D**2

    GasDensity     = IsothermalGasDensity( Phi )
    DMDensity      = IsothermalDMDensity ( Phi )
    return GasDensity, DMDensity

"""
return: gravitational potential normalized by par.Sigma_D**2 but unnormalized by `par.C**2`
"""
def NumericalTotalPotential( rPrime, Psi0, DevPsi0 ):
    # Bundle Parameters
    params = [ par.Kappa, par.Lambda, par.Constant, Psi0, DevPsi0 ]


    # Bundle boundary values
    y0     = [ Psi0, DevPsi0 ]
    # Solve ODE       
    TotalPotential = odeint( f, y0, rPrime, args=(params,) )
    return TotalPotential[:,0]

"""
input: PotInBox: gravitational potential unnormalized by `par.Sigma_D**2` but normalized by `par.C**2`
"""
def NumericalISM( PotInBox, FluidInBox, PresInBox, delta, Center ):

    # create an array stored ISM
    ISM       = np.zeros((5, par.Nz, par.Ny, par.Nx), dtype=par.Precision)

    # unnormalized by `par.C**2`
    PotInBox *= par.C**2

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

    X   = (Idx+0.5)*delta[2]-Center[2]
    Y   = (Jdx+0.5)*delta[1]-Center[1]
    Z   = (Kdx+0.5)*delta[0]-Center[0] 
    R   = np.sqrt(X**2+Y**2)


    # expression below have assumed that the potential at the center of sphere is zero
    ISM[0] = par.ISM0 *  np.exp( -( PotInBox -  PotInBoxExtendZ*par.Epsilon**2 )/par.Sigma_t**2 )
  
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



    # In cylindrical coordinate:
    #  ∂Φ     ∂Φ   ∂R     ∂Φ   ∂R     ∂Φ   x     ∂Φ   y
    # ---- = ---- ---- + ---- ---- = ---- --- + ---- ---, where R = ( x**2 + y**2 )**0.5, and Φ=Φ(x, y, z)
    #  ∂R     ∂x   ∂x     ∂y   ∂y     ∂x   R     ∂y   R

    Diff_Phi_R       = np.abs( np.gradient(PotInBox,axis=2) * X/R + np.gradient(PotInBox,axis=1) * Y/R )
    VelocityPhi      = par.Epsilon * np.sqrt( R * Diff_Phi_R )

    # Vx = Vr * cosθ
    # Vy = Vr * sinθ

    # Vx = Vφ * sinθ
    # Vy = Vφ * cosθ
    CosTheta = X/R
    SinTheta = Y/R

    ISM[1] = VelocityPhi * SinTheta / par.C
    ISM[2] = VelocityPhi * CosTheta / par.C
    ISM[3] = 0
    ISM[4] = PresInBox

    return ISM
