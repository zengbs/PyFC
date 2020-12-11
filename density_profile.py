import numpy as np
from scipy.integrate import odeint
import parameters as par


 
def f(y, rPrime, params):
    Psi, DevPsi = y      # unpack current values of y
    Kappa, Lambda, Constant, Psi0, DevPsi0 = params  # unpack parameters

    derivs = [ DevPsi, -2*DevPsi/rPrime + Constant*( np.exp(-Psi) + Lambda*Lambda/Kappa/Kappa*np.exp(-Kappa*Kappa*Psi) ) ]
    return derivs

def FreePara2DerivedPara( ):
    # Peak density of DM     
    Rho0_D = par.Rho0_g * par.Kappa * par.Kappa / par.Lambda / par.Lambda

    # Core radius of DM
    Radius_D = par.Lambda*par.Radius_g

    # Velocity dispersion of DM
    Sigma_D = par.Kappa * par.Sigma_g
    return Rho0_D, Sigma_D, Radius_D


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
    GasDensityProfile = par.Rho0_g * np.exp(  -( Phi - par.Phi0 ) / par.Sigma_g**2 )
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
def NumericalISM( PotInBox ):

    # create an array stored ISM
    ISM             = np.zeros((5, par.Nx, par.Ny, par.Nz), dtype=par.Precision)

    # unnormalized by `par.C**2`
    PotInBox       *= par.C**2

    # remove ghost zone inside `PotInBox`
    PotInBox        = PotInBox[:, :, par.GRA_GHOST_SIZE:par.Nz-par.GRA_GHOST_SIZE]
    PotInBox        = PotInBox[:, par.GRA_GHOST_SIZE:par.Ny-par.GRA_GHOST_SIZE, :]
    PotInBox        = PotInBox[par.GRA_GHOST_SIZE:par.Nx-par.GRA_GHOST_SIZE, :, :]

    # extract the potential on the equator
    if par.Nz%2 == 0:
       PotOnEquator = 0.5*( PotInBox[:,:,par.Nz/2-1] + PotInBox[:,:,par.Nz/2] )
    else:
       PotOnEquator = PotInBox[:,:,(par.Nz-1)*0.5]


    PotInBoxExtendZ = np.tile(PotOnEquator,(par.Nz,1,1))


    Kdx = np.indices((par.Nx,par.Ny,par.Nz))[2]
    Z   = (Kdx+0.5)*delta[2]-Center[2]

    a = par.a0 * np.exp(-Z/par.z0)

    # expression below assume the potential at the center of sphere is zero
    ISM = par.Rho0_g * np.exp( -( PotInBox -  PotInBoxExtendZ*a**2 )/par.Sigma_g**2 ) 

    return ISM
