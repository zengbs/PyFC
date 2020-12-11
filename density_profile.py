import numpy as np
from scipy.integrate import odeint
import parameters as par


 
def f(y, rPrime, params):
    Psi, DevPsi = y      # unpack current values of y
    Kappa, Lambda, Constant, Psi0, DevPsi0 = params  # unpack parameters

    derivs = [ DevPsi, -2*DevPsi/rPrime + Constant*( np.exp(-Psi) + Lambda*Lambda/Kappa/Kappa*np.exp(-Kappa*Kappa*Psi) ) ]
    return derivs

def FreePara2DerivedPara( Rho0_g, Sigma_g, Radius_g, Lambda, Kappa ):
    # Peak density of DM     
    Rho0_D = Rho0_g * Kappa * Kappa / Lambda / Lambda

    # Core radius of DM
    Radius_D = Lambda*Radius_g

    # Velocity dispersion of DM
    Sigma_D = Kappa * Sigma_g
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
def IsothermalGasDensity( Phi, Phi0, Rho0_g, Sigma_g ):
    GasDensityProfile = Rho0_g * np.exp(  -( Phi - Phi0 ) / Sigma_g / Sigma_g )
    return GasDensityProfile

# Density of isothermal DM sphere as a function of total potential
def IsothermalDMDensity( Phi, Phi0, Rho0_D, Sigma_D ):
    DMDensityProfile = Rho0_D * np.exp(  -( Phi - Phi0 ) / Sigma_D / Sigma_D )
    return DMDensityProfile

def NumericalDensity( rPrime ):
    Psi0           = par.Phi0    / par.Sigma_D / par.Sigma_D
    DevPsi0        = par.DevPhi0 / par.Sigma_D / par.Sigma_D
   
    Psi            = NumericalTotalPotential( rPrime, par.Kappa, par.Lambda, par.Constant, Psi0, DevPsi0 )
    Phi            = Psi * par.Sigma_D * par.Sigma_D

    GasDensity     = IsothermalGasDensity( Phi, [par.Phi0]*Phi.shape[0], par.Rho0_g, par.Sigma_g )
    DMDensity      = IsothermalDMDensity ( Phi, [par.Phi0]*Phi.shape[0], par.Rho0_D, par.Sigma_D )
    return GasDensity, DMDensity
  
def NumericalTotalPotential( rPrime, Kappa, Lambda, Constant, Psi0, DevPsi0 ):
    # Bundle Parameters
    params = [ Kappa, Lambda, Constant, Psi0, DevPsi0 ]

    # Bundle boundary values
    y0     = [ Psi0, DevPsi0 ]
    # Solve ODE       
    TotalPotential = odeint( f, y0, rPrime, args=(params,) )
    return TotalPotential[:,0]


#def NumericalISMDens():
#    Phi = NumericalTotalPotential()
