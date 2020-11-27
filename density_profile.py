import numpy as np
from scipy.integrate import odeint

NEWTON_G = 1.0

 
def f(y, r, params):
    Psi, Omega = y      # unpack current values of y
    rPrime, Kappa, Lambda, Constant = params  # unpack parameters

    derivs = [ Omega, - 2*Omega/rPrime + Constant*( np.exp(-Psi) + Lambda*Lambda/Kappa/Kappa*np.exp(-Kappa*Kappa*Psi) ) ]
    return derivs


# Exat NFW density
def ExactNFWDensity( r, ScaleRadius, Rho0_DM ):
    a = r / ScaleRadius   
    NFWDensityProfile1D = Rho0_DM / ( a * np.square( 1 + a ) )
    return NFWDensityProfile1D

# Exat NFW potential
def ExactNFWPotential(r, Rho0_DM, ScaleRadius):
    Exact = -4*np.pi*NEWTON_G*Rho0_DM*np.power(ScaleRadius,3) * np.log(1+r/ScaleRadius) / r - -4*np.pi*NEWTON_G*Rho0_DM*ScaleRadius**2
    return Exact

# Density of isothermal gas sphere as a function of total potential
def IsothermalGasDensity( Phi, PhiOrigin, Rho0_g, Sigma_g ):
    GasDensityProfile = Rho0_g * np.exp(  -( Phi - PhiOrigin ) / Sigma_g / Sigma_g )
    return GasDensityProfile

# Density of isothermal DM sphere as a function of total potential
def IsothermalDMDensity( Phi, PhiOrigin, Rho0_D, Sigma_D ):
    DMDensityProfile = Rho0_D * np.exp(  -( Phi - PhiOrigin ) / Sigma_D / Sigma_D )
    return DMDensityProfile

def NumericalDensity( rPrime, Kappa, Lambda, Constant, Sigma_g, Rho0_g, PhiOrigin ):
    Psi            = NumericalTotalPotential( rPrime, Kappa, Lambda, Constant, PhiOrigin )
    Rho0_D         = Rho0_g * Kappa * Kappa / Lambda / Lambda
    Sigma_D        = Kappa * Sigma_g

    GasDensity     = IsothermalGasDensity( Psi, [PhiOrigin]*TotalPotential.shape[0], Rho0_g, Sigma_g )
    DMDensity      = IsothermalDMDensity ( Psi, [PhiOrigin]*TotalPotential.shape[0], Rho0_D, Sigma_D )
    return GasDensity, DMDensity
  
def NumericalTotalPotential( rPrime, Kappa, Lambda, Constant, PhiOrigin ):
    # Bundle Parameters
    params = [ rPrime, Kappa, Lambda, Constant ]
                      
    # Bundle boundary values
    y0     = [ Psi0, Phi0 ]
                      
    # Solve ODE       
    TotalPotential = odeint( f, y0, r, args=(params,) )
    return TotalPotential[:,1]
