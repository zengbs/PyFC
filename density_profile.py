import numpy as np
from scipy.integrate import odeint

NEWTON_G = 1.0

 
def f(y, r, params):
    Psi, Phi = y      # unpack current values of y
    CsSqr, PhiGasOrigin, Rho0_g, Rho0_DM, ScaleRadius = params  # unpack parameters
    # Assuming the self-gravity of gases is minor compared to that of DM
    #derivs = [ - 2*Psi/r + 4*np.pi*NEWTON_G*ExactNFWDensity( r, ScaleRadius, Rho0_DM ) , Psi ]

    # Consider the self-gravity of gases as well
    derivs = [ - 2*Psi/r + 4*np.pi*NEWTON_G*( IsothermalDensity( Phi, PhiGasOrigin, Rho0_g, CsSqr ) 
               + ExactNFWDensity( r, ScaleRadius, Rho0_DM ) ) , Psi ]
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

# Density of isothermal sphere as a function of total potential
def IsothermalDensity( Phi, PhiGasOrigin, Rho0_g, CsSqr ):
    GasDensityProfile = Rho0_g * np.exp(  -( Phi - PhiGasOrigin ) / CsSqr )
    return GasDensityProfile

def NumericalGasDensity( r, Rho0_DM, ScaleRadius, Rho0_g, CsSqr ):
    PhiGasOrigin   = ExactNFWPotential(r[0], Rho0_DM, ScaleRadius)
    TotalPotential = NumericalTotalPotential( r, Rho0_DM, ScaleRadius, Rho0_g, CsSqr )
    GasDensity     = IsothermalDensity( TotalPotential, [PhiGasOrigin]*TotalPotential.shape[0], Rho0_g, CsSqr )
    return GasDensity
  
def NumericalTotalPotential( r, Rho0_DM, ScaleRadius, Rho0_g, CsSqr ):
    # Psi0: The derivative of potential at the center of sphere
    # Phi0: The               potential at the center of sphere
    if ( Rho0_g == 0 ):
         Psi0 = 0.5*ScaleRadius**-2
         Phi0 = 0
    elif( Rho0_DM == 0 ):
         Psi0 = 0
         Phi0 = 0
    else:
         Psi0 = 0
         Phi0 = 0

    # Assuming gas potential and total potential have the same value at the center of sphere
    PhiGasOrigin = Phi0

    # Bundle Parameters
    params = [ CsSqr, PhiGasOrigin, Rho0_g, Rho0_DM, ScaleRadius ]
                      
    # Bundle boundary values
    y0     = [ Psi0, Phi0 ]
                      
    # Solve ODE       
    TotalPotential = odeint( f, y0, r, args=(params,) )
    return TotalPotential[:,1]
