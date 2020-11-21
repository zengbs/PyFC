import numpy as np
import matplotlib.pyplot as plt 
from scipy.integrate import odeint

NEWTON_G = 1.0

 
def f(y, r, params):
    Psi, Phi = y      # unpack current values of y
    CsSqr, PhiGasOrigin, Rho0_g, Rho0_DM, ScaleRadius = params  # unpack parameters
    derivs = [ - 2*Psi/r + 4*np.pi*NEWTON_G*( GetGasDensityProfile( Phi, PhiGasOrigin, Rho0_g, CsSqr ) 
               + NFWDensityProfile( r, ScaleRadius, Rho0_DM ) ) , Psi ]
    return derivs




def NFWDensityProfile( r, ScaleRadius, Rho0_DM ):

    a = r / ScaleRadius   

    NFWDensityProfile1D = Rho0_DM / ( a * np.square( 1 + a ) )

    return NFWDensityProfile1D


def GetGasDensityProfile( Phi, PhiGasOrigin, Rho0_g, CsSqr ):

    GasDensityProfile = Rho0_g * np.exp( - ( Phi - PhiGasOrigin ) / CsSqr )

    return GasDensityProfile


    # Exat NFW potential
def GetExactNFWPotential(r, Rho0_DM, ScaleRadius):
    Exact = -4*np.pi*NEWTON_G*Rho0_DM*np.power(ScaleRadius,3) * np.log(1+r/ScaleRadius) / r
    return Exact

def GetTotalDensityProfile( Radius ):
    fig, ax = plt.subplots( 1, 1, sharex=False, sharey=False )
    fig.subplots_adjust( hspace=0.1, wspace=0.1 )                                                             
    fig.set_size_inches( 10, 5 )

    dr = 2e-4
    r = np.arange(1e-4, Radius, dr)


    # Parameters for DM
    Rho0_DM      = 1   # 
    ScaleRadius  = 1.0   # scale radius

    # Plot exact NFW potential profile
    if ( Rho0_DM > 0 ):
         ExactNFWPotential = GetExactNFWPotential(r, Rho0_DM, ScaleRadius)
         ax.plot(r, np.absolute(ExactNFWPotential), '>')

    # The total potential value at the center of sphere
    if ( Rho0_DM > 0 ):
         Phi0         = ExactNFWPotential[0]
    else:
         Phi0         = 1e-6

    Psi0         = 0     # The derivative of potential at the center of sphere

    # Parameters for gas
    Rho0_g       = 1     # Peak density
    CsSqr        = 1e-2  # The square of sound speed
    PhiGasOrigin = Phi0  # Assuming gas potential and DM potential have the same value at the center of sphere

    # Bundle Parameters
    params = [ CsSqr, PhiGasOrigin, Rho0_g, Rho0_DM, ScaleRadius ]
    
    # Bundle boundary values
    y0     = [ Psi0, Phi0 ]

    TotalPotential = odeint( f, y0, r, args=(params,) )




    # Plot total (gas+DM) potential profile
    ax.plot(r, np.absolute(TotalPotential[:,1]), 'x')

    # Plot numerical NFW denisty profile
    if ( Rho0_DM > 0 ):
         ax.plot(r, NFWDensityProfile( r, ScaleRadius, Rho0_DM ), '^')

    # Plot numerical gas density profile
    if ( Rho0_g > 0 ):
         GasDensityProfile = GetGasDensityProfile( TotalPotential[:,1], [PhiGasOrigin]*TotalPotential[:,1].shape[0], Rho0_g, CsSqr )
         ax.plot(r, GasDensityProfile, '*' )


    ax.set_xscale('log')
    ax.set_yscale('log')

    plt.show()

    return GetTotalDensityProfile


Radius = 1e1

GetTotalDensityProfile(Radius)
