import numpy as np
import matplotlib.pyplot as plt 
from scipy.integrate import odeint

NEWTON_G = 1.0

 
def f(y, r, params):
    Psi, Phi = y      # unpack current values of y
    CsSqr, PhiGasOrigin, Rho0_g, Rho0_DM, ScaleRadius = params  # unpack parameters
    derivs = [ - 2*Psi/r + 4*np.pi*NEWTON_G*( GasDensityProfile( Phi, PhiGasOrigin, Rho0_g, CsSqr ) 
               + NFWDensityProfile( r, ScaleRadius, Rho0_DM ) ) , Psi ]
    return derivs




def NFWDensityProfile( r, ScaleRadius, Rho0_DM ):

    a = r / ScaleRadius   

    NFWDensityProfile1D = Rho0_DM / ( a * np.square( 1 + a ) )

    return NFWDensityProfile1D


def GasDensityProfile( Phi, PhiGasOrigin, Rho0_g, CsSqr ):

    GasDensityProfile = Rho0_g * np.exp( - ( Phi - PhiGasOrigin ) / CsSqr )

    return GasDensityProfile


    # Exat NFW potential
def ExactNFWPotential(Rho0_DM, ScaleRadius, ScaleRadius):
    Exact = -4*np.pi*NEWTON_G*Rho0_DM*np.power(ScaleRadius,3) * np.log(1+r/ScaleRadius) / r
    return Exact

def GetTotalDensityProfile( Radius ):

    dr = 2e-4
    r = np.arange(1e-4, Radius, dr)


    # Parameters for DM
    Rho0_DM      = 0.0 # 
    ScaleRadius  = 1.0 # scale radius


    # The total potential value at the center of sphere
    Phi0         = 1e-6
    Psi0         = 0

    # Parameters for gas
    CsSqr        = 1e-5 # The square of sound speed
    PhiGasOrigin = Phi0 # Assuming gas potential and DM potential have the same value at the center of sphere
    Rho0_g       = 1e-2 # Peak density

    # Bundle Parameters
    params = [ CsSqr, PhiGasOrigin, Rho0_g, Rho0_DM, ScaleRadius ]
    
    # Bundle boundary values
    y0     = [ Psi0, Phi0 ]

    TotalPotential = odeint( f, y0, r, args=(params,) )

    fig, ax = plt.subplots( 1, 1, sharex=False, sharey=False )
    fig.subplots_adjust( hspace=0.1, wspace=0.1 )                                                             
    fig.set_size_inches( 10, 5 )



    # Plot total (gas+DM) potential profile
    ax.plot(r, np.absolute(TotalPotential[:,1]), 'x')

    # Plot numerical NFW denisty profile
    ax.plot(r, NFWDensityProfile( r, ScaleRadius, Rho0_DM ), '^')

    # Plot numerical gas density profile
    ax.plot(r, GasDensityProfile( TotalPotential[:,1], [PhiGasOrigin]*TotalPotential[:,1].shape[0],  Rho0_g, CsSqr ), '*' )

    # Plot exact NFW potential profile
    ExactNFWPotential = ExactNFWPotential(Rho0_DM, ScaleRadius, ScaleRadius)
    ax.plot(r, np.absolute(Exact), '>')
     


    ax.set_xscale('log')
    ax.set_yscale('log')

    plt.show()

    return GetTotalDensityProfile


Radius = 1e1

GetTotalDensityProfile(Radius)
