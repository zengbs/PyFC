import numpy as np
import matplotlib.pyplot as plt 
from scipy.integrate import odeint

NEWTON_G = 1.0

 
def f(y, r, params):
    Psi, Phi = y      # unpack current values of y
    CsSqr, PhiOrigin, Rho0_g, Rho0_DM, ScaleRadius = params  # unpack parameters
    derivs = [ - 2*Psi/r + 4*np.pi*NEWTON_G*( GasProfile( Phi, PhiOrigin, Rho0_g, CsSqr ) 
               + NFWProfile( r, ScaleRadius, Rho0_DM ) ) , Psi ]
    return derivs




def NFWProfile( r, ScaleRadius, Rho0_DM ):

    a = r / ScaleRadius   

    NFWDensityProfile1D = Rho0_DM / ( a * np.square( 1 + a ) )

    return NFWDensityProfile1D


def GasProfile( Phi, PhiOrigin, Rho0_g, CsSqr ):

    GasProfile = Rho0_g * np.exp( -1.0/CsSqr * ( Phi - PhiOrigin ) )

    return GasProfile


def DensityProfile( Radius ):

    dr = 5e-4
    r = np.arange(1e-4, Radius, dr)


    # Parameters for DM
    Rho0_DM     = 1e-3
    ScaleRadius = 1.0

    # Exat NFW potential
    Exact = -4*np.pi*NEWTON_G*Rho0_DM*np.power(ScaleRadius,3) * np.log(1+r/ScaleRadius) / r

    # The potential value at the center of sphere
    Phi0        = Exact[0]
    Psi0        = 0

    # Parameters for gas
    CsSqr       = 0.2
    PhiOrigin   = Exact[0]
    Rho0_g      = 0

    # Bundle Parameters
    params = [ CsSqr, PhiOrigin, Rho0_g, Rho0_DM, ScaleRadius ]
    
    # Bundle boundary values
    y0     = [ Psi0, Phi0 ]

    Potential = odeint( f, y0, r, args=(params,) )

    fig, ax = plt.subplots( 1, 1, sharex=False, sharey=False )
    fig.subplots_adjust( hspace=0.1, wspace=0.1 )                                                             
    fig.set_size_inches( 10, 5 )

    #ax.plot(r, Potential[:,0], 'o')
    ax.plot(r, np.absolute(Potential[:,1]), 'x')
    #ax.plot(r, NFWProfile( r, ScaleRadius, Rho0_DM ), '^')


    ax.plot(r, np.absolute(Exact), '>')
    


    ax.set_xscale('log')
    ax.set_yscale('log')

    plt.show()

    return DensityProfile


Radius = 1e1

DensityProfile(Radius)
