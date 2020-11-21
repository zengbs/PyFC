import numpy as np
import matplotlib.pyplot as plt 
from scipy.integrate import odeint

NEWTON_G = 1.0

 
def f(y, r, params):
    Psi, Phi = y      # unpack current values of y
    CsSqr, PhiOrigin, Rho0_g, Rho0_DM, ScaleRadius = params  # unpack parameters
    derivs = [ -2*Psi/r + 4*np.pi*NEWTON_G*( GasProfile( Phi, PhiOrigin, CsSqr ) + NFWProfile( r, ScaleRadius, Rho0_DM ) ) , Psi ]
    return derivs


NPoints = 5000


def NFWProfile( r, ScaleRadius, Rho0_DM ):

    a = r / ScaleRadius   

    NFWDensityProfile1D = Rho0_DM / ( a * np.square( 1 + a ) )

    return NFWDensityProfile1D


def GasProfile( Phi, PhiOrigin, CsSqr ):
    Rho0_g = 1.0

    GasProfile = Rho0_g * np.exp( -1.0/CsSqr * ( Phi - PhiOrigin ) )

    return GasProfile


def DensityProfile( L ):

    Radius = 0.5*np.sqrt(L[0]**2 + L[1]**2 + L[2]**2)

    #r = np.logspace(-3, np.log(Radius), 500)
    r = np.arange(1e-3, Radius, 5e-3)
    print(r)

    # Parameters
    CsSqr       = 0.2
    PhiOrigin   = 0.0
    Rho0_g      = 0.0
    Rho0_DM     = 1.0
    ScaleRadius = 1.0

    # Boundary values at the center of sphere
    Phi0        = 0
    Psi0        = 0

    # Bundle Parameters
    params = [ CsSqr, PhiOrigin, Rho0_g, Rho0_DM, ScaleRadius ]
    
    # Bundle boundary values
    y0     = [ Psi0, Phi0 ]

    DensityProfile = odeint( f, y0, r, args=(params,) )

    fig, ax = plt.subplots( 1, 1, sharex=False, sharey=False )
    fig.subplots_adjust( hspace=0.1, wspace=0.1 )                                                             
    fig.set_size_inches( 10, 5 )

    ax.plot(r, DensityProfile[:,0], 'o')
    ax.plot(r, DensityProfile[:,1], 'o')


    ax.set_xscale('log')
    ax.set_yscale('log')

    plt.show()

    return DensityProfile


LSize = 1e3

DensityProfile([LSize]*3)
