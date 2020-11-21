import numpy as np
import matplotlib.pyplot as plt 
from scipy.integrate import odeint

NEWTON_G = 1.0

 
def f(y, r, params):
    Psi, Phi = y      # unpack current values of y
    CsSqr, PhiOrigin, Rho0_g, Rho0_DM, ScaleRadius = params  # unpack parameters
    derivs = [ -2*Psi/r + 4*np.pi*NEWTON_G*( GasProfile( Phi, CsSqr ) + NFWProfile( r, ScaleRadius, Rho0_DM ) ) , Psi ]
    return derivs


NPoints = 5000


def NFWProfile( r, ScaleRadius, Rho0_DM ):

    a = r / ScaleRadius   

    NFWDensityProfile1D = Rho0_DM / ( a * np.square( 1 + a ) )

    return NFWDensityProfile1D


def GasProfile( Potential, CsSqr ):
    Rho0_g = 1.0

    GasProfile = Rho0_g * np.exp( -1.0/CsSqr * ( Potential - Potential[0] ) )

    return GasProfile


def DensityProfile( L ):

    Radius = 0.5*np.sqrt(L[0]**2 + L[1]**2 + L[3]**2)

    r = np.linspace(0, Radius, num=NPoints)

    Profile = np.linspace(0, Radius, num=NPoints)

    # Parameters
    CsSqr       =
    PhiOrigin   =
    Rho0_g      =
    Rho0_DM     = 
    ScaleRadius = 

    # Boundary values at the center of sphere
    Phi0        = 
    Psi0        =

    # Bundle Parameters
    params = [ CsSqr, PhiOrigin, Rho0_g, Rho0_DM, ScaleRadius ]
    
    # Bundle boundary values
    y0     = [ Psi0, Phi0 ]

    DensityProfile = odeint( f, y0, args=(params,) )

    return DensityProfile


