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
    GasDensityProfile = Rho0_g * np.exp(  -( Phi - PhiGasOrigin ) / CsSqr )
    return GasDensityProfile


    # Exat NFW potential
def GetExactNFWPotential(r, Rho0_DM, ScaleRadius):
    Exact = -4*np.pi*NEWTON_G*Rho0_DM*np.power(ScaleRadius,3) * np.log(1+r/ScaleRadius) / r
    return Exact

def Plot( Radius ):
    # Parameters for DM
    Rho0_DM      = 0   # 
    ScaleRadius  = 1.0   # scale radius
    
    # Parameters for gas
    Rho0_g       = 1     # Peak density
    CsSqr        = 1     # The square of sound speed



    fig, ax = plt.subplots( 1, 1, sharex=False, sharey=False )
    fig.subplots_adjust( hspace=0.1, wspace=0.1 )                                                             
    fig.set_size_inches( 10, 5 )
    handles=[]

    dr = 1e-3  # cell spacing
    r = np.arange(1e-2, Radius, dr)


    # Plot exact NFW potential profile in the absence of gases
    if ( Rho0_DM > 0 ):
         ExactNFWPotential = GetExactNFWPotential(r, Rho0_DM, ScaleRadius)
         A,=ax.plot(r, np.absolute(ExactNFWPotential), ls='-', color='k', label='Exact NFW potential without gases', zorder=10)
         handles.append(A)
         Phi0         = ExactNFWPotential[0] # The total potential value at the center of sphere
    else:
         Phi0         = 1e-3
 

    Psi0         = 0     # The derivative of potential at the center of sphere

    PhiGasOrigin = Phi0  # Assuming gas potential and DM potential have the same value at the center of sphere

    # Bundle Parameters
    params = [ CsSqr, PhiGasOrigin, Rho0_g, Rho0_DM, ScaleRadius ]
    
    # Bundle boundary values
    y0     = [ Psi0, Phi0 ]

    # Solve ODE
    TotalPotential = odeint( f, y0, r, args=(params,) )


    # Plot numerical NFW denisty profile in the absence of gases
    if ( Rho0_DM > 0 ):
         B,=ax.plot(r, NFWDensityProfile( r, ScaleRadius, Rho0_DM ), ls='-', color='g', label='Exact NFW denisty profile without gases', zorder=9)
         handles.append(B)

    # Plot total (gas+DM) potential profile
    if ( Rho0_g == 0 ):
         label='Numerical NFW potential solving from the Poisson solver'
    elif ( Rho0_DM == 0 ):
         label='Numerical gas potential solving from the Poisson solver'
    else:
         label='Numerical total potential solving from the Poisson solver'

    C,=ax.plot(r, np.absolute(TotalPotential[:,1]), 'x', label=label, zorder=8, color='r')
    handles.append(C)

    # Plot numerical gas density profile in the presence of DM
    if ( Rho0_g > 0 and Rho0_DM > 0 ):
         GasDensityProfile = GetGasDensityProfile( TotalPotential[:,1], [PhiGasOrigin]*TotalPotential[:,1].shape[0], Rho0_g, CsSqr )
         D,=ax.plot(r, GasDensityProfile, '*', label='Numerical gas denisty profile with DM', zorder=7 )
         handles.append(D)



    ## Plot numerical gas density profile in the absence of DM ##
    # Bundle Parameters
    Rho0_DM = 0
    params  = [ CsSqr, PhiGasOrigin, Rho0_g, Rho0_DM, ScaleRadius ]
    
    # Bundle boundary values
    y0             = [ Psi0, Phi0 ]
    TotalPotential = odeint( f, y0, r, args=(params,) )

    if ( Rho0_g > 0 ):
         GasDensityProfile = GetGasDensityProfile( TotalPotential[:,1], [PhiGasOrigin]*TotalPotential[:,1].shape[0], Rho0_g, CsSqr )
         E,=ax.plot(r, GasDensityProfile, 'o', label='Numerical gas denisty profile without DM', zorder=6 )
         handles.append(E)


    ax.set_xscale('log')
    ax.set_yscale('log')

    ax.legend(handles=handles,loc='lower left', fontsize=16)
    plt.show()



Radius = 1e2

Plot(Radius)
