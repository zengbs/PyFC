import numpy as np
import matplotlib.pyplot as plt 
from scipy.integrate import odeint
from density_profile import *


def Plot( Radius, dr0, dr, Rho0_g, CsSqr, ScaleRadius, Rho0_DM ):

    fig, ax = plt.subplots( 1, 2, sharex=False, sharey=False )
    fig.subplots_adjust( hspace=0.1, wspace=0.1 )                                                             
    fig.set_size_inches( 15, 5 )
    handles0=[]
    handles1=[]

    r = np.arange(dr0, Radius, dr)

    # Gas-only
    Rho0_DM_Temp = 0
    GasOnlyPotential = NumericalTotalPotential( r, Rho0_DM_Temp, ScaleRadius, Rho0_g, CsSqr )
    GasOnlyDensity   = NumericalGasDensity( r, Rho0_DM_Temp, ScaleRadius, Rho0_g, CsSqr )
    A,=ax[0].plot( r, GasOnlyDensity  , '*', label='Gas density without DM' )
    B,=ax[1].plot( r, GasOnlyPotential, '^', label='Gas potential without DM' )

    # DM-only
    Rho0_g_Temp  = 0
    DMOnlyPotential  = NumericalTotalPotential( r, Rho0_DM, ScaleRadius, Rho0_g_Temp, CsSqr )
    ExactDensity     = ExactNFWDensity( r, ScaleRadius, Rho0_DM )
    ExactPotential   = ExactNFWPotential(r, Rho0_DM, ScaleRadius)
    C,=ax[0].plot( r, ExactDensity   ,  ls='-' ,color='k', zorder=9 , label='Exact DM density without Gas' )
    D,=ax[1].plot( r, ExactPotential ,  ls='-' ,color='k', zorder=10, label='Exact DM potential without Gas' )
    E,=ax[1].plot( r, DMOnlyPotential,     'o'                      , label='DM potential without Gas solving from the Poisson solver' )
    
    # DM and gas
    TotalPotential   = NumericalTotalPotential( r, Rho0_DM, ScaleRadius, Rho0_g, CsSqr )
    GasDensity       = NumericalGasDensity( r, Rho0_DM, ScaleRadius, Rho0_g, CsSqr )
    F,=ax[0].plot( r, GasDensity,   '>' , label='Gas density with DM' )
    G,=ax[1].plot( r, TotalPotential, 'x' , label='Gas+DM potential', zorder=1 )
    

    ax[0].set_xscale('log')
    ax[0].set_yscale('log')

    ax[1].set_xscale('log')
    ax[1].set_yscale('linear')


    ax[0].legend(handles=[C,A,F],loc='lower left', fontsize=12)
    ax[1].legend(handles=[D,E,B,G],loc='upper left', fontsize=12)
    plt.show()


# Plot radius
Radius       = 1e2

# Integration step
dr           = 3e-4

# Initial radius for integration
dr0          = 1e-4

# Peak gas density
Rho0_g       = 1

# Sound speed
CsSqr        = 0.02

# DM scale radius
ScaleRadius  = 1

# DM scale density
Rho0_DM      = 1

Plot( Radius, dr0, dr, Rho0_g, CsSqr, ScaleRadius, Rho0_DM )
