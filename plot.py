import numpy as np
import matplotlib.pyplot as plt 
from scipy.integrate import odeint
from density_profile import *


def Plot( ParaPhy, ParaNum, Radius ):
    # Unbundle physical parameters
    Radius_g, Rho0_g, Sigma_g, Lambda, Kappa, Temp_g, Constant, Phi0, DevPhi0 = ParaPhy

    # Unbundle numerical parameters
    BPoint, CoarseDr, Radius = ParaNum   

    # Get derived parameters
    Rho0_D, Sigma_D, Radius_D  = Free2DerivedPara( Rho0_g, Sigma_g, Radius_g, Lambda, Kappa )

    Psi0           = Phi0    / Sigma_D / Sigma_D
    DevPsi0        = DevPhi0 / Sigma_D / Sigma_D

    fig, ax = plt.subplots( 1, 2, sharex=False, sharey=False )
    fig.subplots_adjust( hspace=0.1, wspace=0.1 )                                                             
    fig.set_size_inches( 15, 5 )
    handles0=[]
    handles1=[]

    rPrime = np.arange(BPoint, Radius, CoarseDr)

 
    ## Gas-only
    #ParaPhy[4]              = 0 # Kappa = 0
    #GasOnlyPotential        = NumericalTotalPotential( rPrime, Kappa, Lambda, Constant, Psi0, DevPsi0 )
    #GasOnlyDensity, Nothing = NumericalDensity( rPrime, ParaPhy )
    #A,=ax[0].plot( rPrime, GasOnlyDensity  , '*', label='Gas density without DM'   )
    #B,=ax[1].plot( rPrime, GasOnlyPotential, '^', label='Gas potential without DM' )

    ## DM-only
    #ParaPhy[1]             = 0 # Rho0_g = 0
    #DMOnlyPotential        = NumericalTotalPotential( rPrime, Kappa, Lambda, Constant, Psi0, DevPsi0 )
    #Nothing, DMOnlyDensity = NumericalDensity( rPrime, ParaPhy )
    #C,=ax[0].plot( rPrime, DMOnlyDensity  , 'o', label='DM density without Gas'   )
    #D,=ax[1].plot( rPrime, DMOnlyPotential, 'x', label='DM potential without Gas' )
    
    # DM and gas
    #ParaPhy[1]            = Rho0_g
    #ParaPhy[4]            = Kappa
    TotalPotential        = NumericalTotalPotential( rPrime, Kappa, Lambda, Constant, Psi0, DevPsi0 )
    GasDensity, DMDensity = NumericalDensity( rPrime, ParaPhy )
    E,=ax[0].plot( rPrime, GasDensity,    '>', label='Gas density with DM' )
    F,=ax[0].plot( rPrime, DMDensity,     '<', label='DM density with Gas' )
    G,=ax[1].plot( rPrime, TotalPotential,'-', label='Gas+DM potential' )
    


    ax[0].set_xscale('log')
    ax[0].set_yscale('log')

    ax[1].set_xscale('log')
    ax[1].set_yscale('linear')

    ax[0].set_xlabel(r'$r/r_{D}$'          , size=20)
    ax[1].set_xlabel(r'$r/r_{D}$'          , size=20)
    ax[0].set_ylabel(r'$\rho$'             , size=20)
    ax[1].set_ylabel(r'$\Phi/\sigma_{D}^2$', size=20)

    #ax[0].legend(handles=[C,A,F,E],loc='lower left', fontsize=12)
    #ax[1].legend(handles=[D,B,G],loc='upper left', fontsize=12)
    ax[0].legend(handles=[F,E],loc='lower left', fontsize=12)
    ax[1].legend(handles=[G],loc='upper left', fontsize=12)
    plt.show()

##########################
###  Free parameters   ###
##########################
# The eq(2) in Sutherland & Bicknell (2007)
Constant = 9

# Core radius of gas sphere (kpc)
Radius_g = 1

# Peak gas density (1/cm^3)
Rho0_g = 0.5

# Velocity dispersion of gas (km/s)
Sigma_g  = 250

# Lambda = r_{D}/r_{g}
Lambda = 5

# Kappa = \sigma_{D}/\sigma_{g}
Kappa = 2

# Temperature of gas (K)
Temp_g = 1e7

# The potential at the center of sphere
Phi0 = 0

# The spatial derivative of potential at the center of sphere
DevPhi0 = 0

##########################
### Numerical parameters ###
##########################

# The boundary point for integration
BPoint   = 1e-4

# The integration step
CoarseDr = 1e-4

##########################
### Derived parameters ###
##########################

Rho0_D, Sigma_D, Radius_D = Free2DerivedPara( Rho0_g, Sigma_g, Radius_g, Lambda, Kappa )


############################
###   Bundle paramters   ###
############################
# plot radius in unit of Radius_D (core radius of DM)
Radius = 6.9

ParaPhy = [ Radius_g, Rho0_g, Sigma_g, Lambda, Kappa, Temp_g, Constant, Phi0, DevPhi0 ]
ParaNum = [ BPoint, CoarseDr, Radius ]

Plot( ParaPhy, ParaNum, Radius )
