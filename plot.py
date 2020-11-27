import numpy as np
import matplotlib.pyplot as plt 
from scipy.integrate import odeint
from density_profile import *


def Plot( ParaPhy, ParaNum, PlotRadius ):
    # Unbundle physical parameters
    Radius_g, Rho0_g, Sigma_g, Lambda, Kappa, Temp_g, Constant, Phi0, DevPhi0 = ParaPhy

    # Unbundle numerical parameters
    BPoint, CoarseDr = ParaNum   

    # Get derived parameters
    Rho0_D, Sigma_D, Radius_D  = Free2DerivedPara( Rho0_g, Sigma_g, Radius_g, Lambda, Kappa )

    Psi0           = Phi0    / Sigma_D / Sigma_D
    DevPsi0        = DevPhi0 / Sigma_D / Sigma_D

    fig, ax = plt.subplots( 1, 2, sharex=False, sharey=False )
    fig.subplots_adjust( hspace=0.1, wspace=0.1 )                                                             
    fig.set_size_inches( 15, 5 )
    handles0=[]
    handles1=[]

    rPrime = np.arange(BPoint, PlotRadius, CoarseDr)

    r = rPrime*Radius_D
 
    # Gas-only
    Rho0_D_Temp = 0
    GasOnlyPotential        = NumericalTotalPotential( rPrime, Kappa, Lambda, Constant, Psi0, DevPsi0 )
    GasOnlyDensity, Nothing = NumericalDensity( rPrime, ParaPhy )
    GasOnlyPotential       *= Sigma_D*Sigma_D
    A,=ax[0].plot( r, GasOnlyDensity  , '*', label='Gas density without DM'   )
    B,=ax[1].plot( r, GasOnlyPotential, '^', label='Gas potential without DM' )

    # DM-only
    Rho0_g_Temp  = 0
    DMOnlyPotential       = NumericalTotalPotential( rPrime, Kappa, Lambda, Constant, Psi0, DevPsi0 )
    Nothing, ExactDensity = ExactNFWDensity( rPrime, ParaPhy )
    DMOnlyPotential       *= Sigma_D*Sigma_D
    C,=ax[0].plot( r, DMOnlyPotential, 'o', label='DM density without Gas'   )
    D,=ax[1].plot( r, DMOnlyPotential, 'x', label='DM potential without Gas' )
    
    # DM and gas
    TotalPotential        = NumericalTotalPotential( rPrime, Kappa, Lambda, Constant, Psi0, DevPsi0 )
    GasDensity, DMDensity = NumericalGasDensity( rPrime, ParaPhy )
    TotalPotential       *= Sigma_D*Sigma_D
    E,=ax[0].plot( r, GasDensity,    '>', label='Gas density with DM' )
    F,=ax[0].plot( r, DMDensity,     '<', label='DM density with Gas' )
    G,=ax[1].plot( r, TotalPotential,'#', label='Gas+DM potential' )
    


    ax[0].set_xscale('log')
    ax[0].set_yscale('log')

    ax[1].set_xscale('log')
    ax[1].set_yscale('linear')


    ax[0].legend(handles=[C,A,F],loc='lower left', fontsize=12)
    ax[1].legend(handles=[D,E,B,G],loc='upper left', fontsize=12)
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
Sigma_g = 250  

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
ParaPhy = [ Radius_g, Rho0_g, Sigma_g, Lambda, Kappa, Temp_g, Constant, Phi0, DevPhi0 ]
ParaNum = [ BPoint, CoarseDr ]

PlotRadius = 10

Plot( ParaPhy, ParaNum, PlotRadius )
