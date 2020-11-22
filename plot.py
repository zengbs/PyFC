import numpy as np
import matplotlib.pyplot as plt 
from scipy.integrate import odeint
from DensityProfile import *


def Plot( Radius ):
    # Parameters for DM
    ScaleRadius  = 0.5   # scale radius
    
    # Parameters for gas
    CsSqr        = 0.2     # The square of sound speed



    fig, ax = plt.subplots( 1, 2, sharex=False, sharey=False )
    fig.subplots_adjust( hspace=0.1, wspace=0.1 )                                                             
    fig.set_size_inches( 15, 5 )
    #handles0=[]
    #handles1=[]

    dr = 1e-4  # cell spacing
    r = np.arange(1e-4, Radius, dr)

    # Gas-only
    Rho0_DM = 0
    Rho0_g  = 1
    GasOnlyPotential = NumericalTotalPotential( r, Rho0_DM, ScaleRadius, Rho0_g, CsSqr )
    GasOnlyDensity   = NumericalGasDensity( r, Rho0_DM, ScaleRadius, Rho0_g, CsSqr )
    ax[0].plot( r, GasOnlyDensity,   '*' )
    ax[1].plot( r, GasOnlyPotential, '^' )

    # DM-only
    Rho0_DM = 1
    Rho0_g  = 0
    DMOnlyPotential  = NumericalTotalPotential( r, Rho0_DM, ScaleRadius, Rho0_g, CsSqr )
    ExactDensity     = ExactNFWDensity( r, ScaleRadius, Rho0_DM )
    ExactPotential   = ExactNFWPotential(r, Rho0_DM, ScaleRadius)
    ax[0].plot( r, ExactDensity,     ls='-' ,color='k', zorder=9)
    ax[1].plot( r, ExactPotential,   ls='-' ,color='k', zorder=10)
    ax[1].plot( r, DMOnlyPotential,  'o'    )
    
    # DM and gas
    Rho0_DM = 1
    Rho0_g  = 1
    DMGasPotential   = NumericalTotalPotential( r, Rho0_DM, ScaleRadius, Rho0_g, CsSqr )
    DMGasDensity     = NumericalGasDensity( r, Rho0_DM, ScaleRadius, Rho0_g, CsSqr )
    ax[0].plot( r, DMGasDensity,   '*' )
    ax[1].plot( r, DMGasPotential, '^' )
    



    ax[0].set_xscale('log')
    ax[0].set_yscale('log')

    ax[1].set_xscale('log')
    ax[1].set_yscale('linear')

    #ax[0].legend(handles=handles0,loc='lower left', fontsize=16)
    #ax[1].legend(handles=handles1,loc='upper left', fontsize=16)
    plt.show()



Radius = 1e1

Plot(Radius)

