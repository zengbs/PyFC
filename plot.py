import numpy as np
import matplotlib.pyplot as plt 
from scipy.integrate import odeint
from density_profile import *
import parameters as par


def Plot( ):
    par.Parameters()

    Psi0           = par.Phi0    / par.Sigma_D**2
    DevPsi0        = par.DevPhi0 / par.Sigma_D**2

    fig, ax = plt.subplots( 1, 2, sharex=False, sharey=False )
    fig.subplots_adjust( hspace=0.1, wspace=0.1 )                                                             
    fig.set_size_inches( 15, 5 )
    handles0=[]
    handles1=[]

    rPrime = np.arange(par.BPoint, par.Radius/par.Radius_D, par.CoarseDr)

 
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

    TotalPotential        = NumericalTotalPotential( rPrime, Psi0, DevPsi0 )
    GasDensity, DMDensity = NumericalDensity( rPrime )
    E,=ax[0].plot( rPrime*par.Radius_D, GasDensity,    '>', label='Gas density with DM' )
    F,=ax[0].plot( rPrime*par.Radius_D, DMDensity,     '<', label='DM density with Gas' )
    G,=ax[1].plot( rPrime*par.Radius_D, TotalPotential*(par.Sigma_D/par.C)**2,'-', label='Gas+DM potential' )
    


    ax[0].set_xscale('log')
    ax[0].set_yscale('log')

    ax[1].set_xscale('log')
    ax[1].set_yscale('linear')

    ax[0].set_xlabel(r'$r$'   , size=20)
    ax[1].set_xlabel(r'$r$'   , size=20)
    ax[0].set_ylabel(r'$\rho$', size=20)
    ax[1].set_ylabel(r'$\Phi$', size=20)

    #ax[0].legend(handles=[C,A,F,E],loc='lower left', fontsize=12)
    #ax[1].legend(handles=[D,B,G],loc='upper left', fontsize=12)
    ax[0].legend(handles=[F,E],loc='lower left', fontsize=12)
    ax[1].legend(handles=[G],loc='upper left', fontsize=12)
    plt.show()


Plot( )
