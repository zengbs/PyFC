import numpy as np
from pri2con import Pri2Con
from sphere import SphericalSphere

def _DumpFile(out):


    # Output file name
    filename = "UM_IC"

    # Number of cells along x/y/z
    Nx = out.shape[0] 
    Ny = out.shape[1] 
    Nz = out.shape[2] 

    # Box size
    Lx = 100
    Ly = 100
    Lz = 100

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
  
    ##########################
    ### Derived parameters ###
    ##########################

    # Peak density of DM     
    Rho0_D = Rho0_g * Kappa * Kappa / Lambda / Lambda

    # Core radius of DM
    Radius_D = Lambda*Radius_g

    # Velocity dispersion of DM
    Sigma_D = Kappa * Sigma_g


    ############################
    ###   Bundle paramters   ###
    ############################
    Para = [ Radius_g, Rho0_g, Sigma_g, Temp_g, Radius_D, Rho0_D, Sigma_D, Lambda, Kappa, Constant ]


    # Number of cells along each dimension of the input grid.
    N = np.array([Nx, Ny, Nz], dtype='int')
    L = np.array([Lx, Ly, Lz], dtype='float64') 

    FluidInBox = SphericalSphere( L, N, Para )

    FluidInBox.tofile(filename)

out=np.zeros((512,512,512))
_DumpFile(out)
