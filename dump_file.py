import numpy as np
from pri2con import Pri2Con
from sphere import SphericalSphere
from density_profile import Free2DerivedPara

def _DumpFile(out):


    # Output file name
    filename = "UM_IC"

    # Number of cells along x/y/z
    Nx = out.shape[0] 
    Ny = out.shape[1] 
    Nz = out.shape[2] 

    # Box size
    Lx = 10
    Ly = 10
    Lz = 10

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

    # Plot radius
    PlotRadius = 0.4*Lx**0.5
 
    ##########################
    ### Derived parameters ###
    ##########################
    Rho0_D, Sigma_D, Radius_D = Free2DerivedPara( Rho0_g, Sigma_g, Radius_g, Lambda, Kappa )


    ############################
    ###   Bundle paramters   ###
    ############################
    ParaPhy = [ Radius_g, Rho0_g, Sigma_g, Lambda, Kappa, Temp_g, Constant, Phi0, DevPhi0 ]
    ParaNum = [ BPoint, CoarseDr, PlotRadius ]
    


    # Number of cells along each dimension of the input grid.
    N = np.array([Nx, Ny, Nz], dtype='int')
    L = np.array([Lx, Ly, Lz], dtype='float64') 

    FluidInBox = SphericalSphere( L, N, ParaPhy, ParaNum )

    FluidInBox.tofile(filename)

out=np.zeros((512,512,512))
_DumpFile(out)
