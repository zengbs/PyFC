
def Parameters():
  from density_profile import Free2DerivedPara

  global Nx, Ny, Nz, Lx, Ly, Lz, kB, C, Constant, NEWTON_G
  global Radius_g, Rho0_g, Sigma_g, Lambda, Kappa
  global Temp_g, Phi0, DevPhi0, Matom, Mu 
  global BPoint, CoarseDr, Radius, Precision
  global Rho0_D, Sigma_D, Radius_D

  
  # Number of cells along x/y/z
  Nx = 64 
  Ny = 64  
  Nz = 64  

  # Box size
  Lx = 3
  Ly = 3
  Lz = 3

  ##########################
  ### Physical constants ###
  ##########################
  # The Boltzmann constant (eV/K)
  kB = 8.6173303e-5

  # Speed of light (km/s)
  C = 299792.458

  # Gravitational constant
  NEWTON_G = 1.0

 
  ##########################
  ###  Free parameters   ###
  ##########################
  # Sphere radius
  Radius = 0.4*Lx

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

  # Atomic mass unit (g)
  Matom = 1.66e-24

  # Mean molecular weight
  Mu = 0.6

  # Derived parameters
  Rho0_D, Sigma_D, Radius_D = Free2DerivedPara( Rho0_g, Sigma_g, Radius_g, Lambda, Kappa )

  ############################
  ### Numerical parameters ###
  ############################

  # The boundary point for integration
  BPoint   = 1e-5

  # The integration step
  CoarseDr = 1e-5

  # Floating-point precision
  Precision = 'float32'
