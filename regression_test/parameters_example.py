
def Parameters():
  from density_profile import FreePara2DerivedPara

  global Nx, Ny, Nz, Lx, Ly, Lz, kB, C, Constant, NEWTON_G
  global CoreRadius_g, Rho0_g, Sigma_g, Lambda, Kappa
  global Temp_g, Phi0, DevPhi0, Matom, Mu, z0, a0
  global BPoint, CoarseDr, SphereRadius, Precision, GRA_GHOST_SIZE
  global Rho0_D, Sigma_D, CoreRadius_D, CriticalTemp, Const_Erg2eV

  
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
  kB       = 8.6173303e-5

  # Speed of light (km/s)
  C        = 299792.458

  # Gravitational constant
  NEWTON_G = 1.0

  # electron volt
  Const_Erg2eV = 6.2415e11
  ##########################
  ###  Free parameters   ###
  ##########################
  # Sphere radius
  SphereRadius = 0.4*Lx

  # The eq(2) in Sutherland & Bicknell (2007)
  Constant = 9

  # Core radius of gas sphere (kpc)
  CoreRadius_g = 1

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

  # Scale height of ISM disk
  z0   = 0.3

  # Rotational coefficient
  a0   = 1

  # Critical temperature for ISM disk (K)
  #CriticalTemp = 1e5
  CriticalTemp = 3e4

  # Derived parameters
  Rho0_D, Sigma_D, CoreRadius_D = FreePara2DerivedPara( )

  ############################
  ### Numerical parameters ###
  ############################

  # The boundary point for integration
  BPoint   = 1e-5

  # The integration step
  CoarseDr = 1e-5

  # ghost zone size for potential
  GRA_GHOST_SIZE = 2

  # Floating-point precision
  Precision = 'float32'
