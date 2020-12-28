import numpy as np


def Parameters():
  from density_profile import FreePara2DerivedPara

  global Nx, Ny, Nz, Lx, Ly, Lz, Const_kB, Const_C, Constant, NEWTON_G
  global CoreRadius_g, Rho0_g, Sigma_g, Lambda, Kappa, Epsilon, a0, z0
  global Temp_g, Phi0, DevPhi0, Const_AtomMass, Const_MeanMolecularWeight, Sigma_t, ISM0
  global BPoint, CoarseDr, SphereRadius, Precision, GRA_GHOST_SIZE
  global Rho0_D, Sigma_D, CoreRadius_D, CriticalTemp, Const_Erg2eV, Case
  global kmin, mean

  # Fractal parameters
  kmin = 10
  mean = 10  

  # Number of cells along x/y/z
  Nx = 512 
  Ny = 512  
  Nz = 512 
  N = np.array([Nx, Ny, Nz])

  if Nx % 16 is not 0 or Ny % 16 is not 0 or Nz % 16 is not 0:
     print("Nx/y/z % 16 != 0")
     exit()

  # Box size
  Lx = 3 
  Ly = 3 
  Lz = 3 

  ##########################
  ### Physical constants ###
  ##########################
  # Atomic mass unit (g)
  Const_AtomMass             = 1.66e-24

  # Mean molecular weight
  Const_MeanMolecularWeight  = 0.6

  # The Boltzmann constant (eV/K)
  Const_kB                   = 8.6173303e-5

  # Speed of light (cm/s)
  Const_C                    = 29979245800.0

  # Gravitational constant
  NEWTON_G                   = 1.0

  # electron volt per erg
  Const_Erg2eV               = 6.2415e11

  ##########################
  ###  Free parameters   ###
  ##########################
  # Sphere radius
  SphereRadius = 0.45*Lx

  # The eq(2) in Sutherland & Bicknell (2007)
  Constant  = 9

  # Core radius of gas sphere (kpc)
  CoreRadius_g = 1

  # Peak gas density (1/cm^3)
  Rho0_g = 0.5

  # Velocity dispersion of gas (cm/s)
  Sigma_g = 25000000

  # Lambda = r_{D}/r_{g}
  Lambda = 5

  # Kappa = \sigma_{D}/\sigma_{g}
  Kappa = 2

  # The potential at the center of sphere
  Phi0 = 0

  # The spatial derivative of potential at the center of sphere
  DevPhi0 = 0

  # Derived parameters
  Rho0_D, Sigma_D, CoreRadius_D = FreePara2DerivedPara( )

  ############################
  ###         Case         ###
  ### Mukherjee, Standard  ###
  ############################
  Case = "Mukherjee"  

  ### If Case = "Mukherjee" ###
  # Rotational coefficient
  Epsilon   = 0.93
  
  # Temperature of gas (K)
  Temp_g = 1e7
  #Temp_g = 2.7055e6 # this sound speed of this temperature is 25000000
 
  # Turbulent velocity dispersion (cm/s)
  # --> This parameter controls disk height
  Sigma_t = 20000000

  # Critical temperature for ISM disk (K)
  #CriticalTemp = 1e-99
  CriticalTemp = 3e4

  # number density at the center of ISM (1/cm^3)
  ISM0   = 200



  ## If Case = "Standard" ###
  #Case = "Standard"  
  #a0 = 1
  #z0 = 0.1
  #Sigma_t = 25000000 # (cm/s)
  ## Critical temperature for ISM disk (K)
  ##CriticalTemp = 4e6
  #CriticalTemp = 1e-99
  ##CriticalTemp = 1e+99


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
