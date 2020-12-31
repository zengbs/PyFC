import numpy as np
import fluid


def Parameters():
  from density_profile import FreePara2DerivedPara

  global Const_kB, Const_C, NEWTON_G, Const_Erg2eV, Const_AtomicMassUnit, Const_MolecularWeight, Const_SolarMass, Const_MolecularWeightPerElectron
  global Nx, Ny, Nz, Lx, Ly, Lz, delta
  global dens_kmin, dens_mean, dens_sigma, dens_beta, dens_fromfile
  global Uxyz_kmin, Uxyz_mean, Uxyz_sigma, Uxyz_beta, Uxyz_fromfile
  global Precision, GRA_GHOST_SIZE
  global SphereRadius, DensRatio
  global PeakElectronNumberDensity, Temperature, PotCenter, DiffPotCenter, PeakGasNumberDensity
  global V_halo, d_halo, DiskMass, a, b, BulgeMass, d_bulge
  global Cs


  ##########################
  ### Physical constants ###
  ##########################

  # atomic mass unit (g)
  Const_AtomicMassUnit       = 1.66e-24

  # molecular weight
  Const_MolecularWeight      = 0.61

  # The Boltzmann constant (eV/K)
  Const_kB                   = 8.6173303e-5

  # Speed of light (cm/s)
  Const_C                    = 29979245800.0

  # Gravitational constant
  NEWTON_G                   = 1.0

  # electron volt per erg
  Const_Erg2eV               = 6.2415e11

  # solar mass (g)
  Const_SolarMass            = 1.9891e33

  ##########################
  #######   PyFC  ##########
  ##########################

  # Fractal parameters for density
  # `None` stands for generating fractal cube by PyFC
  dens_fromfile = None
  dens_kmin  = 6.0
  dens_mean  = 1.0  
  dens_sigma = np.sqrt(5.0)
  dens_beta  = -5.0/3.0

  # Fractal parameters for Ux/y/z
  # `None` stands for generating fractal cube by PyFC
  Uxyz_fromfile = None
  Uxyz_kmin  = 3.0
  Uxyz_mean  = 1.0  
  Uxyz_sigma = 10000000.0 / np.sqrt(3.0) / Const_C
  Uxyz_beta  = -5.0/3.0

  ###############################
  ###  Physical parameters    ###
  ###############################

  # Box size (kpc)
  Lx = 4.0 
  Ly = 4.0 
  Lz = 4.0 

  # Sphere radius
  SphereRadius = 0.45*Lx

  # density ratio on both sides of the surface of the sphere
  DensRatio = 1.

  # peak electron number density (cm**-3)
  PeakElectronNumberDensity   = 2.

  # initial temperature (K)
  Temperature = 2e6

  # potential at the center of sphere
  PotCenter = 0.

  # spatial derivative of potential at the center of sphere
  DiffPotCenter = 0.



  # *** gravitational potential ***

  # velocity dispersion of halo (km/s)
  V_halo            = 131.5

  # distance (kpc)
  d_halo     = 12.

  # disk mass (solar mass)
  DiskMass = 1e11
  a = 6.50           # kpc
  b = 0.26           # kpc

  # bulge mass (solar mass)
  BulgeMass = 3.4e10
  d_bulge = 0.7      # kpc


  ############################
  ### Numerical parameters ###
  ############################
  # ghost zone size for potential
  GRA_GHOST_SIZE = 2

  # Floating-point precision
  Precision = 'float32'

  # Number of cells along x/y/z
  Nx = 512 
  Ny = 512 
  Nz = 512 
  N = np.array([Nx, Ny, Nz])

  if Nx % 16 is not 0 or Ny % 16 is not 0 or Nz % 16 is not 0:
     print("Nx/y/z % 16 != 0")
     exit()

  ############################
  # Derived parameters
  ############################
  # Left edge and right edge coordinates of the desired
  # simulation domain which will be used in GAMER. (Normalized by CoreRadius_D)                                                        
  le = np.array([ 0.0,  0.0,  0.0])
  re = np.array([par.Lz, par.Ly, par.Lx]) / float(par.CoreRadius_D)
      
  # Cell spacing (normalized by CoreRadius_D)
  N = np.array([par.Nz, par.Ny, par.Nx ], dtype=par.Precision)
  delta = (re-le)/N

  # molecular weight per electron
  Const_MolecularWeightPerElectron = 5.0*Const_MolecularWeight/(2.0+Const_MolecularWeight)

  # peak gas number density (cm**-3)
  PeakGasNumberDensity   = PeakElectronNumberDensity*Const_MolecularWeightPerElectron*Const_AtomicMassUnit

  # physical coordinate of center of sphere (origin is at corner) (kpc)
  Center = np.array([par.Lz,par.Ly,par.Lx])*0.5

  # sound speed
  Cs = Tem2Cs(par.Temperature)
