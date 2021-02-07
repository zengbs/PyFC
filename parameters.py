import numpy as np
from fluid import Eta2Cs, Tem2Eta


def Parameters():

  global Const_kB, Const_C, NEWTON_G, Const_Erg2eV, Const_AtomicMassUnit, Const_MolecularWeight, Const_SolarMass, Const_MolecularWeightPerElectron
  global Nx, Ny, Nz, Lx, Ly, Lz, delta
  global dens_kmin, dens_mean, dens_sigma, dens_beta, dens_fromfile
  global Uxyz_kmin, Uxyz_mean, Uxyz_sigma, Uxyz_beta, Uxyz_fromfile
  global Precision, GRA_GHOST_SIZE
  global DensityRatio, Center
  global PeakElectronNumberDensity, Temperature, PeakGasMassDensity
  global V_halo, d_halo, DiskMass, a, b, BulgeMass, d_bulge
  global Cs, Eta
  global UNIT_D, UNIT_V, UNIT_L, UNIT_M, UNIT_P, UNIT_E, UNIT_T

  ##########################
  ###    Unit (cgs)      ###
  ##########################
  UNIT_D = 1e-24            # mass density
  UNIT_V = 29979245800.     # velocity
  UNIT_L = 3.08567758149e21 # length
  UNIT_T = UNIT_L / UNIT_V  # time
  UNIT_M = UNIT_D*UNIT_L**3 # mass
  UNIT_P = UNIT_D*UNIT_V**2 # energy density
  UNIT_E = UNIT_P*UNIT_L**3 # energy

  ##########################
  ### Physical constants ###
  ##########################

  # atomic mass unit (g)
  Const_AtomicMassUnit       = 1.66e-24

  # molecular weight
  Const_MolecularWeight      = 0.61

  # The Boltzmann constant (erg/K)
  Const_kB                   = 1.38064852e-16

  # Speed of light (cm/s)
  Const_C                    = 29979245800.0

  # Gravitational constant
  NEWTON_G                   = 6.6738e-8
  NEWTON_G                  /= UNIT_L**3
  NEWTON_G                  *= UNIT_M
  NEWTON_G                  *= UNIT_T**2

  # electron volt per erg
  Const_Erg2eV               = 6.2415e11

  # solar mass (g)
  Const_SolarMass            = 1.9885e33

  # kpc (cm)
  Const_kpc                  = 3.08567758149e21

  ##########################
  #######   PyFC  ##########
  ##########################

  # Fractal parameters for density
  # `None` stands for generating fractal cube by PyFC
  dens_fromfile = None
  dens_kmin     = 6.0
  dens_mean     = 1.0  
  dens_sigma    = np.sqrt(5.0)
  dens_beta     = -5.0/3.0

  # Fractal parameters for Ux/y/z
  # `None` stands for generating fractal cube by PyFC
  Uxyz_fromfile = None
  Uxyz_kmin     = 5.0
  Uxyz_mean     = 1.0 / Const_C 
  Uxyz_sigma    = 1000000.0 / Const_C
  Uxyz_beta     = -5.0/3.0

  ###############################
  ###  Physical parameters    ###
  ###  (normalized by UNIT)   ###
  ###############################

  # Box size (kpc)
  Lx = 50.
  Ly = 50.
  Lz = 50.

  Lx *= Const_kpc
  Ly *= Const_kpc
  Lz *= Const_kpc
  Lx /= UNIT_L
  Ly /= UNIT_L
  Lz /= UNIT_L

  # density ratio on both sides of the surface of the sphere
  DensityRatio = 45000.

  # peak electron number density (cm**-3)
  PeakElectronNumberDensity = 2.

  # initial temperature (K)
  Temperature = 2e6


  # *** gravitational potential ***

  # velocity dispersion of halo (km/s)
  V_halo  = 131.5
  V_halo *= 1e5 # km -> cm
  V_halo /= UNIT_V

  # distance (kpc)
  d_halo     = 12.

  d_halo *= Const_kpc
  d_halo /= UNIT_L

  # disk mass (solar mass)
  DiskMass = 1e11
  DiskMass *= Const_SolarMass
  DiskMass /= UNIT_M
  a = 6.50           # kpc
  b = 0.26           # kpc

  a *= Const_kpc
  b *= Const_kpc
  a /= UNIT_L     
  b /= UNIT_L     

  # bulge mass (solar mass)
  BulgeMass = 3.4e10
  BulgeMass *= Const_SolarMass
  BulgeMass /= UNIT_M

  d_bulge = 0.7      # kpc
  d_bulge *= Const_kpc
  d_bulge /= UNIT_L

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

  ###############################
  ###  Derived parameters     ###
  ###  (normalized by UNIT)   ###
  ###############################
  # Left edge and right edge coordinates of the desired
  # simulation domain which will be used in GAMER. (Normalized by CoreRadius_D)                                                        
  le = np.array([ 0.0,  0.0,  0.0])
  re = np.array([  Lx,   Ly,   Lz])
      
  # Cell spacing (normalized by CoreRadius_D)
  N      = np.array([Nx, Ny, Nz], dtype=Precision)
  if Precision is 'float32':
     delta  = (re-le)/N.astype(np.float32)
  else:
     delta  = (re-le)/N.astype(np.float64)

  # molecular weight per electron
  Const_MolecularWeightPerElectron = 5.0*Const_MolecularWeight/(2.0+Const_MolecularWeight)

  # peak gas mass density (g/cm**3)
  PeakGasMassDensity  = Const_MolecularWeightPerElectron*PeakElectronNumberDensity*Const_AtomicMassUnit
  PeakGasMassDensity /= UNIT_D

  # physical coordinate of center of sphere (origin is at corner) (kpc)
  Center = np.array([Lx,Ly,Lz])*0.5

  # sound speed
  Eta = Tem2Eta( Temperature )
  Cs = Eta2Cs( Eta )
  Cs /= UNIT_V
