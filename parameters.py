import numpy as np
from fluid import Eta2Cs, Tem2Eta


def Parameters():

  global Const_kB, Const_C, NEWTON_G, Const_Erg2eV, Const_AtomicMassUnit, Const_MolecularWeight, Const_SolarMass, Const_MolecularWeightPerElectron
  global Nx, Ny, Nz, Lx, Ly, Lz, delta
  global Precision, AmbientDens
  global dens_fromfile, dens_kmin, dens_mean, dens_sigma, dens_beta
  global Uxyz_fromfile, Uxyz_kmin, Uxyz_mean, Uxyz_sigma, Uxyz_beta  
  global Center
  global AmbientTemperature
  global Cs, Eta
  global UNIT_D, UNIT_V, UNIT_L, UNIT_M, UNIT_P, UNIT_E, UNIT_T
  global dens_FractalOn, Uxyz_FractalOn

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
  dens_FractalOn = True
  dens_fromfile  = None
  dens_kmin      = 16.66
  dens_mean      = 1.
  dens_sigma     = np.sqrt(5.0)
  dens_beta      = -5.0/3.0

  # Fractal parameters for Ux/y/z
  # `None` stands for generating fractal cube by PyFC
  Uxyz_FractalOn = False
  Uxyz_fromfile  = None
  Uxyz_kmin      = 3.0
  Uxyz_mean      = 1.0  
  Uxyz_sigma     = 10000000.0 / np.sqrt(3.0) / Const_C
  Uxyz_beta      = -5.0/3.0

  ###############################
  ###  Physical parameters    ###
  ###  (normalized by UNIT)   ###
  ###############################

  # Box size (kpc)
  Lx = 5.0 
  Ly = 5.0 
  Lz = 5.0 

  Lx *= Const_kpc
  Ly *= Const_kpc
  Lz *= Const_kpc
  Lx /= UNIT_L
  Ly /= UNIT_L
  Lz /= UNIT_L


  # initial temperature (K)
  AmbientTemperature = 1e7
  AmbientDens = 1e-25
  AmbientDens /= UNIT_D

  ############################
  ### Numerical parameters ###
  ############################
  # Floating-point precision
  Precision = 'float32'

  # Number of cells along x/y/z
  Nx = 256 
  Ny = 256 
  Nz = 256 
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


  # physical coordinate of center of sphere (origin is at corner) (kpc)
  Center = np.array([Lx,Ly,Lz])*0.5

  # sound speed
  Eta = Tem2Eta( AmbientTemperature )
  Cs = Eta2Cs( Eta )
  Cs /= UNIT_V
