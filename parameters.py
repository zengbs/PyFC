import numpy as np
from fluid import Eta2Cs, Tem2Eta


def Parameters():

  global Const_kB, Const_C, NEWTON_G, Const_Erg2eV, Const_AtomicMassUnit, Const_MolecularWeight, Const_SolarMass
  global Nx, Ny, Nz, Lx, Ly, Lz, delta
  global dens_kmin, dens_mean, dens_sigma, dens_beta, dens_fromfile
  global Uxyz_kmin, Uxyz_mean, Uxyz_sigma, Uxyz_beta, Uxyz_fromfile
  global Precision, GRA_GHOST_SIZE
  global DensityRatio, Center
  global Cs, Eta, AmbientTemperature
  global UNIT_D, UNIT_V, UNIT_L, UNIT_M, UNIT_P, UNIT_E, UNIT_T
  global Bulge_Rho0, Bulge_a, Bulge_r, Bulge_alpha, Bulge_q
  global Halo_Rho0, Halo_a, Halo_alpha, Halo_beta, Halo_q
  global Disk_Sigma, Disk_Rd, Disk_z0, Disk_z1, Disk_alpha0, Disk_alpha1
  global ISM_Sigma, ISM_Rg, ISM_Rm, ISM_zg

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
  Const_pc                   = 3.08567758149e18
  Const_kpc                  = 1e3*Const_pc

  ##########################
  #######   PyFC  ##########
  ##########################

  # Fractal parameters for density
  # `None` stands for generating fractal cube by PyFC
  dens_fromfile = 'off'
  dens_kmin     = 6.0
  dens_mean     = 1.0  
  dens_sigma    = np.sqrt(5.0)
  dens_beta     = -5.0/3.0

  # Fractal parameters for Ux/y/z
  # `None` stands for generating fractal cube by PyFC
  Uxyz_fromfile = 'off'
  Uxyz_kmin     = 1.0
  Uxyz_mean     = 1.0 / Const_C
  Uxyz_sigma    = 1e5 / Const_C
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

  AmbientTemperature = 1e10

  # Bulge
  Bulge_Rho0  = 0.427 # M_sun * pc**-3
  Bulge_a     = 1.0   # kpc
  Bulge_r     = 1.9   # kpc
  Bulge_alpha = 1.8
  Bulge_q     = 0.6

  Bulge_a    *= Const_kpc
  Bulge_r    *= Const_kpc
  Bulge_a    /= UNIT_L
  Bulge_r    /= UNIT_L

  Bulge_Rho0 *= Const_SolarMass*Const_pc**-3
  Bulge_Rho0 /= UNIT_D

  #  Dark halo
  Halo_Rho0   = 0.711 # M_sun * pc**-3
  Halo_a      = 3.83  # kpc
  Halo_alpha  = -2.0
  Halo_beta   = 2.96
  Halo_q      = 0.8

  Halo_a     *= Const_kpc
  Halo_a     /= UNIT_L

  Halo_Rho0  *= Const_SolarMass*Const_pc**-3
  Halo_Rho0  /= UNIT_D

  # Stellar disk
  Disk_Sigma  = 1905*0.75 # M_sun * pc**-2
  Disk_Rd     = 2.0       # kpc
  Disk_z0     = 0.3       # kpc
  Disk_z1     = 1.0       # kpc
  Disk_alpha0 = 0.5
  Disk_alpha1 = 1 - Disk_alpha0

  Disk_z0     *= Const_kpc
  Disk_z0     /= UNIT_L
  Disk_z1     *= Const_kpc
  Disk_z1     /= UNIT_L
  Disk_Rd     *= Const_kpc
  Disk_Rd     /= UNIT_L

  # ISM
  ISM_Sigma  = 1905-Disk_Sigma # M_sun * pc**-2
  ISM_Rg     = 2*Disk_Rd       # kpc
  ISM_Rm     = 4               # kpc
  ISM_zg     = 0.8             # kpc

  ISM_Rm     *= Const_kpc
  ISM_Rm     /= UNIT_L
  ISM_Rg     *= Const_kpc
  ISM_Rg     /= UNIT_L
  ISM_zg     *= Const_kpc
  ISM_zg     /= UNIT_L

  Disk_Sigma *= Const_SolarMass*Const_pc**-2
  Disk_Sigma /= UNIT_D*UNIT_L
  ISM_Sigma  *= Const_SolarMass*Const_pc**-2
  ISM_Sigma  /= UNIT_D*UNIT_L

  ############################
  ### Numerical parameters ###
  ############################
  # ghost zone size for potential
  GRA_GHOST_SIZE = 2

  # Floating-point precision
  Precision = 'float64'

  # Number of cells along x/y/z
  Nx = 32 
  Ny = 32 
  Nz = 32 
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


  # physical coordinate of center of sphere (origin is at corner) (kpc)
  Center = np.array([Lx,Ly,Lz])*0.5

  # sound speed
  Eta = Tem2Eta( AmbientTemperature )
  Cs = Eta2Cs( Eta )
  Cs /= UNIT_V
