import numpy as np
from sphere import SetIC
import parameters as par
import pyFC
import os
import sys


def DumpFile():

    ############################
    ###   Fractal denisty   ####
    ############################
    
    if par.dens_FractalOn:
       if par.dens_fromfile is None:
          fc = pyFC.LogNormalFractalCube(ni=par.Nx, nj=par.Ny, nk=par.Nz,
                                         kmin=par.dens_kmin, mean=par.dens_mean,
                                         sigma=par.dens_sigma, beta=par.dens_beta)
          fc.gen_cube()                           
          FractalDensity = pyFC.write_cube(fc=fc, app=True, prec='single')
          FractalDensity.tofile("FractalDensity")
       else:
          if  os.path.isfile("FractalDensity"):
              FractalDensity = np.fromfile("FractalDensity",dtype=par.Precision)
              FractalDensity = FractalDensity.reshape(par.Nx,par.Ny,par.Nz)
          else:
              print("FractalDensity does not exist!!")
              exit(0)
    else:
       FractalDensity = None


    ############################
    ###   Fractal Ux/y/z    ####
    ############################

    if par.Uxyz_FractalOn:
       if par.Uxyz_fromfile is None:
          fc = pyFC.LogNormalFractalCube(ni=par.Nx, nj=par.Ny, nk=par.Nz,
                                         kmin=par.Uxyz_kmin, mean=par.Uxyz_mean,
                                         sigma=par.Uxyz_sigma, beta=par.Uxyz_beta)
          fc.gen_cube()                           
          FractalUxyz = pyFC.write_cube(fc=fc, app=True, prec='single')
          FractalUxyz.tofile("FractalUxyz")
       else:
          if  os.path.isfile("FractalUxyz"):
              FractalUxyz = np.fromfile("FractalUxyz",dtype=par.Precision)
              FractalUxyz = FractalUxyz.reshape(par.Nx,par.Ny,par.Nz)
          else:
              print("FractalUxyz does not exist!!")
              exit(0)


       FractalU = np.power(FractalUxyz, 2.0)*3

       print( "The varience of the fractal density = %e" % np.var (FractalDensity) )
       print( "The     mean of the fractal density = %e" % np.mean(FractalDensity) )
       print( "The varience of the fractal Uxyz    = %e" % np.var (FractalUxyz   ) )
       print( "The     mean of the fractal Uxyz    = %e" % np.mean(FractalUxyz   ) )
       print( "The varience of the fractal U       = %e" % np.var (FractalU      ) )
       print( "The     mean of the fractal U       = %e" % np.mean(FractalU      ) )

       sys.stdout.flush()
    else:
       FractalUxyz = None


    #############################
    ####     Dump density    ####
    #############################

    GasDens,   GasMomX,  GasMomY,  GasMomZ,  GasEngy, Potential = SetIC( FractalDensity, FractalUxyz )

    FileName = "UM_IC"

    Fluid3D = np.zeros((5, par.Nx, par.Ny, par.Nz),dtype=par.Precision)

    Fluid3D[0] = GasDens * par.UNIT_D
    Fluid3D[1] = GasMomX * par.UNIT_D * par.UNIT_V  
    Fluid3D[2] = GasMomY * par.UNIT_D * par.UNIT_V
    Fluid3D[3] = GasMomZ * par.UNIT_D * par.UNIT_V
    Fluid3D[4] = GasEngy * par.UNIT_P

    Fluid3D.tofile(FileName)

    #############################
    ####     Dump potential   ###
    #############################
    FileName = "ExtPotTable"

    Potential *= par.UNIT_V**2
    Potential.tofile(FileName)

par.Parameters()
DumpFile()

print("Nx                               = %d" % par.Nx                              )
print("Ny                               = %d" % par.Ny                              )
print("Nz                               = %d" % par.Nz                              )
print("Lx                               = %d" % par.Lx                              )
print("Ly                               = %d" % par.Ly                              )
print("Lz                               = %d" % par.Lz                              )
print("Const_kB                         = %e" % par.Const_kB                        )
print("Const_C                          = %e" % par.Const_C                         )
print("NEWTON_G                         = %e" % par.NEWTON_G                        )
print("Const_AtomicMassUnit             = %e" % par.Const_AtomicMassUnit            )
print("Const_MolecularWeight            = %e" % par.Const_MolecularWeight           )
print("Precision                        = %s" % par.Precision                       )
print("dens_kmin                        = %e" % par.dens_kmin                       )
print("dens_mean                        = %e" % par.dens_mean                       )
print("dens_sigma                       = %e" % par.dens_sigma                      )
print("dens_beta                        = %e" % par.dens_beta                       )
print("Uxyz_kmin                        = %e" % par.Uxyz_kmin                       )
print("Uxyz_mean                        = %e" % par.Uxyz_mean                       )
print("Uxyz_sigma                       = %e" % par.Uxyz_sigma                      )
print("Uxyz_beta                        = %e" % par.Uxyz_beta                       )
print("UNIT_D                           = %e" % par.UNIT_D                          )
print("UNIT_V                           = %e" % par.UNIT_V                          )
print("UNIT_L                           = %e" % par.UNIT_L                          )
print("UNIT_T                           = %e" % par.UNIT_T                          )
print("UNIT_M                           = %e" % par.UNIT_M                          )
print("UNIT_P                           = %e" % par.UNIT_P                          )
print("UNIT_E                           = %e" % par.UNIT_E                          )
print("DensityRatio                     = %e" % par.DensityRatio                    )
print("PeakElectronNumberDensity        = %e" % par.PeakElectronNumberDensity       )
print("Temperature                      = %e" % par.Temperature                     )
print("V_halo                           = %e" % par.V_halo                          )
print("d_halo                           = %e" % par.d_halo                          )
print("DiskMass                         = %e" % par.DiskMass                        )
print("a                                = %e" % par.a                               )
print("b                                = %e" % par.b                               )
print("BulgeMass                        = %e" % par.BulgeMass                       )
print("d_bulge                          = %e" % par.d_bulge                         )
print("Const_MolecularWeightPerElectron = %e" % par.Const_MolecularWeightPerElectron)
print("PeakGasMassDensity               = %e" % par.PeakGasMassDensity              )
print("CenterX                          = %e" % par.Center[0]                       )
print("CenterY                          = %e" % par.Center[1]                       )
print("CenterZ                          = %e" % par.Center[2]                       )
print("Cs                               = %e" % par.Cs                              )
