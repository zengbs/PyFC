import numpy as np
from sphere import SphericalSphere
from density_profile import FreePara2DerivedPara
import parameters as par
import pyFC
import os
import sys


def DumpFile():
    par.Parameters()

    ############################
    ###   Fractal denisty   ####
    ############################

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


    ############################
    ###   Fractal Ux/y/z    ####
    ############################

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

    #FractalDensity = 1
    #FractalUxyz    = 1


    FractalU = np.power(FractalUxyz, 2.0)*3

    print( "The varience of the fractal density = %e" % np.var (FractalDensity) )
    print( "The     mean of the fractal density = %e" % np.mean(FractalDensity) )
    print( "The varience of the fractal Uxyz    = %e" % np.var (FractalUxyz   ) )
    print( "The     mean of the fractal Uxyz    = %e" % np.mean(FractalUxyz   ) )
    print( "The varience of the fractal U       = %e" % np.var (FractalU      ) )
    print( "The     mean of the fractal U       = %e" % np.mean(FractalU      ) )

    sys.stdout.flush()


    #############################
    ####     Dump density    ####
    #############################
    FluidInBox, PotInBox = SphericalSphere( FractalDensity, FractalUxyz )

    FileName = "UM_IC"

    FluidInBox[0] *= par.Const_AtomicMassUnit*par.Const_MolecularWeight 
    FluidInBox[1] *= par.Const_C*par.Const_AtomicMassUnit*par.Const_MolecularWeight 
    FluidInBox[2] *= par.Const_C*par.Const_AtomicMassUnit*par.Const_MolecularWeight 
    FluidInBox[3] *= par.Const_C*par.Const_AtomicMassUnit*par.Const_MolecularWeight 
    FluidInBox[4] *= par.Const_AtomicMassUnit*par.Const_MolecularWeight*par.Const_C**2

    FluidInBox.tofile(FileName)

    #############################
    ####     Dump potential   ###
    #############################
    FileName = "ExtPotTable"

    PotInBox *= par.Const_AtomicMassUnit*par.Const_MolecularWeight*par.Const_C**2
    PotInBox.tofile(FileName)

par.Parameters()
DumpFile()

print("Nx                        = %d" % par.Nx                        )
print("Ny                        = %d" % par.Ny                        )
print("Nz                        = %d" % par.Nz                        )
print("Lx                        = %d" % par.Lx                        )
print("Ly                        = %d" % par.Ly                        )
print("Lz                        = %d" % par.Lz                        )
print("Const_kB                  = %e" % par.Const_kB                  )
print("Const_C                   = %e" % par.Const_C                   )
print("NEWTON_G                  = %e" % par.NEWTON_G                  )
print("SphereRadius              = %e" % par.SphereRadius              )
print("Constant                  = %e" % par.Constant                  )
print("CoreRadius_g              = %e" % par.CoreRadius_g              )
print("PeakNumberDensity                    = %e" % par.PeakNumberDensity                    )
print("Sigma_g                   = %e" % par.Sigma_g                   )
print("Lambda                    = %e" % par.Lambda                    )
print("Kappa                     = %e" % par.Kappa                     )
print("PotCenter                      = %e" % par.PotCenter                      )
print("DiffPotCenter                   = %e" % par.DiffPotCenter                   )
print("Const_AtomicMassUnit            = %e" % par.Const_AtomicMassUnit            )
print("Const_MolecularWeight = %e" % par.Const_MolecularWeight )
print("Rho0_D                    = %e" % par.Rho0_D                    )
print("Sigma_D                   = %e" % par.Sigma_D                   )
print("CoreRadius_D              = %e" % par.CoreRadius_D              )
print("BPoint                    = %e" % par.BPoint                    )
print("CoarseDr                  = %e" % par.CoarseDr                  )
print("Case                      = %s" % par.Case                      )
print("Precision                 = %s" % par.Precision                 )
print("Critical temperature      = %e" % par.CriticalTemp              )
print("kmin                      = %e" % par.dens_kmin                 )
print("mean                      = %e" % par.dens_mean                 )
print("sigma                     = %e" % par.dens_sigma                )
print("beta                      = %e" % par.dens_beta                 )
print("kmin                      = %e" % par.Uxyz_kmin                 )
print("mean                      = %e" % par.Uxyz_mean                 )
print("sigma                     = %e" % par.Uxyz_sigma                )
print("beta                      = %e" % par.Uxyz_beta                 )

if par.Case == "Mukherjee":
   print("Epsilon                   = %e" % par.Epsilon                )
   print("Temp_g                    = %e" % par.Temp_g                 )
   print("Sigma_t                   = %e" % par.Sigma_t                )
if par.Case == "Standard":
   print("a0                        = %e" % par.a0                     )
   print("z0                        = %e" % par.z0                     )
