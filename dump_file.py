import numpy as np
from sphere import SetIC
import parameters as par
import pyFC
import os
import sys
from multiprocessing import Pool
from threading import Thread

par.Parameters()

def GetFractalDensity(x):
    fc = pyFC.LogNormalFractalCube(ni=par.Nx, nj=par.Ny, nk=par.Nz,
                                   kmin=par.dens_kmin, mean=par.dens_mean,
                                   sigma=par.dens_sigma, beta=par.dens_beta)
    fc.gen_cube()                           
    FractalDensity = pyFC.write_cube(fc=fc, app=True, prec='single')
    FractalDensity.tofile("FractalDensity")
    return FractalDensity


    # Ux
def GetFractalUx(x):
    fc = pyFC.LogNormalFractalCube(ni=par.Nx, nj=par.Ny, nk=par.Nz,
                                   kmin=par.Uxyz_kmin, mean=par.Uxyz_mean,
                                   sigma=par.Uxyz_sigma, beta=par.Uxyz_beta)
    fc.gen_cube()                           
    FractalUx = pyFC.write_cube(fc=fc, app=True, prec='single')
    FractalUx.tofile("FractalUx")
    return FractalUx

    # Uy
def GetFractalUy(x):
    fc = pyFC.LogNormalFractalCube(ni=par.Nx, nj=par.Ny, nk=par.Nz,
                                   kmin=par.Uxyz_kmin, mean=par.Uxyz_mean,
                                   sigma=par.Uxyz_sigma, beta=par.Uxyz_beta)
    fc.gen_cube()                           
    FractalUy = pyFC.write_cube(fc=fc, app=True, prec='single')
    FractalUy.tofile("FractalUy")
    return FractalUy

    # Uz
def GetFractalUz(x):
    fc = pyFC.LogNormalFractalCube(ni=par.Nx, nj=par.Ny, nk=par.Nz,
                                   kmin=par.Uxyz_kmin, mean=par.Uxyz_mean,
                                   sigma=par.Uxyz_sigma, beta=par.Uxyz_beta)
    fc.gen_cube()                           
    FractalUz = pyFC.write_cube(fc=fc, app=True, prec='single')
    FractalUz.tofile("FractalUz")
    return FractalUz

pool = Pool(processes=4)

############################
###   Fractal denisty   ####
############################
if par.dens_fromfile == None:
   result1 = pool.map_async(GetFractalDensity,range(1))

############################
###   Fractal Ux/y/z    ####
############################

if par.Uxyz_fromfile == None:
   result2 = pool.map_async(GetFractalUx,range(1))
   result3 = pool.map_async(GetFractalUy,range(1))
   result4 = pool.map_async(GetFractalUz,range(1))

if par.dens_fromfile == None:
   FractalDensity = np.array(result1.get())
elif par.dens_fromfile == 'off':
   FractalDensity = np.full( [par.Nx, par.Ny, par.Nz ], 1.0, dtype=par.Precision )

if par.Uxyz_fromfile == None:
   FractalUx      = np.array(result2.get())
   FractalUy      = np.array(result3.get())
   FractalUz      = np.array(result4.get())
elif par.Uxyz_fromfile == 'off':
   FractalUx      = np.full( FractalDensity.shape, 0.0, dtype=par.Precision )
   FractalUy      = np.full( FractalDensity.shape, 0.0, dtype=par.Precision )
   FractalUz      = np.full( FractalDensity.shape, 0.0, dtype=par.Precision )

FractalDensity.tofile("FractalDensity")
FractalUx.tofile("FractalUx")
FractalUy.tofile("FractalUy")
FractalUz.tofile("FractalUz")

print( "The varience of the fractal density = %e" % np.var (FractalDensity) )
print( "The     mean of the fractal density = %e" % np.mean(FractalDensity) )
print( "The varience of the fractal Ux      = %e" % np.var (FractalUx     ) )
print( "The     mean of the fractal Ux      = %e" % np.mean(FractalUx     ) )
print( "The varience of the fractal Uy      = %e" % np.var (FractalUy     ) )
print( "The     mean of the fractal Uy      = %e" % np.mean(FractalUy     ) )
print( "The varience of the fractal Uz      = %e" % np.var (FractalUz     ) )
print( "The     mean of the fractal Uz      = %e" % np.mean(FractalUz     ) )

sys.stdout.flush()


#############################
####     Dump density    ####
#############################
GasDens,   GasMomX,  GasMomY,  GasMomZ,  GasEngy, Potential = SetIC( FractalDensity, FractalUx, FractalUy, FractalUz )

if par.Precision == 'float64':
   FileName = "UM_IC_double"
elif par.Precision == 'float32':
   FileName = "UM_IC_single"

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
if par.Precision == 'float64':
   FileName = "ExtPotTable_double"
elif par.Precision == 'float32':
   FileName = "ExtPotTable_float"

Potential *= par.UNIT_V**2
Potential.tofile(FileName)


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
print("kmin                             = %e" % par.dens_kmin                       )
print("mean                             = %e" % par.dens_mean                       )
print("sigma                            = %e" % par.dens_sigma                      )
print("beta                             = %e" % par.dens_beta                       )
print("kmin                             = %e" % par.Uxyz_kmin                       )
print("mean                             = %e" % par.Uxyz_mean                       )
print("sigma                            = %e" % par.Uxyz_sigma                      )
print("beta                             = %e" % par.Uxyz_beta                       )
print("UNIT_D                           = %e" % par.UNIT_D                          )
print("UNIT_V                           = %e" % par.UNIT_V                          )
print("UNIT_L                           = %e" % par.UNIT_L                          )
print("UNIT_T                           = %e" % par.UNIT_T                          )
print("UNIT_M                           = %e" % par.UNIT_M                          )
print("UNIT_P                           = %e" % par.UNIT_P                          )
print("UNIT_E                           = %e" % par.UNIT_E                          )
print("TrunDensityRatio                 = %e" % par.TrunDensityRatio                )
print("FracDensityRatio                 = %e" % par.TrunDensityRatio                )
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
