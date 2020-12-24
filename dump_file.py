import numpy as np
from pri2con import Pri2Con
from sphere import SphericalSphere
from density_profile import FreePara2DerivedPara
import parameters as par

def DumpFile():
    par.Parameters()

    ############################
    ###     Dump density    ####
    ############################
    FluidInBox, PotInBox = SphericalSphere()
    FileName = "%s" % (par.Case)+"_UM_IC"

    FluidInBox[0] *= par.Const_AtomMass*par.Const_MeanMolecularWeight 
    FluidInBox[1] *= par.Const_C*1e5*par.Const_AtomMass*par.Const_MeanMolecularWeight 
    FluidInBox[2] *= par.Const_C*1e5*par.Const_AtomMass*par.Const_MeanMolecularWeight 
    FluidInBox[3] *= par.Const_C*1e5*par.Const_AtomMass*par.Const_MeanMolecularWeight 
    FluidInBox[4] *= par.Const_AtomMass*par.Const_MeanMolecularWeight
    FluidInBox.tofile(FileName)

    ############################
    ###     Dump potential   ###
    ############################
    FileName = "%s" % (par.Case)+"_ExtPotTable"

    PotInBox *= par.Const_AtomMass*par.Const_MeanMolecularWeight*(par.Const_C*1e5)**2
    PotInBox.tofile(FileName)

par.Parameters()
out=np.zeros((par.Nx, par.Ny, par.Nz))
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
print("Rho0_g                    = %e" % par.Rho0_g                    )
print("Sigma_g                   = %e" % par.Sigma_g                   )
print("Lambda                    = %e" % par.Lambda                    )
print("Kappa                     = %e" % par.Kappa                     )
print("Phi0                      = %e" % par.Phi0                      )
print("DevPhi0                   = %e" % par.DevPhi0                   )
print("Const_AtomMass            = %e" % par.Const_AtomMass            )
print("Const_MeanMolecularWeight = %e" % par.Const_MeanMolecularWeight )
print("Rho0_D                    = %e" % par.Rho0_D                    )
print("Sigma_D                   = %e" % par.Sigma_D                   )
print("CoreRadius_D              = %e" % par.CoreRadius_D              )
print("BPoint                    = %e" % par.BPoint                    )
print("CoarseDr                  = %e" % par.CoarseDr                  )
print("Case                      = %s" % par.Case                      )
print("Precision                 = %s" % par.Precision                 )
print("Critical temperature      = %e" % par.CriticalTemp              )

if par.Case == "Mukherjee":
   print("Epsilon                = %e" % par.Epsilon                   )
   print("Temp_g                 = %e" % par.Temp_g                    )
   print("Sigma_t                = %e" % par.Sigma_t                   )
if par.Case == "Standard":
   print("a0                     = %e" % par.a0                        )
   print("z0                     = %e" % par.z0                        )
