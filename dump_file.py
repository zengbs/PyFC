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
    #FluidInBox, PotInBox = SphericalSphere()
    #FluidInBox.tofile("UM_IC")
    ISM, PotInBox = SphericalSphere()
    ISM.tofile("UM_IC")

    ############################
    ###     Dump potential   ###
    ############################
    PotInBox.tofile("ExtPotTable")

par.Parameters()
out=np.zeros((par.Nx, par.Ny, par.Nz))
DumpFile()

print("Nx        = %d" % par.Nx        )
print("Ny        = %d" % par.Ny        )
print("Nz        = %d" % par.Nz        )
print("Lx        = %d" % par.Lx        )
print("Ly        = %d" % par.Ly        )
print("Lz        = %d" % par.Lz        )
print("kB        = %e" % par.kB        )
print("C         = %e" % par.C         )
print("NEWTON_G  = %e" % par.NEWTON_G  )
print("Radius    = %e" % par.Radius    )
print("Constant  = %e" % par.Constant  )
print("Radius_g  = %e" % par.Radius_g  )
print("Rho0_g    = %e" % par.Rho0_g    )
print("Sigma_g   = %e" % par.Sigma_g   )
print("Lambda    = %e" % par.Lambda    )
print("Kappa     = %e" % par.Kappa     )
print("Temp_g    = %e" % par.Temp_g    )
print("Phi0      = %e" % par.Phi0      )
print("DevPhi0   = %e" % par.DevPhi0   )
print("Matom     = %e" % par.Matom     )
print("Mu        = %e" % par.Mu        )
print("Rho0_D    = %e" % par.Rho0_D    )
print("Sigma_D   = %e" % par.Sigma_D   )
print("Radius_D  = %e" % par.Radius_D  )
print("BPoint    = %e" % par.BPoint    )
print("CoarseDr  = %e" % par.CoarseDr  )
print("Precision = %s" % par.Precision )
