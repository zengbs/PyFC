import numpy as np


print("Computing gas density ...")
sys.stdout.flush()
Rho                      = TotGasDensity()
print("Computing gas density ... done !!")
sys.stdout.flush()

GradientRho_axis0        = np.gradient(Rho, axis=0)
GradientRho_axis1        = np.gradient(Rho, axis=1)
GradientRho_axis2        = np.gradient(Rho, axis=2)


print("Computing potential ...")
sys.stdout.flush()
Pot                      = TotPotential()
print("Computing potential ... done !!")
sys.stdout.flush()

GradientPot_axis0        = np.gradient(Pot, axis=0)
GradientPot_axis1        = np.gradient(Pot, axis=1)
GradientPot_axis2        = np.gradient(Pot, axis=2)


PoissonSource            = -4*np.pi*par.NEWTON_G*Rho**2
PoissonSource           -= GradientRho_GradientPot
PoissonSource           -= GradientRho_axis0*GradientPot_axis0
PoissonSource           -= GradientRho_axis1*GradientPot_axis1
PoissonSource           -= GradientRho_axis2*GradientPot_axis2
