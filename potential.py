import numpy as np
from scipy.integrate import quad
import parameters as par
import density_profile


def Potential_Spheroidal (R, z, fun, c):
    a0 = 1
    def phi (mSqr, func):
        phi = quad(func, 0, mSqr)[0]
        return phi

    def SecondTerm( tau ):
        c0 = np.sqrt( 1 - e**2 )*a0
        mSqr  = R**2 / ( tau + a0**2 )
        mSqr += z**2 / ( tau + c0**2 )
        mSqr *= a0**2
        integrand  = phi(mSqr, fun)
        integrand /= tau + a0**2
        integrand /= np.sqrt(tau + c0**2)
        return integrand

    e = np.sqrt( 1 - c**2 )

    Pot  =  -2.0*np.pi*par.NEWTON_G * np.sqrt( 1 - e**2 ) / e

    Pot1 = phi(np.inf, fun)*np.arcsin(e)
  
    Pot2 = 0.5*a0*e*quad(SecondTerm, 0, np.inf)[0]

    Pot *= Pot1 - Pot2

    return Pot




par.Parameters()
PotBulge    = Potential_Spheroidal(1, 1, density_profile.Bulge   , par.Bulge_q)
PotDarkHalo = Potential_Spheroidal(1, 1, density_profile.DarkHalo, par.Halo_q )
print(PotBulge, PotDarkHalo)
