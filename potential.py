import numpy as np
from scipy.integrate import quad
import parameters as par

import density_profile

def phi (mSqr):
    phi = quad(density_profile.Bulge, 0, mSqr)[0]
    return phi

def Pot_Bulge (R, z):
    a0 = 1
    def SecondTerm( tau ):
        c0 = np.sqrt( 1 - e**2 )*a0
        mSqr  = R**2 / ( tau + a0**2 )
        mSqr += z**2 / ( tau + c0**2 )
        mSqr *= a0**2
        integrand  = phi(mSqr)
        integrand /= tau + a0**2
        integrand /= np.sqrt(tau + c0**2)
        return integrand

    e = np.sqrt( 1 - par.Bulge_q**2 )

    Pot  =  -2.0*np.pi*par.NEWTON_G * np.sqrt( 1 - e**2 ) / e

    Pot_1 = phi(np.inf)*np.arcsin(e)
  
    Pot_2 = 0.5*a0*e*quad(SecondTerm, 0, np.inf)[0]

    Pot *= Pot_1 - Pot_2

    return Pot


par.Parameters()
Pot = Pot_Bulge(1, 1)
print(Pot)
