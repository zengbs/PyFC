import numpy as np
from scipy.integrate import quad
import parameters as par
import density_profile
from scipy.special import kn


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

def Potential_Disk(R, z):
    def Integrand1( zp ):
        integrand1  = 0.5*par.Disk_alpha0 / par.Disk_z0 * np.exp( -np.absolute(zp)/par.Disk_z0 ) 
        integrand1 += 0.5*par.Disk_alpha1 / par.Disk_z1 * np.exp( -np.absolute(zp)/par.Disk_z1 ) 
        return integrand1
    def Integrand2( a ):
        PlusSqr     = np.sqrt( z**2 + ( a + R )**2 )
        MinuSqr     = np.sqrt( z**2 + ( a - R )**2 )
        integrand2  = a*kn(0, a/par.Disk_Rd)
        integrand2 *= np.arcsin( 2*a / ( PlusSqr + MinuSqr ) )
        return integrand2
    integral1 = quad(Integrand1, np.NINF, np.inf)[0]
    integral2 = quad(Integrand2,       0, np.inf)[0]
    Pot  = - 4*par.NEWTON_G*par.Disk_Sigma/par.Disk_Rd
    Pot *= integral1 * integral2
    return Pot


par.Parameters()
PotBulge    = Potential_Spheroidal(1, 1, density_profile.Bulge   , par.Bulge_q)
PotDarkHalo = Potential_Spheroidal(1, 1, density_profile.DarkHalo, par.Halo_q )
PotDisk     = Potential_Disk(1, 1)
print(PotBulge, PotDarkHalo, PotDisk)
