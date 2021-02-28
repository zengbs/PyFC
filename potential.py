import numpy as np
from scipy.integrate import quad, tplquad
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

def Check(x):
    if not np.isfinite(x):
       print("Inf!!!")
       return True
    if np.isnan(x):
       print("NaN!!!")
       return True


def Potential_ThickDisk(R, z, Sigma, Zeta):
    def OutterIntegral(zp):
        def MediumIntegral( a ):
            def InnerIntegral( Rp ):
                inner  = Rp*Sigma(Rp)
                inner /= np.sqrt( Rp**2 - a**2 )
                return inner
            zpp = z - zp
            if zpp != 0:
               PlusSqr = np.sqrt( zpp**2 + ( a + R )**2 )
               MinuSqr = np.sqrt( zpp**2 + ( a - R )**2 )
               medium  =  ( a + R ) / PlusSqr
               medium -=  ( a - R ) / MinuSqr
               medium /=  np.sqrt( R**2 - zpp**2 - a**2 + PlusSqr*MinuSqr )
            else:
               if a > R:
                  medium = 0
               elif a < R:
                  medium = np.sqrt(2/(R**2-a**2))
               else:
                  print(a, R, zpp)
                  exit(0)
            medium *= quad(InnerIntegral, a, np.inf)[0]
            medium *= -2**1.5*par.NEWTON_G
            return medium
        outter  = quad(MediumIntegral, 0, np.inf)[0]
        outter *= Zeta(zp)
        return outter
     
    Pot = quad(OutterIntegral, -np.inf, np.inf)[0]
    return Pot



def Potential_Disk(R, z):
    def Zeta(z):
        Zeta  = 0.5*par.Disk_alpha0 / par.Disk_z0 * np.exp( -np.absolute(z)/par.Disk_z0 ) 
        Zeta += 0.5*par.Disk_alpha1 / par.Disk_z1 * np.exp( -np.absolute(z)/par.Disk_z1 ) 
        return Zeta
    def Sigma(R):
        return par.Disk_Sigma * np.exp( -R/par.Disk_Rd )
    return Potential_ThickDisk(R, z, Sigma, Zeta)

def Potential_ISM(R, z):
    def Zeta(z):
        return 0.5*par.ISM_Sigma/par.ISM_zg * np.exp(- np.absolute(z)/par.ISM_zg)
    def Sigma(R):
        return np.exp( -R/par.ISM_Rg - par.ISM_Rm / R )
    return Potential_ThickDisk(R, z, Sigma, Zeta)
   

par.Parameters()
#PotBulge    = Potential_Spheroidal(1, 1, density_profile.Bulge   , par.Bulge_q)
#PotDarkHalo = Potential_Spheroidal(1, 1, density_profile.DarkHalo, par.Halo_q )
#print(PotBulge, PotDarkHalo)
PotDisk     = Potential_Disk(10, 1)
#PotISM      = Potential_ISM(3, 2)
print(PotDisk)
