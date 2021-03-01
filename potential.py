import numpy as np
from scipy.integrate import quad, tplquad
from scipy.misc import derivative
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
       return True
    if np.isnan(x):
       return True


def Potential_ThickDisk1(R, z, Sigma, Zeta):
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
               if Check(medium):
                  return np.nan
            else:
               if a > R:
                  medium = 0
               elif a < R:
                  medium = np.sqrt(2/(R**2-a**2))
               else:
                  return np.nan
            medium *= quad(InnerIntegral, a, np.inf)[0]
            medium *= -2**1.5*par.NEWTON_G
            return medium
           
        outter  = quad(MediumIntegral, 0, np.inf)[0]
        outter *= Zeta(zp)
        return outter
     
    Pot = quad(OutterIntegral, -np.inf, np.inf)[0]
    return Pot

def Potential_ThickDisk2(R, z, Sigma, Zeta):
    def OutterIntegral(zp):
        def MediumIntegral(a):
            def InnerIntegral(a):
                def fun(Rp):
                    fun  = Rp*Sigma(Rp)
                    fun /= np.sqrt( Rp**2 - a**2 )
                    return fun
                return quad(fun, a, np.inf)[0]
            PlusSqr = np.sqrt( (z-zp)**2 + ( a + R )**2 )
            MinuSqr = np.sqrt( (z-zp)**2 + ( a - R )**2 )
            medium  = derivative(InnerIntegral, a, dx=1e-4)
            if 2*a/( PlusSqr + MinuSqr ) > 1:
               arg = 1
            elif 2*a/( PlusSqr + MinuSqr ) < -1:
               arg = -1
            else:
               arg = 2*a/( PlusSqr + MinuSqr )
            medium *= np.arcsin(arg)
            medium *= 4*par.NEWTON_G
            return medium 
        
        outter = quad(MediumIntegral, 0, np.inf)[0]
        outter *= Zeta(zp)
        return outter
         
    Pot = quad(OutterIntegral, -np.inf, np.inf)[0]
    return Pot


from scipy.special import kn

def Potential_ExponentialDisk(R, z, Sigma, Zeta):
    def OutterIntegral(zp):
        def MediumIntegral(a):
            PlusSqr = np.sqrt( (z-zp)**2 + ( a + R )**2 )
            MinuSqr = np.sqrt( (z-zp)**2 + ( a - R )**2 )
            medium  = a*kn(0, a/par.Disk_Rd)
            if 2*a/( PlusSqr + MinuSqr ) > 1:
               arg = 1
            elif 2*a/( PlusSqr + MinuSqr ) < -1:
               arg = -1
            else:
               arg = 2*a/( PlusSqr + MinuSqr )
            medium *= np.arcsin(arg)
            return medium 
        
        outter = quad(MediumIntegral, 0, np.inf)[0]
        outter *= Zeta(zp)
        return outter
         
    Pot  = quad(OutterIntegral, -np.inf, np.inf)[0]
    Pot *= -4*par.NEWTON_G*par.Disk_Sigma/par.Disk_Rd
    return Pot

def Potential_Disk(R, z):
    def Zeta(z):
        Zeta  = 0.5*par.Disk_alpha0 / par.Disk_z0 * np.exp( -np.absolute(z)/par.Disk_z0 ) 
        Zeta += 0.5*par.Disk_alpha1 / par.Disk_z1 * np.exp( -np.absolute(z)/par.Disk_z1 ) 
        return Zeta
    def Sigma(R):
        return par.Disk_Sigma * np.exp( -R/par.Disk_Rd )
    Pot = Potential_ExponentialDisk(R, z, Sigma, Zeta)
    return Pot


def Potential_ISM(R, z):
    def Zeta(z):
        return 0.5*par.ISM_Sigma/par.ISM_zg * np.exp(- np.absolute(z)/par.ISM_zg)
    def Sigma(R):
        return np.exp( -R/par.ISM_Rg - par.ISM_Rm / R )
    Pot = Potential_ThickDisk2(R, z, Sigma, Zeta)
    if Pot != Pot:
       Pot = Potential_ThickDisk1(R, z, Sigma, Zeta)
    return Pot
   



par.Parameters()
PotDisk     = Potential_Disk(2, 3)
PotISM      = Potential_ISM (2, 3)
