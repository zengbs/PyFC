import numpy as np
from scipy.integrate import quad
from scipy.interpolate import interp2d
from scipy.misc import derivative
import parameters as par
import density_profile
from scipy.special import kn
from density_profile import Create3DCoordinateArray, Bulge, DarkHalo
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import sys

par.Parameters()

def Check(x):
    if not np.isfinite(x):
       return True
    if np.isnan(x):
       return True


# (2.125b)
def Potential_Spheroidal (R, z, DensFun, e):
    a0 = 1 # arbitrary constant
    def phi (mSqr, dens_func): # (2.117)
        phi = quad(dens_func, 0, mSqr)[0]
        return phi

    def SecondTerm( tau ): # the last term in (2.125b)
        c0 = np.sqrt( 1 - e**2 )*a0
         
        # m^2 = (2.125a)
        mSqr  = R**2 / ( tau + a0**2 )
        mSqr += z**2 / ( tau + c0**2 )
        mSqr *= a0**2
        integrand  = phi(mSqr, DensFun)
        integrand /= tau + a0**2
        integrand /= np.sqrt(tau + c0**2)
        return integrand

    Pot  =  -2.0*np.pi*par.NEWTON_G * np.sqrt( 1 - e**2 ) / e
    Pot1 = phi(np.inf, DensFun)*np.arcsin(e)
    Pot2 = 0.5*a0*e*quad(SecondTerm, 0, np.inf)[0]
    Pot *= Pot1 - Pot2
    return Pot


def Potential_Bulge(R, z):
    Pot = Potential_Spheroidal(R, z, Bulge, np.sqrt(1-par.Bulge_q**2))
    return Pot

def Potential_DarkHalo(R, z):
    Pot = Potential_Spheroidal(R, z, DarkHalo, np.sqrt(1-par.Halo_q**2))
    return Pot

# (2.154) and (2.169)
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

# (2.153a) and (2.169)
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




# (2.170)
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

# (2.210)
def Potential_Disk(R, z):
    def Zeta(z):
        Zeta  = 0.5*par.Disk_alpha0 / par.Disk_z0 * np.exp( -np.absolute(z)/par.Disk_z0 ) 
        Zeta += 0.5*par.Disk_alpha1 / par.Disk_z1 * np.exp( -np.absolute(z)/par.Disk_z1 ) 
        return Zeta
    def Sigma(R):
        return par.Disk_Sigma * np.exp( -R/par.Disk_Rd )
    Pot = Potential_ExponentialDisk(R, z, Sigma, Zeta)
    return Pot

# (2.211)
def Potential_ISM(R, z):
    def Zeta(z):
        return 0.5*par.ISM_Sigma/par.ISM_zg * np.exp(- np.absolute(z)/par.ISM_zg)
    def Sigma(R):
        return np.exp( -R/par.ISM_Rg - par.ISM_Rm / R )
    Pot = Potential_ThickDisk2(R, z, Sigma, Zeta)
    if Pot != Pot:
       Pot = Potential_ThickDisk1(R, z, Sigma, Zeta)
    return Pot

# (2.69a)
def Potential_Miyamoto(R, z):
    Temp = np.sqrt(z**2 + par.Miyamoto_b**2)
    Pot = np.power(R**2 + ( par.Miyamoto_a + Temp )**2, -0.5)
    return -Pot*par.NEWTON_G*par.Miyamoto_M


def TotPotential():
    Nx = par.Nx + 2*par.GRA_GHOST_SIZE
    Ny = par.Ny + 2*par.GRA_GHOST_SIZE
    Nz = par.Nz + 2*par.GRA_GHOST_SIZE

    Idx_Extended = np.indices((int(1.5*Nx/2),int(1.5*Ny/2),int(1.5*Nz/2)))[0]
    Idx          = np.indices((int(    Nx/2),int(    Ny/2),int(    Nz/2)))[0]
    Jdx          = np.indices((int(    Nx/2),int(    Ny/2),int(    Nz/2)))[1]
    Kdx          = np.indices((int(    Nx/2),int(    Ny/2),int(    Nz/2)))[2]

    X3D_Extended = (Idx_Extended+0.5)*par.delta[0]
    X3D          = (Idx         +0.5)*par.delta[0]
    Y3D          = (Jdx         +0.5)*par.delta[1]
    Z3D          = (Kdx         +0.5)*par.delta[2]

    X1D_Extended = X3D_Extended[:,0,0]
    X1D          = X3D         [:,0,0]
    Y1D          = Y3D         [0,:,0]
    Z1D          = Z3D         [0,0,:]

    IJdxSqr = (Idx**2 + Jdx**2)[:,:,0]

    IJdxSqr_1D, index_indices, inverse_indices = np.unique(IJdxSqr, return_index=True, return_inverse=True) 
   
    R1D = np.sqrt(Y3D**2+Z3D**2)
    R1D = R1D.flatten()[index_indices]
    Pot2D_Extended = np.zeros((len(X1D_Extended), len(Z1D)))

    # potential calculation
    for i in range(len(X1D_Extended)):
       for k in range(len(Z1D)):
           Pot2D_Extended[i][k]  = Potential_Miyamoto ( X1D_Extended[i], Z1D[k] )
           #Pot2D_Extended[i][k]  = Potential_Bulge   ( X1D_Extended[i], Z1D[k] )
           #Pot2D_Extended[i][k] += Potential_DarkHalo( X1D_Extended[i], Z1D[k] )
           #Pot2D_Extended[i][k] += Potential_Disk    ( X1D_Extended[i], Z1D[k] )
           #Pot2D_Extended[i][k] += Potential_ISM     ( X1D_Extended[i], Z1D[k] )
       print("i=%03d/%03d" % (i, len(X1D_Extended)))
       sys.stdout.flush()
     
    # interpolation
    Pot2DFun = interp2d( Z1D, X1D_Extended, Pot2D_Extended, kind='cubic' )
    Pot2D = Pot2DFun(Z1D, R1D)
    # Pot2D.shape = (128,5839)

    # convert 2D -> 3D
    Pot3D = np.zeros((int(Nx/2),int(Ny/2),int(Nz/2)))

    for k in range(0,len(Z1D)):
      Pot3D[k,:,:] = Pot2D[:,k][inverse_indices].reshape(int(Nx/2),int(Ny/2))

    

    # flip 
    Pot3D = np.concatenate((np.flip(Pot3D, axis=2), Pot3D),axis=2)
    Pot3D = np.concatenate((np.flip(Pot3D, axis=0), Pot3D),axis=0)
    Pot3D_up = np.concatenate((np.flip(Pot3D, axis=1), Pot3D),axis=1)

    #############################################################
    if __name__ == '__main__':

      Idx = np.indices((int(Nx),int(Ny),int(Nz)))[1]
      Jdx = np.indices((int(Nx),int(Ny),int(Nz)))[2]                                                                         
      Kdx = np.indices((int(Nx),int(Ny),int(Nz)))[0]

      X3D = (Idx + 0.5)*par.delta[0]
      Y3D = (Jdx + 0.5)*par.delta[1]
      Z3D = (Kdx + 0.5)*par.delta[2]
   
      CenterX = (X3D[0][0][0]+X3D[ 0][-1][ 0])*0.5
      CenterY = (Y3D[0][0][0]+Y3D[ 0][ 0][-1])*0.5
      CenterZ = (Z3D[0][0][0]+Z3D[-1][ 0][ 0])*0.5

      X3D -= CenterX
      Y3D -= CenterY
      Z3D -= CenterZ
 
      R3D = np.sqrt(X3D**2+Y3D**2)

      Pot3D_down = Potential_Miyamoto( R3D, Z3D )

       
      return Pot3D_up, Pot3D_down
    else:
      return Pot3D_up

if __name__ == '__main__':

  fig1 = plt.figure()
  
  Pot3D_up, Pot3D_down = TotPotential()
  
  Nx = par.Nx + par.GRA_GHOST_SIZE
  Ny = par.Ny + par.GRA_GHOST_SIZE
  Nz = par.Nz + par.GRA_GHOST_SIZE
  
  pos = plt.imshow(Pot3D_up[:,:,int(Nz/3)], norm=LogNorm(), cmap='nipy_spectral')
  fig1.colorbar(pos)
  fig1.savefig('image_up.png')
  
  fig2 = plt.figure()
  pos = plt.imshow(Pot3D_down[:,:,int(Nz/3)], norm=LogNorm(), cmap='nipy_spectral')
  fig2.colorbar(pos)
  fig2.savefig('image_down.png')
  
  
  #// profile
  f, ax = plt.subplots(1,1)
   
  f.subplots_adjust( hspace=0.05, wspace=0.35 )
  f.set_size_inches( 7, 5 ) 
   
  ax.plot(Pot3D_up  [:,int(Ny/3),int(Nz/3)])
  ax.plot(Pot3D_down[:,int(Ny/3),int(Nz/3)])
  
  print(np.amax(np.absolute(1-np.divide(Pot3D_up,Pot3D_down)),axis=2))
  
  f.savefig('profile.png')
