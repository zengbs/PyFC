import sys
import numpy as np
import parameters as par

# **********************************************************************
# description: create four 3D array storing 
#              (1) x-coordinate
#              (2) y-coordinate
#              (3) z-coordinate
#              (4) the distance between cell itself and center
#              (5) 
# input      : None
# output     : 3D array stroing x, y, z, and radial distance
# **********************************************************************
def Create3DCoordinateArray(Nx, Ny, Nz):
    Idx       = np.indices((Nx, Ny, Nz))[0]
    Jdx       = np.indices((Nx, Ny, Nz))[1]
    Kdx       = np.indices((Nx, Ny, Nz))[2]
              
    X         = Idx+0.5-par.Nx*0.5
    Y         = Jdx+0.5-par.Ny*0.5
    Z         = Kdx+0.5-par.Nz*0.5
    X        *= par.delta[0]
    Y        *= par.delta[1]
    Z        *= par.delta[2]
    r         = np.sqrt(X**2+Y**2+Z**2)
    R         = np.sqrt(X**2+Y**2)
    return X, Y, Z, r, R

#par.Parameters()
#X, Y, Z, r, R   = Create3DCoordinateArray(3, 3, 3)
#Idx = np.indices((3,3,3), dtype=par.Precision)[0]
#Jdx = np.indices((3,3,3), dtype=par.Precision)[1]
#Kdx = np.indices((3,3,3), dtype=par.Precision)[2]
#print("Idx------------------------------")
#for i in range(3):
# for j in range(3):
#  for k in range(3):
#    print("X[%d][%d][%d] = %f" % ( i, j, k, X[i][j][k]))
#print("Jdx------------------------------")
#for i in range(3):
# for j in range(3):
#  for k in range(3):
#    print("Y[%d][%d][%d] = %f" % ( i, j, k, Y[i][j][k]))
#print("Kdx------------------------------")
#for i in range(3):
# for j in range(3):
#  for k in range(3):
#    print("Z[%d][%d][%d] = %f" % ( i, j, k, Z[i][j][k]))
#exit(0)

def Bulge(mSqr):
    Rho  = par.Bulge_Rho0
    Rho *= np.power( np.sqrt(mSqr)/par.Bulge_a, -par.Bulge_alpha )
    Rho *= np.exp( -mSqr / par.Bulge_r**2 )
    return Rho

def DarkHalo(mSqr):
    Rho  = par.Halo_Rho0
    Rho *= np.power( mSqr / par.Halo_a, -par.Halo_alpha )
    Rho *= np.power( 1 + mSqr/par.Halo_a, par.Halo_alpha - par.Halo_beta )
    return Rho

def StellarDisk(R, z):
    Rho   = 0.5*par.Disk_alpha0 / par.Disk_z0 * np.exp( -np.absolute(z)/par.Disk_z0 ) 
    Rho  += 0.5*par.Disk_alpha1 / par.Disk_z1 * np.exp( -np.absolute(z)/par.Disk_z1 ) 
    Rho  *= par.Disk_Sigma * np.exp( -R/par.Disk_Rd )
    return Rho

def ISM(R, z):
    Rho  = 0.5*par.ISM_Sigma/par.ISM_zg
    Rho *= np.exp( -R/par.ISM_Rg - par.ISM_Rm / R - np.absolute(z)/par.ISM_zg )
    #Array=np.where(Rho**2>0.0, False, True)
    #print(np.any(Array))
    #exit(0)
    return Rho

def Miyamoto(R, z):
    Temp = np.sqrt(z**2 + par.Miyamoto_b**2)
    Rho  = par.Miyamoto_b**2 * par.Miyamoto_M
    Rho /= 4 * np.pi
    Rho *= par.Miyamoto_a*R**2 + ( par.Miyamoto_a + 3*Temp )*( par.Miyamoto_a + Temp )**2
    Rho /= np.power( R**2 + ( par.Miyamoto_a + Temp )**2, 2.5 )
    Rho /= np.power( Temp , 3)
    return Rho

def TotGasDensity():
    x, y, z, r, R = Create3DCoordinateArray(par.Nx, par.Ny, par.Nz )
    #mSqrBulge = R**2 + z**2/par.Bulge_q
    #mSqrHalo  = R**2 + z**2/par.Halo_q
    #TotalDensity  = Bulge(mSqrBulge) 
    #TotalDensity += DarkHalo(mSqrHalo) 
    #TotalDensity += StellarDisk(R, z) 
    #TotalDensity += ISM(R, z)
    TotalDensity = Miyamoto(R, z)
    return TotalDensity
