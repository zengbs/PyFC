import numpy as np
import h5py

def _HDF5(self, out):

    # Number of cells along x/y/z
    Nx = self.ni
    Ny = self.nj
    Nz = self.nk

    # Number of cells along each dimension of the input grid.
    ddims = np.array([Nx, Ny, Nz], dtype='int')
    
    # Left edge and right edge coordinates of the desired
    # simulation domain which will be used in GAMER.
    le = np.zeros(3)
    re = np.ones(3)
    
    # Cell spacing
    delta = (re-le)/ddims
    
    # Construct the grid cell edge coordinates
    x = np.linspace(le[0], re[0], ddims[0]+1)
    y = np.linspace(le[1], re[1], ddims[1]+1)
    z = np.linspace(le[2], re[2], ddims[2]+1)
    
    # Find the grid cell midpoints
    x = 0.5*(x[1:]+x[:-1])
    y = 0.5*(y[1:]+y[:-1])
    z = 0.5*(z[1:]+z[:-1])
    print(Nx, Ny, Nz)
    #  Use the 1-D coordinate arrays to consruct 3D coordinate arrays                                                                        
    # that we will use to compute an analytic vector potential
    xx, yy, zz = np.meshgrid(x, y, z, sparse=False, indexing='ij')

    Rho  = out.T
    Ux   = np.full_like( Rho, 0.0 )
    Uy   = np.full_like( Rho, 0.0 )
    Uz   = np.full_like( Rho, 0.0 )
    Temp = np.full_like( Rho, 1   )
    Pres = Temp*Rho

    Dens, MomX, MomY, MomZ, Engy = Pri2Con( Rho, Ux, Uy, Uz, Pres )

    filename = "Hydro_IC"

    f = h5py.File(filename,"w")

    f.create_dataset('Dens', data=Dens)
    f.create_dataset('MomX', data=MomX)
    f.create_dataset('MomY', data=MomY)
    f.create_dataset('MomZ', data=MomZ)
    f.create_dataset('Engy', data=Engy)

    f.create_dataset('x', data=x)
    f.create_dataset('y', data=y)
    f.create_dataset('z', data=z)

    f.flush()
    f.close()

def Pri2Con( Rho, Ux, Uy, Uz, Pres ):

    LorentzFactor = np.sqrt( 1.0 + np.square(Ux) + np.square(Uy) + np.square(Uz) );
    Temperature = Pres/Rho;
    HTilde = EoS_Temp2HTilde( Temperature );
    Dens = Rho*LorentzFactor;
    Factor0 = Dens*HTilde + Dens;
  
    MomX = Ux*Factor0;
    MomY = Uy*Factor0;
    MomZ = Uz*Factor0;
    MSqr_DSqr  = np.square(MomX)+np.square(MomY)+np.square(MomZ);
    MSqr_DSqr /= np.square(Dens);
    HTildeFunction = SRHD_HTildeFunction( HTilde, MSqr_DSqr, Temperature );
    Engy  = MSqr_DSqr + HTildeFunction;
    Engy /= 1.0 + np.sqrt( 1.0 + MSqr_DSqr + HTildeFunction );
    Engy *= Dens;

    return Dens, MomX, MomY, MomZ, Engy

def EoS_Temp2HTilde( Temp ):
    TempSqr = Temp*Temp;
    Factor = 2.25*TempSqr;
  
    HTilde  = Factor;
    HTilde /= 1.0 + np.sqrt( Factor + 1.0 );
    HTilde += 2.5*Temp;
  
    return HTilde;

def SRHD_HTildeFunction( HTilde, MSqr_DSqr, Temp ):

    H =  HTilde + 1.0;
    Factor0 = np.square( H ) + MSqr_DSqr;
                                            
    Constant = 0.0                                
    HTildeFunction = np.square( HTilde ) + 2.0*HTilde - 2.0*Temp - 2.0*Temp*HTilde
    + np.square( Temp * H ) / Factor0 + Constant;

    return HTildeFunction 
