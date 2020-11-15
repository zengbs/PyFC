import numpy as np

def _DumpFile(out):

    # Number of cells along x/y/z
    Nx = out.shape[0] 
    Ny = out.shape[1] 
    Nz = out.shape[2] 


    Rho  = out
    Ux   = np.full_like( Rho, 0.0 )
    Uy   = np.full_like( Rho, 0.0 )
    Uz   = np.full_like( Rho, 0.0 )
    Temp = np.full_like( Rho, 1   )
    Pres = Temp*Rho

    Dens, MomX, MomY, MomZ, Engy = Pri2Con( Rho, Ux, Uy, Uz, Pres )
  
    filename = "UM_IC"

    FluidVar = np.zeros((5, Nx, Ny, Nz), dtype=np.float64)

    FluidVar[0] = Rho
    FluidVar[1] = Ux
    FluidVar[2] = Uy
    FluidVar[3] = Uz
    FluidVar[4] = Temp

    FluidVar.tofile(filename)


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
