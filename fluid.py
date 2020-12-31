import parameters as par
import numpy as np

# *************************************************
# description: temperature to sound speed
# input      : temperature (K)
# output     : sound speed (cm/s)
# **************************************************


def Tem2Cs( Temperature ):
    eta  = par.Const_kB*Temperature
    eta /= par.Const_AtomicMassUnit*par.Const_MolecularWeight
    eta /= par.Const_C**2
    eta /= par.Const_Erg2eV

    h = 2.5*eta + np.sqrt( 2.25 * eta**2 + 1.0 )                                                                                          
    hTilde = 2.5*eta + 2.25*eta*eta/(1+np.sqrt(1+2.25*eta**2))
    a = eta *( 5.0*h - 8.0*eta )
    b = 3.0 * h * (h - eta)
    Cs_sq = a / b
    Cs = np.sqrt( Cs_sq )
    return Cs*par.Const_C


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
