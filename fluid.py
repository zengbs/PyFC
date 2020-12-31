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
