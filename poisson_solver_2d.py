import numpy as np


# Poisson solver with periodic boundary condition

def PoissonSolver(source, Lx, Ly, Nx, Ny):
    Nx *= 1
    Ny *= 1 
    i           = np.indices((int(Nx), int(Ny)))[0]
    j           = np.indices((int(Nx), int(Ny)))[1]

    source_FFT  = np.fft.fft2(source)
    print(source_FFT.shape)
   
    DC               = source_FFT[0][0] 
    source_FFT      *= (Lx/Nx)*(Ly/Ny)
    source_FFT      /= 2*( np.cos(2*np.pi*i/Nx) + np.cos(2*np.pi*j/Ny) - 2 )
    source_FFT[0][0] = DC

    potential   = np.fft.ifft2(source_FFT)

    return potential

# exact solution:

def ExactSolution(Lx, Ly, Nx, Ny):
    pi = np.pi


    delta_x     = Lx/Nx
    delta_y     = Ly/Ny

    x           = ((np.indices((int(Nx), int(Ny)))[0])+0.5)*delta_x
    y           = ((np.indices((int(Nx), int(Ny)))[1])+0.5)*delta_y

    # Case 1
    source    =  2*pi**2*np.sin(pi*x)*np.sin(pi*y)
    potential =  np.sin(pi*x)*np.cos(pi*y+pi/2)

    # Case 2
    #source    =  5*pi**2*np.sin(pi*y)*np.sin(2*pi*x)
    #potential =  np.sin(pi*y)*np.cos(2*pi*x+pi/2)

    # Case 3
    #source    =  5*pi**2*np.sin(pi*x)*np.sin(2*pi*y)
    #potential =  np.sin(pi*x)*np.cos(2*pi*y+pi/2)

    # Case 4
    #source    =  8*pi**2*np.sin(2*pi*y)*np.sin(2*pi*x)
    #potential =  np.sin(2*pi*y)*np.cos(2*pi*x+pi/2)

    # Case 5
    #source    = 8*pi**2*( np.sin(2*pi*y)**2*np.cos(4*pi*x) + np.cos(4*pi*y)*np.sin(2*pi*x)**2 ) 
    #potential =  np.sin(2*pi*y)**2*np.cos(2*pi*x+pi/2)**2

    # Case 6
    #source     = y*(1-y)*(-x**2+5*x-4)*np.exp(-(x+y)) + x*(1-x)*(-y**2+5*y-4)*np.exp(-(x+y))
    #potential  = x*y*(1-x)*(1-y)*np.exp(-(x+y))

    # Case 7
    #source     = -2*y*(1-y)-2*x*(1-x)
    #potential  = x*y*(1-x)*(1-y)

    # Case 8
    #source     = 2*pi**2*( np.sin(pi*y)**2*np.cos(2*pi*x) + np.cos(2*pi*y)*np.sin(pi*x)**2 ) 
    #potential  =  np.sin(pi*y)**2*np.cos(pi*x+pi/2)**2

    #source = np.concatenate((source,source), axis=0)
    #source = np.concatenate((source,source), axis=1)
    #source = np.concatenate((source,source), axis=0)
    #source = np.concatenate((source,source), axis=1)

    #potential = np.concatenate((potential,potential), axis=0)
    #potential = np.concatenate((potential,potential), axis=1)
    #potential = np.concatenate((potential,potential), axis=0)
    #potential = np.concatenate((potential,potential), axis=1)

    return potential, source


if __name__ == '__main__':
    Lx = 1.0
    Ly = 1.0

    Nx = 256 
    Ny = 256 


    AnalyticalPotential, Source = ExactSolution(Lx, Ly, Nx, Ny)

    NumericalPotential = PoissonSolver(Source, Lx, Ly, Nx, Ny)

 

    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm


    f, ax = plt.subplots(1,3)
 
    f.subplots_adjust( hspace=0.05, wspace=0.35 )
    f.set_size_inches( 10, 5 ) 


    pos1 = ax[0].imshow(Source, norm=None, cmap='nipy_spectral')
    pos1 = ax[1].imshow(np.real(AnalyticalPotential), norm=None, cmap='nipy_spectral')
    pos2 = ax[2].imshow(np.real(NumericalPotential),  norm=None, cmap='nipy_spectral')
    ax[0].tick_params ( axis='both', which='both', colors='k', direction='in', top=True  , bottom=True, labelsize=14 ) 
    ax[1].tick_params ( axis='both', which='both', colors='k', direction='in', top=True  , bottom=True, labelsize=14 ) 
    ax[2].tick_params ( axis='both', which='both', colors='k', direction='in', top=True  , bottom=True, labelsize=14 ) 

    ax[0].set_xticklabels([])
    ax[1].set_xticklabels([])
    ax[2].set_xticklabels([])

    ax[0].set_yticklabels([])
    ax[1].set_yticklabels([])
    ax[2].set_yticklabels([])
    f.savefig('compare.png', bbox_inches='tight')
    plt.show()
     
     
    ##// profile
    #f, ax = plt.subplots(1,1)
    # 
    #f.subplots_adjust( hspace=0.05, wspace=0.35 )
    #f.set_size_inches( 7, 5 ) 
    # 
    #ax.plot(Pot3D_up  [:,int(Ny/3),int(Nz/3)])
    #ax.plot(Pot3D_down[:,int(Ny/3),int(Nz/3)])
    # 
    #print(np.amax(np.absolute(1-np.divide(Pot3D_up,Pot3D_down)),axis=2))
    # 
    #f.savefig('profile.png')
