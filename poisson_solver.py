import numpy as np



def PoissonSolver(source, Lx, Ly, Lz, Nx, Ny, Nz):

    return potential

# exact solution:
def ExactSolution(Lx, Ly, Lz, Nx, Ny, Nz):

    delta_x     = Lx/Nx 
    delta_y     = Ly/Ny
    delta_z     = Lz/Nz

    x           = (np.indices((int(Nx), int(Ny), int(Nz)))[0] - Nx/2 + 0.5)*delta_x
    y           = (np.indices((int(Nx), int(Ny), int(Nz)))[1] - Ny/2 + 0.5)*delta_y
    z           = (np.indices((int(Nx), int(Ny), int(Nz)))[2] - Nz/2 + 0.5)*delta_z

    r           = np.sqrt(x**2+y**2+z**2)

    source      = 1
    source     /= r * np.power( 1 + r, 3 )

    potential   = -1/(2*(1+r))

    return potential, source


if __name__ == '__main__':
    Lx = 1.0
    Ly = 1.0
    Lz = 1.0

    Nx = 128 
    Ny = 128 
    Nz = 128 

    AnalyticalPotential, Source = ExactSolution(Lx, Ly, Lz, Nx, Ny, Nz)

    NumericalPotential = PoissonSolver(Source, Lx, Ly, Lz, Nx, Ny, Nz)

 

    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm

    fig1 = plt.figure()
    pos = plt.imshow(np.real(AnalyticalPotential[:,:,int(Nz/2)]), norm=None, cmap='nipy_spectral')
    fig1.colorbar(pos)
    fig1.savefig('analytical.png')
     
    fig2 = plt.figure()
    pos = plt.imshow(np.real(NumericalPotential[:,:,int(Ny/2)]), norm=None, cmap='nipy_spectral')
    fig2.colorbar(pos)
    fig2.savefig('numerical.png')
     
     
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
