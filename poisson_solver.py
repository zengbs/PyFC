import numpy as np


# Poisson solver with periodic boundary condition

def PoissonSolver(source, Lx, Ly, Lz, Nx, Ny, Nz):
#def PoissonSolver(source, Lx, Ly, Lz, Nx, Ny, Nz):
    I           = np.indices((Nx, Ny, Nz))[0]/Lx
    J           = np.indices((Nx, Ny, Nz))[1]/Ly
    K           = np.indices((Nx, Ny, Nz))[2]/Lz

    # change 0 to -1 to avoid division by zero, 
    I           = np.where(I==0, -1, I)
    J           = np.where(J==0, -1, J)
    K           = np.where(K==0, -1, K)

    answer_FFT  = np.fft.fftn(source)

    # do not take into account the zero terms, which lead to division by zero
    answer_FFT  = np.where(np.logical_and(I==-1,J==-1,K==-1), 0, answer_FFT/(I**2 + J**2 + K**2))

    answer_FFT /= 4*np.pi**2

    answer      = np.fft.ifftn(answer_FFT)

    return answer

# exact solution:

def ExactSolution(Lx, Ly, Lz, Nx, Ny, Nz):

    delta_x     = Lx/Nx
    delta_y     = Ly/Ny
    delta_z     = Lz/Nz

    x           = np.indices((int(Nx), int(Ny), int(Nz)))[0]*delta_x
    y           = np.indices((int(Nx), int(Ny), int(Nz)))[1]*delta_y
    z           = np.indices((int(Nx), int(Ny), int(Nz)))[2]*delta_z

    source = 2*np.exp(x+y+z)*( z*y*(1-y)*(1-z) + x*z*(1-x)*(1-z) + x*y*(1-x)*(1-y) )
  
    potential = np.exp(-(x+y+z))*x*y*z*(1-x)*(1-y)*(1-z)

    return potential, source


if __name__ == '__main__':
    Lx = 1.0
    Ly = 1.0
    Lz = 1.0

    Nx = 256
    Ny = 256
    Nz = 256

    AnalyticalPotential, Source = ExactSolution(Lx, Ly, Lz, Nx, Ny, Nz)

    NumericalPotential = PoissonSolver(Source, Lx, Ly, Lz, Nx, Ny, Nz)

    print(np.max(1-AnalyticalPotential/NumericalPotential))      
