import numpy as np


def PoissonSolver(source):
    I           = np.indices(Nx,Ny,Nz)[0]*2*np.pi/par.Lx
    J           = np.indices(Nx,Ny,Nz)[1]*2*np.pi/par.Ly
    K           = np.indices(Nx,Ny,Nz)[2]*2*np.pi/par.Lz

    source_FFT  = np.fft.rfftn(source)
    answer_FFT  = -source_FFT
    answer_FFT /= (I**2 + J**2 + K**2)

    answer      = np.fft.ifftn(answer_FFT)

    return answer
