import numpy as np

UM_IC_double = np.fromfile("UM_IC_double", dtype='float64')
UM_IC_float = UM_IC_double.astype(dtype='float32')
UM_IC_float.tofile("UM_IC_float")

ExtPotTable_double = np.fromfile("ExtPotTable_double", dtype='float64')
ExtPotTable_float = ExtPotTable_double.astype(dtype='float32')
ExtPotTable_float.tofile("ExtPotTable_float")

FractalDensity_double = np.fromfile("FractalDensity_double", dtype='float64')
FractalDensity_float = FractalDensity_double.astype(dtype='float32')
FractalDensity_float.tofile("FractalDensity_float")

FractalUx_double = np.fromfile("FractalUx_double", dtype='float64')
FractalUx_float = FractalUx_double.astype(dtype='float32')
FractalUx_float.tofile("FractalUx_float")

FractalUy_double = np.fromfile("FractalUy_double", dtype='float64')
FractalUy_float = FractalUy_double.astype(dtype='float32')
FractalUy_float.tofile("FractalUy_float")

FractalUz_double = np.fromfile("FractalUz_double", dtype='float64')
FractalUz_float = FractalUz_double.astype(dtype='float32')
FractalUz_float.tofile("FractalUz_float")
