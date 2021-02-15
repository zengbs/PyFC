import numpy as np

UM_IC_double = np.fromfile("UM_IC_double", dtype='float64')
UM_IC_float = UM_IC_double.astype(dtype='float32')
UM_IC_float.tofile("UM_IC_float")

ExtPotTable_double = np.fromfile("ExtPotTable_double", dtype='float64')
ExtPotTable_float = ExtPotTable_double.astype(dtype='float32')
ExtPotTable_float.tofile("ExtPotTable_float")
