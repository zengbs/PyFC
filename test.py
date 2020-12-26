import numpy as np
from multiprocessing import Pool

def fill_array(start_val):
    return list(range(start_val, start_val+10))

if __name__=='__main__':
    pool = Pool(processes=4)
    list_start_vals = range(40, 60)
    array_2D = np.array(pool.map(fill_array, list_start_vals))
    print(pool.map(fill_array, list_start_vals))
    pool.close() # ATTENTION HERE
    #print array_2D
