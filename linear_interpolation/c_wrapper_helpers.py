import ctypes
import numpy as np
import numpy.ctypeslib as npct

DOUBLE = ctypes.c_double
DOUBLE_P = ctypes.POINTER(ctypes.c_double)
DOUBLE_PP = ctypes.POINTER(DOUBLE_P)

def convert_to_float_arr_2D(ptr, n, m):
    arr = np.zeros(shape=(n, m))
    
    for i in range(n):
        for j in range(m):
            arr[i, j] = ptr[i][j]
            
    return arr

def convert_to_ptr_arr_2D(arr):
    n = np.size(arr, axis = 0)
    m = np.size(arr, axis = 1)
    
    dim_X = n * DOUBLE_P
    dim_Y = m * DOUBLE

    arr_ptr = dim_X()

    for i, row in enumerate(arr):
        arr_ptr[i] = dim_Y()

        for j, val in enumerate(row):
            arr_ptr[i][j] = val
            
    return arr_ptr

def convert_to_ptr_arr_1D(arr):            
    return npct.as_ctypes(arr)