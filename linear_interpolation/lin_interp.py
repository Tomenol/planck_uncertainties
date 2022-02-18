# -*- coding: utf-8 -*-
"""
Created on Wed Jun  9 13:53:42 2021

@author: Thomas Maynadié
"""

import ctypes
import numpy.ctypeslib as npct
import numpy as np

import linear_interpolation.c_wrapper_helpers as c_wrapper_helpers

lib = ctypes.cdll.LoadLibrary("./linear_interpolation/src/lin_interp.so")


def interp2D(data, grid, x, y):    
    n = np.size(data, axis=0)
    m = np.size(data, axis=1)
    
    lib.interp2D.restype = ctypes.c_double
    lib.interp2D.argtypes = [ctypes.c_double,  
                             ctypes.c_double, 
                             np.ctypeslib.ndpointer(dtype=ctypes.c_double, ndim=data.ndim, shape=data.shape), 
                             np.ctypeslib.ndpointer(dtype=ctypes.c_double, ndim=grid.ndim, shape=grid.shape), 
                             ctypes.c_int,  
                             ctypes.c_int]
    
    res = lib.interp2D(ctypes.c_double(x),
                       ctypes.c_double(y),
                       data,
                       grid,
                       ctypes.c_int(n),
                       ctypes.c_int(m))
    
    return res

def interpArr2D(new_data, new_grid, old_data, old_grid):
    new_data.astype(np.float64)
    old_data.astype(np.float64)
    new_grid.astype(np.float64)
    old_grid.astype(np.float64)
    
    new_n = np.size(new_data, axis=0)
    old_n = np.size(old_data, axis=0)
    
    new_m = np.size(new_data, axis=1)
    old_m = np.size(old_data, axis=1)
    
    lib.interp_array_2D.argtypes = [np.ctypeslib.ndpointer(dtype=ctypes.c_double, ndim=new_data.ndim, shape=new_data.shape), 
                                   np.ctypeslib.ndpointer(dtype=ctypes.c_double, ndim=new_grid.ndim, shape=new_grid.shape), 
                                   np.ctypeslib.ndpointer(dtype=ctypes.c_double, ndim=old_data.ndim, shape=old_data.shape), 
                                   np.ctypeslib.ndpointer(dtype=ctypes.c_double, ndim=old_grid.ndim, shape=old_grid.shape), 
                                   ctypes.c_int, 
                                   ctypes.c_int, 
                                   ctypes.c_int, 
                                   ctypes.c_int]
    
    lib.interp_array_2D(new_data, new_grid, old_data, old_grid, ctypes.c_int(new_n), ctypes.c_int(new_m), ctypes.c_int(old_n), ctypes.c_int(old_m))