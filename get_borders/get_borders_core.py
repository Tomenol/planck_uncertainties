# -*- coding: utf-8 -*-
"""
Created on Thu Jun 10 19:43:12 2021

@author: Thomas Maynadié
"""

import ctypes
import numpy.ctypeslib as npct
import numpy as np

import linear_interpolation.c_wrapper_helpers as c_wrapper_helpers

lib = ctypes.cdll.LoadLibrary("./get_borders/src/get_borders.so")

def get_contours(data, grid, level, normalize=False):   
    data.astype(np.float64)
    grid.astype(np.float64)
    
    n = np.size(data, axis=0)
    m = np.size(data, axis=1)
            
    lib.get_contour.argtypes = [np.ctypeslib.ndpointer(dtype=ctypes.c_double, ndim=data.ndim, shape=data.shape), 
                                   np.ctypeslib.ndpointer(dtype=ctypes.c_double, ndim=grid.ndim, shape=grid.shape), 
                                   ctypes.c_double,
                                   ctypes.c_int, 
                                   ctypes.c_int]
    
    lib.get_contour(data, grid, ctypes.c_double(level), ctypes.c_int(n), ctypes.c_int(m))
    
    
    lib.restype = ctypes.c_int
    
    size = int(lib.get_contour_size())
    
    x = np.ndarray(size, dtype=np.float64)    
    y = np.ndarray(size, dtype=np.float64)  
    
    x_i = np.ndarray(size, dtype=np.int32)    
    y_i = np.ndarray(size, dtype=np.int32)  
    
    lib.get_contour_values.argtypes = [np.ctypeslib.ndpointer(dtype=ctypes.c_double, ndim=x.ndim, shape=x.shape), 
                                   np.ctypeslib.ndpointer(dtype=ctypes.c_double, ndim=y.ndim, shape=y.shape)]
    
    lib.get_contour_indices.argtypes = [np.ctypeslib.ndpointer(dtype=ctypes.c_int, ndim=x_i.ndim, shape=x_i.shape), 
                                   np.ctypeslib.ndpointer(dtype=ctypes.c_int, ndim=y_i.ndim, shape=y_i.shape)]
    
    
    lib.get_contour_values(x, y)
    lib.get_contour_indices(x_i, y_i)
    
    return x, y, x_i, y_i

def normalize_pdf(data, level):
    data.astype(np.float64)
    
    n = np.size(data, axis=0)
    m = np.size(data, axis=1)
    
    lib.get_norm_factor.argtypes = [np.ctypeslib.ndpointer(dtype=ctypes.c_double, ndim=data.ndim, shape=data.shape), 
                                   ctypes.c_int, 
                                   ctypes.c_int]
    
    lib.get_norm_factor.restype = ctypes.c_double
    
    normalization_factor = lib.get_norm_factor(data, ctypes.c_int(n), ctypes.c_int(m))
    
    data = data * normalization_factor
    
    return data, level * normalization_factor