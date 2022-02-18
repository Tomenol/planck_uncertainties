# -*- coding: utf-8 -*-
"""
Created on Fri Jun 11 10:41:52 2021

@author: Thomas Maynadié
"""

import ctypes
import numpy.ctypeslib as npct
import numpy as np

import linear_interpolation.c_wrapper_helpers as c_wrapper_helpers

lib = ctypes.cdll.LoadLibrary("./optimizer/src/optimizerv2.so")

def optimize_evs(previous_values, X_names, pdf, grid, contours, coefficients, initstep=False, maximize=True, get_history=False):       
    pdf.astype(np.float64)
    grid.astype(np.float64)
    
    n = np.size(pdf, axis=0)
    m = np.size(pdf, axis=1)
    contour_arr_size = len(contours[0])
    
    X_0 = np.copy(previous_values)
    
    X_names = np.char.encode(X_names, encoding='utf-8')
    X_names_ptr = (ctypes.c_char_p * len(X_names))()
    X_names_ptr[:] = X_names
    
    DOUBLE = ctypes.c_double
    DOUBLE_P = ctypes.POINTER(DOUBLE)
    DOUBLE_PP = ctypes.POINTER(DOUBLE_P)
    
    INT = ctypes.c_int
    INT_P = ctypes.POINTER(INT)

    X = DOUBLE_P()
    Y = DOUBLE_P()
    size = INT()
            
    lib.parameters_optimisation.argtypes = [np.ctypeslib.ndpointer(dtype=DOUBLE, ndim=X_0.ndim, shape=X_0.shape),
                                            ctypes.POINTER(ctypes.c_char_p),
                                            np.ctypeslib.ndpointer(dtype=DOUBLE, ndim=pdf.ndim, shape=pdf.shape), 
                                            np.ctypeslib.ndpointer(dtype=DOUBLE, ndim=grid.ndim, shape=grid.shape), 
                                            np.ctypeslib.ndpointer(dtype=DOUBLE, ndim=contours.ndim, shape=contours.shape), 
                                            np.ctypeslib.ndpointer(dtype=DOUBLE, ndim=coefficients.ndim, shape=coefficients.shape), 
                                            INT, 
                                            INT, 
                                            INT, 
                                            ctypes.c_bool, 
                                            ctypes.c_bool, 
                                            DOUBLE_PP, 
                                            DOUBLE_PP, 
                                            INT_P]
        
    lib.parameters_optimisation(X_0, 
                                X_names_ptr, 
                                pdf,
                                grid, 
                                contours, 
                                coefficients, 
                                ctypes.c_int(n), 
                                ctypes.c_int(m), 
                                ctypes.c_int(contour_arr_size), 
                                ctypes.c_bool(maximize), 
                                ctypes.c_bool(initstep),
                                ctypes.byref(X),
                                ctypes.byref(Y),
                                ctypes.byref(size))
                    
    if get_history is True:
        X_ret = np.ndarray([2, int(size.value)], dtype=np.float64)
        
        for i in range(int(size.value)):
            X_ret[0][i] = X[i]
            X_ret[1][i] = Y[i]
            
    else:
        X_ret = np.ndarray(2, dtype=np.float64)
        
        X_ret[0] = X[int(size.value) - 1]
        X_ret[1] = Y[int(size.value) - 1]

    lib.clear_optimisation_cache_dat(X, Y)
    
    return X_ret

def optimize_evs_2(X_mean, x_value, X_names, contours, coefficients, initstep=False, maximize=True, get_history=False):       
    X_mean.astype(np.float64)
    contour_arr_size = len(contours[0])
    
    X_names = np.char.encode(X_names, encoding='utf-8')
    X_names_ptr = (ctypes.c_char_p * len(X_names))()
    X_names_ptr[:] = X_names
    
    DOUBLE = ctypes.c_double
    DOUBLE_P = ctypes.POINTER(DOUBLE)
    DOUBLE_PP = ctypes.POINTER(DOUBLE_P)
    
    INT = ctypes.c_int
    INT_P = ctypes.POINTER(INT)

    X = DOUBLE_P()
    Y = DOUBLE_P()
    size = INT()
    
    lib.parameters_optimisation_2.argtypes = [np.ctypeslib.ndpointer(dtype=DOUBLE, ndim=X_mean.ndim, shape=X_mean.shape),
                                              DOUBLE,
                                              ctypes.POINTER(ctypes.c_char_p),
                                              np.ctypeslib.ndpointer(dtype=DOUBLE, ndim=contours.ndim, shape=contours.shape), 
                                              np.ctypeslib.ndpointer(dtype=DOUBLE, ndim=coefficients.ndim, shape=coefficients.shape), 
                                              INT, 
                                              ctypes.c_bool, 
                                              DOUBLE_PP, 
                                              DOUBLE_PP, 
                                              INT_P]
        
    lib.parameters_optimisation_2(X_mean,
                                  DOUBLE(x_value), 
                                  X_names_ptr, 
                                  contours, 
                                  coefficients, 
                                  INT(contour_arr_size), 
                                  ctypes.c_bool(maximize), 
                                  ctypes.byref(X), 
                                  ctypes.byref(Y), 
                                  ctypes.byref(size))
    
    print("optimization done !")
    
    if get_history is True:
        X_ret = np.ndarray([2, int(size.value)], dtype=np.float64)
        
        for i in range(int(size.value)):
            X_ret[0][i] = X[i]
            X_ret[1][i] = Y[i]
    else:
        X_ret = np.ndarray(2, dtype=np.float64)
        
        X_ret[0] = X[int(size.value) - 1]
        X_ret[1] = Y[int(size.value) - 1]
        
    lib.clear_optimisation_cache_dat(X, Y)
    
    return X_ret

def split_contours_along_y(contours):
    contour_arr_size = len(contours[0])
    
    INT = ctypes.c_int
    INT_P = ctypes.POINTER(INT)
    INT_PP = ctypes.POINTER(INT_P)
    
    contour_inf = INT_P()
    contour_sup = INT_P()
    contour_inf_size = INT()
    contour_sup_size = INT()

    lib.splitContoursAlongY.argtypes = [np.ctypeslib.ndpointer(dtype=ctypes.c_double, ndim=contours.ndim, shape=contours.shape),
                                        ctypes.c_int, INT_PP, INT_PP, INT_P, INT_P]
        
    lib.splitContoursAlongY(contours, 
                            ctypes.c_int(contour_arr_size), 
                            ctypes.byref(contour_inf), 
                            ctypes.byref(contour_sup), 
                            ctypes.byref(contour_inf_size), 
                            ctypes.byref(contour_sup_size))
    
    contour_inf_ret = np.ndarray([int(contour_inf_size.value)], dtype=np.int32)
    contour_sup_ret = np.ndarray([int(contour_sup_size.value)], dtype=np.int32)
        
    for i in range(int(contour_inf_size.value)):
        contour_inf_ret[i] = contour_inf[i]    
        
    for i in range(int(contour_sup_size.value)):
        contour_sup_ret[i] = contour_sup[i]

    lib.clear_contour_cache_dat(contour_inf, contour_sup)
    
    return contour_inf_ret, contour_sup_ret