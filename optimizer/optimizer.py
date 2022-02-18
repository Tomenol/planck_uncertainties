# -*- coding: utf-8 -*-
"""
Created on Fri Jun 11 10:41:52 2021

@author: Thomas Maynadié
"""

import ctypes
import numpy.ctypeslib as npct
import numpy as np

import linear_interpolation.c_wrapper_helpers as c_wrapper_helpers

lib = ctypes.cdll.LoadLibrary("./optimizer/src/optimizer.so")

def optimize_evs(X_0, X_names, pdf, grid, level, maximize=True, get_history=False):       
    pdf.astype(np.float64)
    grid.astype(np.float64)
    
    n = np.size(pdf, axis=0)
    m = np.size(pdf, axis=1)
        
    lib.init_computation_cache_dat()
            
    lib.optimize_evs_gradient_descent.argtypes = [np.ctypeslib.ndpointer(dtype=ctypes.c_double, ndim=X_0.ndim, shape=X_0.shape),
                                                  ctypes.c_char_p,
                                                  ctypes.c_char_p,
                                                  np.ctypeslib.ndpointer(dtype=ctypes.c_double, ndim=pdf.ndim, shape=pdf.shape), 
                                                  np.ctypeslib.ndpointer(dtype=ctypes.c_double, ndim=grid.ndim, shape=grid.shape), 
                                                  ctypes.c_double,
                                                  ctypes.c_int, 
                                                  ctypes.c_int, 
                                                  ctypes.c_bool]
    
    print("level : " + str(level))
    
    lib.optimize_evs_gradient_descent(X_0, 
                    X_names[0].encode('utf-8'), 
                    X_names[1].encode('utf-8'), 
                    pdf, 
                    grid, 
                    ctypes.c_double(level), 
                    ctypes.c_int(n), 
                    ctypes.c_int(m), 
                    ctypes.c_bool(maximize))
    
    lib.get_convergence_results_size.restype = ctypes.c_int
    
    if get_history is True:
        conv_history_size = int(lib.get_convergence_results_size())
        
        print("array size : " + str(conv_history_size))
        
        X = np.ndarray([2, conv_history_size], dtype=np.float64)
    else:
        X = np.ndarray(2, dtype=np.float64)
    
    lib.get_convergence_results.argtypes = [np.ctypeslib.ndpointer(dtype=ctypes.c_double, ndim=X.ndim, shape=X.shape),
                                            ctypes.c_bool]
    
    lib.get_convergence_results(X, ctypes.c_bool(get_history))
    
    lib.clear_computation_cache_dat()
    
    return X

def optimize_evs2(X_0, X_names, pdf, grid, contours, maximize=True, get_history=False):       
    pdf.astype(np.float64)
    grid.astype(np.float64)
    
    n = np.size(pdf, axis=0)
    m = np.size(pdf, axis=1)
    contour_arr_size = len(contours[0])
            
    lib.init_computation_cache_dat()
            
    lib.optimize_evs_gradient_descent2.argtypes = [np.ctypeslib.ndpointer(dtype=ctypes.c_double, ndim=X_0.ndim, shape=X_0.shape),
                                                  ctypes.c_char_p,
                                                  ctypes.c_char_p,
                                                  np.ctypeslib.ndpointer(dtype=ctypes.c_double, ndim=pdf.ndim, shape=pdf.shape), 
                                                  np.ctypeslib.ndpointer(dtype=ctypes.c_double, ndim=grid.ndim, shape=grid.shape), 
                                                  np.ctypeslib.ndpointer(dtype=ctypes.c_double, ndim=contours.ndim, shape=contours.shape), 
                                                  ctypes.c_int, 
                                                  ctypes.c_int, 
                                                  ctypes.c_int, 
                                                  ctypes.c_bool]
        
    lib.optimize_evs_gradient_descent2(X_0, 
                    X_names[0].encode('utf-8'), 
                    X_names[1].encode('utf-8'), 
                    pdf,
                    grid, 
                    contours, 
                    ctypes.c_int(n), 
                    ctypes.c_int(m), 
                    ctypes.c_int(contour_arr_size), 
                    ctypes.c_bool(maximize))
    
    lib.get_convergence_results_size.restype = ctypes.c_int
    
    if get_history is True:
        conv_history_size = int(lib.get_convergence_results_size())
        
        print("array size : " + str(conv_history_size))
        
        X = np.ndarray([2, conv_history_size], dtype=np.float64)
    else:
        X = np.ndarray(2, dtype=np.float64)
    
    lib.get_convergence_results.argtypes = [np.ctypeslib.ndpointer(dtype=ctypes.c_double, ndim=X.ndim, shape=X.shape),
                                            ctypes.c_bool]
    
    lib.get_convergence_results(X, ctypes.c_bool(get_history))
    
    lib.clear_computation_cache_dat()
    
    return X

def optimize_evs3(previous_values, X_names, pdf, grid, contours, coefficients, initstep=False, maximize=True, get_history=False):       
    pdf.astype(np.float64)
    grid.astype(np.float64)
    
    n = np.size(pdf, axis=0)
    m = np.size(pdf, axis=1)
    contour_arr_size = len(contours[0])
            
    lib.init_computation_cache_dat()
    
    X = np.copy(previous_values)
            
    lib.optimize_evs_gradient_descent3.argtypes = [np.ctypeslib.ndpointer(dtype=ctypes.c_double, ndim=X.ndim, shape=X.shape),
                                                  ctypes.c_char_p,
                                                  ctypes.c_char_p,
                                                  np.ctypeslib.ndpointer(dtype=ctypes.c_double, ndim=pdf.ndim, shape=pdf.shape), 
                                                  np.ctypeslib.ndpointer(dtype=ctypes.c_double, ndim=grid.ndim, shape=grid.shape), 
                                                  np.ctypeslib.ndpointer(dtype=ctypes.c_double, ndim=contours.ndim, shape=contours.shape), 
                                                  np.ctypeslib.ndpointer(dtype=ctypes.c_double, ndim=coefficients.ndim, shape=coefficients.shape), 
                                                  ctypes.c_int, 
                                                  ctypes.c_int, 
                                                  ctypes.c_int, 
                                                  ctypes.c_bool, 
                                                  ctypes.c_bool]
        
    lib.optimize_evs_gradient_descent3(X, 
                    X_names[0].encode('utf-8'), 
                    X_names[1].encode('utf-8'), 
                    pdf,
                    grid, 
                    contours, 
                    coefficients, 
                    ctypes.c_int(n), 
                    ctypes.c_int(m), 
                    ctypes.c_int(contour_arr_size), 
                    ctypes.c_bool(maximize), 
                    ctypes.c_bool(initstep))
    
    lib.get_convergence_results_size.restype = ctypes.c_int
    
    if get_history is True:
        conv_history_size = int(lib.get_convergence_results_size())
        
        print("array size : " + str(conv_history_size))
        
        X = np.ndarray([2, conv_history_size], dtype=np.float64)
    else:
        X = np.ndarray(2, dtype=np.float64)
    
    lib.get_convergence_results.argtypes = [np.ctypeslib.ndpointer(dtype=ctypes.c_double, ndim=X.ndim, shape=X.shape),
                                            ctypes.c_bool]
    
    lib.get_convergence_results(X, ctypes.c_bool(get_history))
    
    lib.clear_computation_cache_dat()
        
    return X