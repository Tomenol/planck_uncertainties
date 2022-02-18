# -*- coding: utf-8 -*-
"""
Created on Tue May 25 11:02:56 2021

@author: Thomas Maynadié
"""

import numpy as np
from numpy.ctypeslib import ndpointer 
import ctypes 
import os

evs_lib_path = "EVS\\helpers\\libevs.so"

class EVSWrapper(object):
    def __init__(self):
        self.lib = ctypes.cdll.LoadLibrary(evs_lib_path)
        
        #init c++ types
        self.double_p = ctypes.POINTER(ctypes.c_double)
        self.double_pp = ndpointer(dtype=np.uintp, ndim=1, flags='C')
    
        #def wrapper functions
        self.lib.setCosmology.argtypes = [
            ctypes.c_double, 
            ctypes.c_double, 
            ctypes.c_double, 
            ctypes.c_double, 
            ctypes.c_double, 
            ctypes.c_double, 
            ctypes.c_double]
        
        self.lib.setObsParams.argtypes = [
            ctypes.c_double, 
            ctypes.c_double, 
            ctypes.c_double, 
            ctypes.c_double, 
            ctypes.c_double, 
            ctypes.c_double, 
            ctypes.c_char_p, 
            ctypes.c_char_p]
        
        self.lib.setup_evs_computation.argtypes = []
        
    def setup_evs_computation(self):
        self.lib.setup_evs_computation()
        
    def integration_over_zbins(self, IZ_BIN, z_bin_ref, halo_count_ref, LNM_SUP, LNM_INF, MASS_SUP, MASS_INF, NM_INT, DELTA_MINT):
        self.lib.integration_over_zbins(
            ctypes.c_int(IZ_BIN),
            z_bin_ref,
            halo_count_ref,
            ctypes.c_double(LNM_SUP), 
            ctypes.c_double(LNM_INF), 
            ctypes.c_double(MASS_SUP), 
            ctypes.c_double(MASS_INF), 
            ctypes.c_int(NM_INT), 
            ctypes.c_double(DELTA_MINT))
        
    def integration_over_Mmax(self, PDFHALOS_1, CUMULATIVEHALOS_1, N_CL, Z_BIN, MASS_INF, MASS_SUP, NM_INT, LNM_INF, LNM_SUP, DELTA_MINT):
        self.lib.integration_over_Mmax(
            PDFHALOS_1, 
            CUMULATIVEHALOS_1, 
            ctypes.c_double(N_CL), 
            ctypes.c_double(Z_BIN), 
            ctypes.c_double(MASS_INF), 
            ctypes.c_double(MASS_SUP), 
            ctypes.c_int(NM_INT), 
            ctypes.c_double(LNM_INF), 
            ctypes.c_double(LNM_SUP), 
            ctypes.c_double(DELTA_MINT))
        
    def find_mass_peak(self, NM_MAX, PDF, M, Mmax, pdfMax):
        DoubleArray_p = NM_MAX * ctypes.c_double
        
        PDF_p = DoubleArray_p(*(PDF.astype(np.float64)))
        M_p = DoubleArray_p(*(M.astype(np.float64)))
        
        self.lib.findMassPeak(
            ctypes.c_int(NM_MAX),
            PDF_p,
            M_p,
            Mmax,
            pdfMax)
        
    def find_confidence_interval(self, NM_MAX, cumulDistrib, M, s1_p, s1_m, s2_p, s2_m):
        DoubleArray_p = NM_MAX * ctypes.c_double
        
        cumulDistrib_p = DoubleArray_p(*(cumulDistrib.astype(np.float64)))
        M_p = DoubleArray_p(*(M.astype(np.float64)))
        
        self.lib.findMassConfIntervals(
            ctypes.c_int(NM_MAX),
            cumulDistrib_p,
            M_p,
            s1_p,
            s1_m,
            s2_p,
            s2_m)
        
    def setEVSCosmology(self, Om, Obh2, h, ns, sigma8, w_0, w_a):
        self.lib.setCosmology(
            ctypes.c_double(Om), 
            ctypes.c_double(Obh2), 
            ctypes.c_double(h), 
            ctypes.c_double(ns), 
            ctypes.c_double(sigma8), 
            ctypes.c_double(w_0), 
            ctypes.c_double(w_a))
    
    def setEVSObsParams(self, zmin, zmax, dz, delta_c, fsky, Mlim, select_mf_option, select_delta_units):
        self.lib.setObsParams(
            ctypes.c_double(zmin), 
            ctypes.c_double(zmax), 
            ctypes.c_double(dz), 
            ctypes.c_double(delta_c), 
            ctypes.c_double(fsky), 
            ctypes.c_double(Mlim), 
            select_mf_option.encode('utf-8'), 
            select_delta_units.encode('utf-8'))
