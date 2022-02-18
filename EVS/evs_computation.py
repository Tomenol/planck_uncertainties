# -*- coding: utf-8 -*-
"""
Created on Wed May 26 13:09:48 2021

@author: Thomas Maynadié
"""

from EVS import evs_core
from EVS import evs_threading
from EVS.evs_debug import debug, EVSErrCode

import numpy as np
from numpy import ndarray
import ctypes
import math

import os

# placeholder for evs computation results
# todo add id or informations relative to the computation  (cosmo params, etc)
class EVSResults(object):
    def __init__(self): pass
        
    def getResults(self): return (self.z_bins, self.cluster_counts, self.pdf, self.cumul)

    def setRawResults(self, _z_bins, _cluster_counts, _pdf, _cumul, _m): 
        self.cluster_counts = self.__convert2float(_cluster_counts)
        self.z_bins = self.__convert2float(_z_bins)
        
        self.M = self.__convert2float(_m)
        
        self.pdf = self.__convert2float(_pdf)
        self.cumul = self.__convert2float(_cumul)
        
        self.z_size = len(self.z_bins)
        self.m_size = len(self.M)
                
    def setProcessedResults(self, _Mmax, _PDFmax, _s1m, _s1p, _s2m, _s2p): 
        self.Mmax =  self.__convert2float(_Mmax)
        self.PDFmax = self.__convert2float( _PDFmax)
        
        self.s1m =  self.__convert2float(_s1m)
        self.s1p =  self.__convert2float(_s1p)
        self.s2m =  self.__convert2float(_s2m)
        self.s2p =  self.__convert2float(_s2p)
        
        self.PDFnorm = self.pdf
        
        for i in range(len(self.PDFnorm)):
            for j in range(len(self.PDFnorm[i])):
                self.PDFnorm[i][j] = self.PDFnorm[i][j] / self.PDFmax[i]

    def __convert2float(self, _val):
        if not hasattr(_val, "__len__"): 
            if hasattr(_val, "value"): 
                return float(_val.value)
        else:
            for i in range(len(_val)):
                _val[i] = self.__convert2float(_val[i])
                 
        return _val
        
    def setParameters(self, params):
        self.params = params
    
    def printResultsToOutputFile(self, path, filename):
        iflag = True
    
        if not (hasattr(self, "z_bins") and hasattr(self, "cluster_counts") and hasattr(self, "Mmax") and hasattr(self, "s1m") and hasattr(self, "s1p") and hasattr(self, "s2m") and hasattr(self, "s2p") and hasattr(self, "pdf") and hasattr(self, "cluster_counts") and hasattr(self, "z_bins") and hasattr(self, "M") and hasattr(self, "cumul")): 
            debug("EVS computation results are not initialized, please use setResults(...) before trying to print them.", evs_core.EVSErrCode.EVS_ERROR)
            iflag = False
        elif not os.path.exists(path): 
            debug("output destination doesn't exist '{0}'".format(path), evs_core.EVSErrCode.EVS_ERROR)
            iflag = False
        else:
            try: output_file = open(path + "\\" + filename, "w")
            except OSError: 
                debug("could not create new result file '{0}' in '{1}'".format(filename, path), evs_core.EVSErrCode.EVS_ERROR)
                iflag = False
        
        if iflag == True:
            output_file.write("EVS Results : " + filename + "\n\n")
            output_file.write("zbin number : " + str(self.z_size) + "\n")
            output_file.write("mass number : " + str(self.m_size) + "\n\n")
            
            output_file.write("Parameters :\n")
            
            for parameter in self.params.items():
                if isinstance(parameter[1], str):
                    output_file.write(parameter[0] + " " + parameter[1] + "\n")
                else:
                    output_file.write(parameter[0] + " " + '{:<10.5f}'.format(parameter[1]) + "\n")
            
            output_file.write("\n\n\nRAW data :\n\n")
            
            output_file.write("ZBINS[] :\n")
            for i in range(self.z_size):
                output_file.write('{:<30.20e}'.format(self.z_bins[i]) + " ")
                
            output_file.write("\n\nNCL[] :\n")
            for i in range(self.z_size):
                output_file.write('{:<30.20e}'.format(self.cluster_counts[i]) + " ")
            
            output_file.write("\n\nM[] :\n")
            for i in range(self.m_size):   
                output_file.write('{:<30.20e}'.format(self.M[i]) + " ")
                                    
            output_file.write("\n\nPDF[][] :\n")
            for i in range(self.z_size):
                for j in range(self.m_size):                                          
                    output_file.write('{:<30.20e}'.format(self.pdf[i][j]) + " ")
                    
                output_file.write("\n")   
                    
            output_file.write("\n\nCumulativeDistribution[][] :\n")
            for i in range(self.z_size):
                for j in range(self.m_size):
                    output_file.write('{:<30.20e}'.format(self.cumul[i][j]) + " ")
                    
                output_file.write("\n")
                
            output_file.write("\n\n\nProcessed data :\n\n")
                            
            output_file.write("\n\nPDFMax[] :\n")
            for i in range(self.z_size):
                output_file.write('{:<30.20e}'.format(self.PDFmax[i]) + " ")
                
            output_file.write("\n\nMmax[] :\n")
            for i in range(self.z_size):
                output_file.write('{:<30.20e}'.format(self.Mmax[i]) + " ")
                                    
            output_file.write("\n\ns1m[] :\n")
            for i in range(self.z_size):
                output_file.write('{:<30.20e}'.format(self.s1m[i]) + " ")
                
            output_file.write("\n\ns1p[] :\n")
            for i in range(self.z_size):
                output_file.write('{:<30.20e}'.format(self.s1p[i]) + " ")
                
            output_file.write("\n\ns2m[] :\n")
            for i in range(self.z_size):
                output_file.write('{:<30.20e}'.format(self.s2m[i]) + " ")
                
            output_file.write("\n\ns2p[] :\n")
            for i in range(self.z_size):
                output_file.write('{:<30.20e}'.format(self.s2p[i]) + " ")
                
            output_file.write("\n\nPDF Normalized[][] :\n")
            for i in range(self.z_size):
                for j in range(self.m_size):
                    output_file.write('{:<30.20e}'.format(self.pdf[i][j]/self.PDFmax[i]) + " ")
                
                output_file.write("\n")
                
            debug("Results exported to file '{0}' at '{1}'".format(filename, path), evs_core.EVSErrCode.EVS_STATUS)
            output_file.close()
            
    def importResultsFromFile(self, filepath, filename):
        
        if not os.path.exists(filepath): 
            debug("file destination doesn't exist '{0}'".format(filepath), evs_core.EVSErrCode.EVS_ERROR)
        else:
            try: output_file = open(filepath + "\\" + filename, "r")
            except OSError: 
                debug("Results exported to file '{0}' at '{1}'".format(filename, filepath), evs_core.EVSErrCode.EVS_ERROR)
            
            file_data = output_file.readlines()
            
            output_file.close()
            
            self.params = {}
            
            self.z_bins = []
            self.cluster_counts = []
            self.M = []                            
            self.pdf = []
            self.cumul = []
            
            self.Mmax = []
            self.PDFmax = []
            
            self.s1m =  []
            self.s1p =  []
            self.s2m =  []
            self.s2p =  []
            
            self.pdfNorm = []
            
            for i in range(len(file_data)):
                line = file_data[i].split()
                
                n = len(line)
                
                if n > 0:                    
                    if line[0] == "zbin":           self.zsize = int(line[n - 1])
                    if line[0] == "mass":           self.msize = int(line[n - 1])
                    
                    if line[0] == "Omega_m":        self.params["Omega_m"] = float(line[n - 1])
                    if line[0] == "Omega_b_h2":     self.params["Omega_b_h2"] = float(line[n - 1])
                    if line[0] == "h":              self.params["h"] = float(line[n - 1])
                    if line[0] == "sigma_8":        self.params["sigma_8"] = float(line[n - 1])
                    if line[0] == "n_s":            self.params["n_s"] = float(line[n - 1])
                    if line[0] == "w_0":            self.params["w_0"] = float(line[n - 1])
                    if line[0] == "w_a":            self.params["w_a"] = float(line[n - 1])
                    if line[0] == "normalization":  self.params["normalization"] = str(line[n - 1])
                    if line[0] == "dZ":             self.params["dZ"] = float(line[n - 1])
                    if line[0] == "MFformula":      self.params["MFformula"] = str(line[n - 1])
                    if line[0] == "Mlim":           self.params["Mlim"] = float(line[n - 1])
                    if line[0] == "z_min":          self.params["z_min"] = float(line[n - 1])
                    if line[0] == "z_max":          self.params["z_max"] = float(line[n - 1])
                    
                    if line[0] == "nz_outfile":     self.params["nz_outfile"] = str(line[n - 1])
                    if line[0] == "F_outfile":      self.params["F_outfile"] = str(line[n - 1])
                    if line[0] == "PDF_outfile":    self.params["PDF_outfile"] = str(line[n - 1])
                    
                    if line[0] == "ZBINS[]": 
                        line = file_data[i + 1].split()
                        
                        for j in range(len(line)):
                            self.z_bins.append(float(line[j]))
                            
                        i = i + 1
                        
                        
                    if line[0] == "NCL[]":
                        line = file_data[i + 1].split()
                        
                        for j in range(len(line)):
                            self.cluster_counts.append(float(line[j]))
                            
                        i = i + 1
                    
                    if line[0] == "M[]":
                        line = file_data[i + 1].split()
                        
                        for j in range(len(line)):
                            self.M.append(float(line[j]))
                            
                        i = i + 1
                    
                    if line[0] == "PDF[][]":
                        for j in range(self.zsize):
                            line = file_data[i + j + 1].split()
                            
                            self.pdf.append([])
                            
                            for k in range(len(line)):
                                self.pdf[j].append(float(line[k]))
                            
                        i = i + self.zsize
                    
                    if line[0] == "CumulativeDistribution[][]": 
                        for j in range(self.zsize):
                            line = file_data[i + j + 1].split()
                            
                            self.cumul.append([])
                            
                            for k in range(len(line)):
                                self.cumul[j].append(float(line[k]))
                
                        i = i + self.zsize
 
                    if line[0] == "PDFMax[]": 
                        line = file_data[i + 1].split()
                        
                        for j in range(len(line)):
                            self.PDFmax.append(float(line[j]))
                            
                        i = i + 1
                        
                        
                    if line[0] == "Mmax[]":
                        line = file_data[i + 1].split()
                        
                        for j in range(len(line)):
                            self.Mmax.append(float(line[j]))
                            
                        i = i + 1
                    
                    if line[0] == "s1m[]":
                        line = file_data[i + 1].split()
                        
                        for j in range(len(line)):
                            self.s1m.append(float(line[j]))
                            
                        i = i + 1
                        
                    if line[0] == "s1p[]":
                        line = file_data[i + 1].split()
                        
                        for j in range(len(line)):
                            self.s1p.append(float(line[j]))
                            
                        i = i + 1
                    
                    if line[0] == "s2m[]":
                        line = file_data[i + 1].split()
                        
                        for j in range(len(line)):
                            self.s2m.append(float(line[j]))
                            
                        i = i + 1
                    
                    if line[0] == "s2p[]":
                        line = file_data[i + 1].split()
                        
                        for j in range(len(line)):
                            self.s2p.append(float(line[j]))
                            
                        i = i + 1
                    
                    if line[0] == "PDF":
                        for j in range(self.zsize):
                            line = file_data[i + j + 1].split()
                            
                            self.pdfNorm.append([])
                            
                            for k in range(len(line)):
                                self.pdfNorm[j].append(float(line[k]))
                            
                        i = i + self.zsize
                
            debug("Results imported from file '{0}' at '{1}'".format(filename, filepath), evs_core.EVSErrCode.EVS_STATUS)
            
# evs computation class
class EVSComputation(object):
    def __init__(self): pass

    def compute(self, ctx): pass
    def get_results(self): pass

    #definitions    
    
    # NM_MAX sets the number of points for F(M, z)
    def compute(self, ctx, NM_MAX=1000):
        debug("Initializing evs computation...", evs_core.EVSErrCode.EVS_STATUS)
        
        threadLimit = 300
        
        MDELTA_MIN = ctx.cosmo.Mlim
        ZMAX = ctx.cosmo.z_max
        ZMIN = ctx.cosmo.z_min
        DZ = ctx.cosmo.dZ
        
        NM_INT = 5000 # NM_MAX sets the number of computation points for the integral of N(z)
        NM_ARRAY_DIM = NM_MAX 

        MASS_SUP_MAX = 1e+16
        
        NZ_BIN = int(math.floor((ZMAX - ZMIN)/DZ)) + 1

        threadManager = evs_threading.EVSThreadManager(NZ_BIN)
            
        results = EVSResults()
        
        debug("Done.", evs_core.EVSErrCode.EVS_STATUS)
        
        # ----------------------------------------------
        #               COMPUTE N_cl(Z) 
        # ----------------------------------------------
        debug("Computing cluster number count N_cl(z)...", evs_core.EVSErrCode.EVS_STATUS)
        
        MASS_INF = MDELTA_MIN * 0.8
        MASS_SUP = MASS_SUP_MAX
        LNM_INF = np.log(MASS_INF)
        LNM_SUP = np.log(MASS_SUP)
    
        DELTA_MINT = (LNM_SUP - LNM_INF)/float(NM_INT - 1)
        
        z_bins = np.empty(NZ_BIN, dtype=object)
        cluster_counts = np.empty(NZ_BIN, dtype=object)
        M = np.empty(NM_MAX, dtype=object)
    
        for IZ_BIN in range(NZ_BIN):
            z_bin = ctypes.c_double()
            cluster_count = ctypes.c_double()
            
            z_bin_ref = ctypes.byref(z_bin)
            halo_count_ref = ctypes.byref(cluster_count)
            
            threadManager.createTask(IZ_BIN, ctx.wrapper.integration_over_zbins, (IZ_BIN, z_bin_ref, halo_count_ref, LNM_SUP, LNM_INF, MASS_SUP, MASS_INF, NM_INT, DELTA_MINT))
            
            z_bins[IZ_BIN] = z_bin
            cluster_counts[IZ_BIN] = cluster_count
            
        threadManager.startSession(True)
        debug("Done.", evs_core.EVSErrCode.EVS_STATUS)
        
        # ----------------------------------------------
        #               COMPUTE F(M, Z) 
        # ----------------------------------------------
        
        debug("Computing F(M, z)...", evs_core.EVSErrCode.EVS_STATUS)
        
        pdf = np.empty((NZ_BIN, NM_MAX), dtype=object)
        cumulative_distribution = np.empty((NZ_BIN, NM_MAX), dtype=object)
    
        MASS_SUP_MIN = MDELTA_MIN
        
        LGM_SUP_MIN = np.log10(MASS_SUP_MIN)
        LGM_SUP_MAX = np.log10(MASS_SUP_MAX)
        
        del threadManager
        
        threadManager = evs_threading.EVSThreadManager(threadLimit)
        
        for IM_MAX in range(NM_MAX):
            MASS_SUP = 10.0**(LGM_SUP_MIN + (LGM_SUP_MAX-LGM_SUP_MIN) * float(IM_MAX)/float(NM_MAX - 1))
            M[IM_MAX] = MASS_SUP
        
            LNM_SUP = np.log(MASS_SUP)
            DELTA_MINT = (LNM_SUP - LNM_INF)/float(NM_INT-1)
            
            for IZ_BIN in range(NZ_BIN):
                PDF_z = ctypes.c_double()
                cumulDistribution_z = ctypes.c_double()
                
                PDF_z_ref = ctypes.byref(PDF_z)
                cumulDistribution_z_ref = ctypes.byref(cumulDistribution_z)
                
                threadManager.createTask(IZ_BIN, ctx.wrapper.integration_over_Mmax, (PDF_z_ref, cumulDistribution_z_ref, cluster_counts[IZ_BIN].value, z_bins[IZ_BIN].value, MASS_INF, MASS_SUP, NM_INT, LNM_INF, LNM_SUP, DELTA_MINT))
                
                pdf[IZ_BIN][IM_MAX] = PDF_z
                cumulative_distribution[IZ_BIN][IM_MAX] = cumulDistribution_z
                
        threadManager.startSession(True)
        
        debug("Done.", evs_core.EVSErrCode.EVS_STATUS)
        
        results.setRawResults(z_bins, cluster_counts, pdf, cumulative_distribution, M) 

        # ----------------------------------------------
        #               PROCESSING 
        # ----------------------------------------------
        
        debug("Processing results...", evs_core.EVSErrCode.EVS_STATUS)
        
        Mmax = []
        PDFmax = []
        s1m = []
        s1p = []
        s2m = []
        s2p = []
        
        del threadManager
        
        threadManager = evs_threading.EVSThreadManager(2 * NZ_BIN)

        for i in range(NZ_BIN):
            pdfMax_z = ctypes.c_double()
            pdfMax_z_ref = ctypes.byref(pdfMax_z)
            
            Mmax_z = ctypes.c_double()
            Mmax_z_ref = ctypes.byref(Mmax_z)
            
            s1m_z = ctypes.c_double()
            s1m_z_ref = ctypes.byref(s1m_z)
            
            s1p_z = ctypes.c_double()
            s1p_z_ref = ctypes.byref(s1p_z)
            
            s2m_z = ctypes.c_double()
            s2m_z_ref = ctypes.byref(s2m_z)
            
            s2p_z = ctypes.c_double()
            s2p_z_ref = ctypes.byref(s2p_z)

            threadManager.createTask(IZ_BIN, ctx.wrapper.find_mass_peak, (NM_MAX, results.pdf[i], results.M, Mmax_z_ref, pdfMax_z_ref))
            threadManager.createTask(IZ_BIN, ctx.wrapper.find_confidence_interval, (NM_MAX, results.cumul[i], results.M, s1p_z_ref, s1m_z_ref, s2p_z_ref, s2m_z_ref))
            
            Mmax.append(Mmax_z)
            PDFmax.append(pdfMax_z)
            s1m.append(s1m_z)
            s1p.append(s1p_z)
            s2m.append(s2m_z)
            s2p.append(s2p_z)
            
        threadManager.startSession(True)
        del threadManager
        
        debug("Done.", evs_core.EVSErrCode.EVS_STATUS)
        
        results.setProcessedResults(Mmax, PDFmax, s1m, s1p, s2m, s2p) 
        results.setParameters(ctx.cosmo.getParameters())
        
        return results