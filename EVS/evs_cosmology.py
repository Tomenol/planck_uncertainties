# -*- coding: utf-8 -*-
"""
Created on Tue May 25 11:09:51 2021

@author: Thomas Maynadié
"""

from EVS import evs_core

class EVSCosmology(object):
    def __init__(self, core): 
        self.loadCosmology()
        self.loadSurveyCharacteristic()
        
        self.evscore = core

    def setCosmology(self, **params):
        if "Omega_m" in params:         self.Omega_m         = params["Omega_m"]
        if "Omega_b_h2" in params:      self.Omega_b_h2      = params["Omega_b_h2"]
        if "h" in params:               self.h               = params["h"]
        if "sigma_8" in params:         self.sigma_8         = params["sigma_8"]
        if "n_s" in params:             self.n_s             = params["n_s"]
        if "w_0" in params:             self.w_0             = params["w_0"]
        if "w_a" in params:             self.w_a             = params["w_a"]
        
    def setSurveyCharacteristics(self, **params):
        if "normalization" in params:   self.normalization   = params["normalization"]
        if "DELTA_C" in params:         self.DELTA_C         = params["DELTA_C"]
        if "fsky" in params:            self.fsky            = params["fsky"]
        if "dZ" in params:              self.dZ              = params["dZ"]
        if "MF_formula" in params:      self.MFformula       = params["MFformula"]
        if "M_Limit" in params:         self.Mlim            = params["Mlim"]
        if "z_min" in params:           self.z_min           = params["z_min"]
        if "z_max" in params:           self.z_max           = params["z_max"]
    
    def getCosmology(self):
        params = {
            "Omega_m": self.Omega_m,
            "Omega_b_h2": self.Omega_b_h2,
            "h": self.h,
            "sigma_8": self.sigma_8,
            "n_s": self.n_s,
            "w_0": self.w_0,
            "w_a": self.w_a
        }
        
        return params
    
    def getSurveyCharacteristics(self):
        params = {
            "normalization": self.normalization,
            "DELTA_C": self.DELTA_C,
            "dZ": self.dZ,
            "MFformula": self.MFformula,
            "Mlim": self.Mlim,
            "z_min": self.z_min,
            "z_max": self.z_max
        }
        
        return params
         
    def getParameters(self):
        params = {
            "Omega_m": self.Omega_m,
            "Omega_b_h2": self.Omega_b_h2,
            "h": self.h,
            "sigma_8": self.sigma_8,
            "n_s": self.n_s,
            "w_0": self.w_0,
            "w_a": self.w_a,
            "normalization": self.normalization,
            "DELTA_C": self.DELTA_C,
            "dZ": self.dZ,
            "MFformula": self.MFformula,
            "Mlim": self.Mlim,
            "z_min": self.z_min,
            "z_max": self.z_max
        }
        
        return params
        
    def loadSurveyCharacteristic(self, path="EVS\\config\\Survey characteristics\\default.survey"):
        try:
            config_file = open(path, "r")
        except IOError:
            evs_core.debug("Error: File does not appear to exist.", evs_core.EVSErrCode.EVS_ERROR)
            return 0
    
        data = config_file.readlines()
        
        for line in data:
            line = line.split() 
            
            if line[0] == "DELTA_UNITS" :       self.normalization          = str(line[1])                      # Omega_m def 
            if line[0] == "DELTA" :             self.DELTA_C                = float(line[1])                    # Omega_b h2
            if line[0] == "F_SKY" :             self.fsky                   = float(line[1])                    # h
            if line[0] == "DELTA_Z" :           self.dZ                     = float(line[1])                    # sigma_8
            if line[0] == "MF_FORMULA" :        self.MFformula              = str(line[1])                      # n_s
            if line[0] == "M_LIM" :             self.Mlim                   = float(line[1])                    # w_0
            if line[0] == "Z_MIN" :             self.z_min                  = float(line[1])                    # w_a
            if line[0] == "Z_MAX" :             self.z_max                  = float(line[1])                    # w_a

    def loadCosmology(self, path="EVS\\config\\Cosmologies\\planck2018.cosmo"):
        try:
            config_file = open(path, "r")
        except IOError:
            evs_core.debug("Error: File does not appear to exist.", evs_core.EVSErrCode.EVS_ERROR)
            return 0
        
        data = config_file.readlines()
        
        for line in data:
            line = line.split() 
            
            if line[0] == "OMEGA_M" :       self.Omega_m        = float(line[1])                    # Omega_m def 
            if line[0] == "OMEGA_B_HH" :    self.Omega_b_h2     = float(line[1])                    # Omega_b h2
            if line[0] == "H" :             self.h              = float(line[1])                    # h
            if line[0] == "SIGMA_8" :       self.sigma_8        = float(line[1])                    # sigma_8
            if line[0] == "ANS" :           self.n_s            = float(line[1])                    # n_s
            if line[0] == "W_0" :           self.w_0            = float(line[1])                    # w_0
            if line[0] == "W_A" :           self.w_a            = float(line[1])                    # w_a
        
    def writeParametersToFile(self, path, createnew=True):
        try:
            config_file = open(path, "r")
        except IOError:
            evs_core.debug("Error: File does not appear to exist.", evs_core.EVSErrCode.EVS_ERROR)
            
            if (createnew == True):
                evs_core.debug("Creating new file (ow = true)", evs_core.EVSErrCode.EVS_STATUS)
                 
                config_file = open(path, "w+")
            else:
                return -1
                 
        config_file.write(str(  self.Omega_m)         + "\n")
        config_file.write(str(  self.Omega_b_h2)      + "\n")
        config_file.write(str(  self.h)               + "\n")
        config_file.write(str(  self.sigma_8)         + "\n")
        config_file.write(str(  self.n_s)             + "\n")
        config_file.write(str(  self.w_0)             + "\n")
        config_file.write(str(  self.w_a)             + "\n")
        config_file.write(      self.normalization    + "\n")
        config_file.write(str(  self.DELTA_C)         + "\n")
        config_file.write(str(  self.fsky)            + "\n")
        config_file.write(str(  self.dZ)              + "\n")
        config_file.write(      self.MFformula        + "\n")
        config_file.write(      self.Mlim             + "\n")
        config_file.write(str(  self.z_min)           + "\n")
        config_file.write(str(  self.z_max)           + "\n")
        config_file.write(      self.nz_outfile       + "\n")
        config_file.write(      self.F_outfile        + "\n")
        config_file.write(      self.PDF_outfile      + "\n")
                
    def __setEVSCosmology(self):
        self.evscore.wrapper.setEVSCosmology(
            self.Omega_m, 
            self.Omega_b_h2, 
            self.h, 
            self.n_s, 
            self.sigma_8, 
            self.w_0, 
            self.w_a)
        
    def __setEVSObsParameters(self):
        self.evscore.wrapper.setEVSObsParams(
            self.z_min, 
            self.z_max, 
            self.dZ, 
            self.DELTA_C, 
            self.fsky,
            self.Mlim,
            self.MFformula, 
            self.normalization)
        
    def updateParameters(self):
        self.__setEVSObsParameters()
        self.__setEVSCosmology()        