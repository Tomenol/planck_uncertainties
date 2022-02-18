# -*- coding: utf-8 -*-
"""
Created on Wed Jun 16 08:02:16 2021

@author: Thomas Maynadié
"""

import matplotlib.pyplot as plt

from optimizer import optimizerv2
from get_borders import get_borders_core
from linear_interpolation import lin_interp

# EVS core module
from EVS import evs_core

import getdist
import numpy as np

# general purpose array transpose function
def transpose(arr):
    n = np.size(arr, axis=0)
    m = np.size(arr, axis=1)
    
    arr_tmp = np.ndarray([m, n], dtype=arr.dtype)
    
    for i in range(m):
        for j in range(n):
            arr_tmp[i][j] = arr[j][i]
            
    return arr_tmp

# compute the partial derivatives of the EVS structure formation model at a given point
def compute_simple_evs_derivatives(param_values, param_indices):
    evs_ctx = evs_core.EVS()
    
    derivatives = np.ndarray([len(param_values)], dtype=np.float64)
    
    params = evs_ctx.cosmo.getParameters()
    
    getdist_to_evs = {"omegam":"Omega_m", "w":"w_0", "omegabh2":"Omega_b_h2", "sigma8":"sigma_8", "ns":"n_s", "wa":"w_a", "H0":"h"}
    
    NM_Max = 1000
    alpha = 0.1
    
    for i, param in enumerate(param_indices.keys()):
        Mmax = []
        
        print("Computing dF/d{0} : ".format(param))
        
        evs_ctx.cosmo.loadCosmology()
        old_params = evs_ctx.cosmo.getParameters()
        
        for epsilon in range(-1, 2, 2):
            new_params = old_params.copy()

            for j, parname in enumerate(param_indices.keys()):
                if parname == "H0": new_params[getdist_to_evs[parname]] = param_values[j] / 100.0
                else: new_params[getdist_to_evs[parname]] = param_values[j]

            new_params[getdist_to_evs[param]] = old_params[getdist_to_evs[param]] * (1 + alpha * epsilon)
                    
            evs_ctx.cosmo.setCosmology(**new_params)
            del new_params
                
            result = evs_ctx.evs_calculation(NM_MAX=NM_Max)

            Mmax.append(result.Mmax)
            
        Mmax = np.array(Mmax)
        
        derivatives_raw = (Mmax[1] - Mmax[0]) / (2 * old_params[getdist_to_evs[param]] * alpha)
        derivatives_raw = np.ma.masked_equal(derivatives_raw, 0)
        
        print(derivatives_raw)
        
        derivatives[i] = derivatives_raw.mean()
        
    return derivatives

# compute the partial derivatives of the EVS structure formation model
def compute_evs_derivatives(min_values, max_values, param_indices):    
    if not False in (min_values == max_values):
        print("Initial calculations. param values min=max")
        
        derivative = compute_simple_evs_derivatives(min_values, param_indices)
        derivatives = np.array([derivative, derivative], dtype=np.float64)
        
    else:
        print("Computing derivatives for min values.")
        min_derivative = compute_simple_evs_derivatives(min_values, param_indices)
        
        print("Computing derivatives for max values.")
        max_derivative = compute_simple_evs_derivatives(max_values, param_indices)
        
        derivatives = np.array([min_derivative, max_derivative], dtype=np.float64)
        
    return derivatives

def get_min_max_paramvalues(param_list, param_name_latex, samples, previous_min_value, previous_max_value, coefficients, initstep=False, conf_level=0.68):
    min_values = np.ndarray(len(param_list), dtype=np.float64)
    max_values = np.ndarray(len(param_list), dtype=np.float64)
    
    for i in range(0, len(param_list)-1):
        extracted_param_list = (param_list[i], param_list[i+1])
        
        if previous_min_value[i]    == "H0": previous_min_value[i]      = previous_min_value[i] / 100.0
        if previous_min_value[i+1]  == "H0": previous_min_value[i+1]    = previous_min_value[i+1] / 100.0
        if previous_max_value[i]    == "H0": previous_max_value[i]      = previous_max_value[i] / 100.0
        if previous_max_value[i+1]  == "H0": previous_max_value[i+1]    = previous_max_value[i+1] / 100.0
        
        if i == 0: 
            min_value, max_value = optimize1(extracted_param_list, param_name_latex, samples, previous_min_value, previous_max_value, coefficients, initstep=initstep, conf_level=conf_level)
        
        else:
            print(min_values[i])
            print(max_values[i])
            min_value, max_value = optimize2(extracted_param_list, param_name_latex, samples, min_values[i], max_values[i], coefficients, initstep=initstep, conf_level=conf_level)
        
        min_values[i]   = min_value[0]
        max_values[i]   = max_value[0]
        min_values[i+1] = min_value[1]
        max_values[i+1] = max_value[1]
        
    return min_values, max_values
        
def getPDFandContours(param_index_list, extracted_param_list, samples, conf_level=0.68):
    density = samples.get2DDensity(param_index_list[0], param_index_list[1], normalized=True)
    
    prob = np.array(transpose(density.P), dtype=np.float64)
    param1 = np.array(density.x, dtype=np.float64)
    param2 = np.array(density.y, dtype=np.float64)
    
    if extracted_param_list[0] == "H0": param1 = param1 / 100.0
    if extracted_param_list[1] == "H0": param2 = param2 / 100.0
            
    level = density.getContourLevels(contours=[conf_level]) # 1 sigma = 0.68 / 2 sigma = 0.95
    
    del density

    grid = np.array([param1, param2], dtype=np.float64)
    
    new_size = 1000
    
    new_prob = np.ndarray([new_size, new_size])
    
    new_x = np.linspace(min(param1), max(param1), new_size)
    new_y = np.linspace(min(param2), max(param2), new_size)
    new_grid = np.array([new_x, new_y])
    
    del param1
    del param2        
    
    lin_interp.interpArr2D(new_prob, new_grid, prob, grid)
    
    del grid
    del prob
    
    contour_x, contour_y, contour_xi, contour_yi = get_borders_core.get_contours(new_prob, new_grid, level)
    
    del contour_xi
    del contour_yi
    
    return new_prob, new_grid, contour_x, contour_y

def optimize1(extracted_param_list, param_name_latex, samples, previous_min_value, previous_max_value, coefficients, initstep=False, conf_level=0.68):
    extracted_param_min = (previous_min_value[0], previous_min_value[1])
    extracted_param_max = (previous_max_value[0], previous_max_value[1])
    
    print(extracted_param_list)
    
    param_index_list = list(samples.index[parname] for parname in extracted_param_list)
    
    new_prob, new_grid, contour_x, contour_y = getPDFandContours(param_index_list, extracted_param_list, samples, conf_level)

    min_value = optimizerv2.optimize_evs(extracted_param_min, extracted_param_list, new_prob, new_grid, np.array([contour_x, contour_y], dtype=np.float64), coefficients[0], initstep=initstep, maximize=False, get_history=False)
    max_value = optimizerv2.optimize_evs(extracted_param_max, extracted_param_list, new_prob, new_grid, np.array([contour_x, contour_y], dtype=np.float64), coefficients[1], initstep=initstep, maximize=True, get_history=False)
    
    print("new values param {0} = {1} (min) / {2} (max)".format(extracted_param_list[0], min_value[0], max_value[0]))
    print("new values param {0} = {1} (min) / {2} (max)".format(extracted_param_list[1], min_value[1], max_value[1]))
    
    return min_value, max_value

def optimize2(extracted_param_list, param_name_latex, samples, x_min, x_max, coefficients, initstep=False, conf_level=0.68):
    print(extracted_param_list)
    
    param_index_list = list(samples.index[parname] for parname in extracted_param_list)
    
    mean = samples.mean(param_index_list)
    
    if extracted_param_list[0] == "H0": mean[0] = mean[0] / 100.0
    if extracted_param_list[1] == "H0": mean[1] = mean[1] / 100.0
    
    new_prob, new_grid, contour_x, contour_y = getPDFandContours(param_index_list, extracted_param_list, samples, conf_level)                               
    
    del new_prob 
    del new_grid
    
    print("x values : {0} / {1}".format(x_min, x_max))
    
    print("Running optimisation proedure... (min)")
    min_value = optimizerv2.optimize_evs_2(mean.copy(), x_min, extracted_param_list, np.array([contour_x, contour_y], dtype=np.float64), coefficients[0], initstep=initstep, maximize=False, get_history=True)
 
    print("Running optimisation proedure... (max)")
    max_value = optimizerv2.optimize_evs_2(mean.copy(), x_max, extracted_param_list, np.array([contour_x, contour_y], dtype=np.float64), coefficients[1], initstep=initstep, maximize=True, get_history=True)
        
    size = len(min_value[0])
    
    print("new values param {0} = {1} (min) / {2} (max)".format(extracted_param_list[0], min_value[0][size-1], max_value[0][size-1]))
    print("new values param {0} = {1} (min) / {2} (max)".format(extracted_param_list[1], min_value[1][size-1], max_value[1][size-1]))
            
    fig = plt.figure(figsize = (1, 1), dpi=400)
    ax = fig.add_axes(rect = [0, 0, 2, 2])
    
    ax.set_xlabel(r'$' + param_name_latex[extracted_param_list[0]] + '$')
    ax.set_ylabel(r'$' + param_name_latex[extracted_param_list[1]] + '$')
            
    ax.scatter(contour_x, contour_y, s = 0.1, c="red")
    ax.scatter(min_value[0], min_value[1], s = 10, c="blue")
    ax.scatter([x_min, x_max], [mean[1], mean[1]], s = 10, c="purple")
    ax.scatter(max_value[0], max_value[1], s = 10, c="green")
    
    min_value = np.array([min_value[0][size-1], min_value[1][size-1]])
    max_value = np.array([max_value[0][size-1], max_value[1][size-1]])
    
    ax.scatter(min_value[0], min_value[1], s = 10, c="green")
    ax.scatter(max_value[0], max_value[1], s = 10, c="blue")

    ax.annotate("min", (min_value[0], min_value[1]))
    ax.annotate("max", (max_value[0], max_value[1]))
    
    return min_value, max_value

def compute_uncertainties(param_list, param_name_latex, samples, conf_level=0.68):
    param_indices = dict((parname, samples.index[parname]) for parname in param_list)
    
    min_values = samples.mean(list(param_indices.values()))
    max_values = min_values
    
    initstep = True
    
    for i in range(1):
        #coefficients = compute_evs_derivatives(min_values.copy(), max_values.copy(), param_indices)
        
        # average derivatives from previous EVS computations
        coefficients = np.array([[1678907591254780.5, 
                                 1306984500819551.2, 
                                 -637900441563705.1, 
                                 1558341652281234.2, 
                                 -183215429581535.0, 
                                 -215733125839615.06,
                                 1558341652281234.2], 
                                 [1678907591254780.5, 
                                 1306984500819551.2, 
                                 -637900441563705.1, 
                                 1558341652281234.2, 
                                 -183215429581535.0, 
                                 -215733125839615.06,
                                 1558341652281234.2]]) 
        
        for i in range(len(param_list)):
            print("dF/d{0} = {1} (min) / {2} (max)".format(param_list[i], coefficients[0][i], coefficients[1][i]))
        
        new_min_values, new_max_values = get_min_max_paramvalues(param_list, param_name_latex, samples, min_values.copy(), max_values.copy(), coefficients, initstep=initstep, conf_level=conf_level)
        
        for i in range(len(param_list)):
            print("new values param {0} = {1} (min) / {2} (max)".format(param_list[i], new_min_values[i], new_max_values[i]))
        
        if initstep is True: initstep = False
        
        if (not False in (min_values == new_min_values) and not False in (max_values == new_max_values)):
            print("No changes made during last iteration, exiting...")
            break
        else:
            min_values = new_min_values.copy()
            max_values = new_max_values.copy()
            
            print(min_values)
            print(max_values)
    
    # use coefficients in C code to set relative importance between parameters with regard to uncertainties
    # change the interface between EVS code and getDist code -> cleaner one
    # clean project roots
    # make self actualizing extremum finding algorithm
    # test algorithm and find correct ND distribution for cosmological parameters

def main():
    # script variables 
    param_list = ("omegam", "omegabh2", "H0", "sigma8", "ns", "w", "wa")
    
    param_name_latex = {"omegam":"\Omega_m", "w":"w_0", "omegabh2":"\Omega_b h^2", "sigma8":"\sigma_8", "ns":"n_s", "wa":"w_a", "H0":"H_0"}

    DATA_PATH = "./data/COM_CosmoParams_fullGrid_R3.01/base_w_wa/plikHM_TTTEEE_lowl_lowE_BAO/"
    
    samples = getdist.mcsamples.loadMCSamples(DATA_PATH + '/base_w_wa_plikHM_TTTEEE_lowl_lowE_BAO')

    min_max = compute_uncertainties(param_list, param_name_latex, samples, conf_level=0.953)
    
if __name__ == '__main__':
    main()