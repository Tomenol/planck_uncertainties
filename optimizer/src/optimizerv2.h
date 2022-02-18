#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>

#include "2D_interp.h"

#define NDIM								2

#define OMEGAM_INDEX 						0
#define OMEGABH2_INDEX 						1
#define H0_INDEX 							2
#define SIGMA8_INDEX 						3
#define NS_INDEX 							4
#define W_INDEX 							5
#define WA_INDEX 							6

double f_minim(double x);
double evs_cost_function(double *X, double *X_mean, char** param_names, double *coeffs, bool maximize);
double sgn(double x);

void evs_cost_function_grad(double *X, char** param_names, double *dX, double *coeffs, bool maximize);
void parameters_optimisation(double *X, char **X_names, double *f, double *grid, double* contours, double* coefficients, int n, int m, int contour_arr_size, bool maximize, bool firststep,double **x_ret, double **y_ret, int *ret_size);
void parameters_optimisation_2(double *X, double x_value, char **X_names, double* contours, double* coefficients, int contour_arr_size, bool maximize, double **x_ret, double **y_ret, int *ret_size);
int findClosestPointAlongDirection(int index, double *contour, int contour_arr_size, int x_sign, int y_sign);
void splitContoursAlongY(double *contour, int contour_arr_size, int **_contour_inf, int **_contour_sup, int *_size_inf, int *_size_sup);

void add_to_convergence_history(double *X, double **x_ch, double **y_ch, int *ch_size);
void add_value_to_array(int **array, int value, int *new_size);
void clear_contour_cache_dat(int *contour_sup, int *contour_inf);
void clear_optimisation_cache_dat(double *X, double *Y);
void compute_initial_step(double *X, double *gradC, char **param_names, double* contours, int contour_arr_size, bool maximize);
