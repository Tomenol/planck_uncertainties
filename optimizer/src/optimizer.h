#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>

#include "2D_interp.h"

#define NDIM								2

#define MAX_STEPS 							(int)1e6
#define MAX_EVS_STEPS 						(int)1e4

#define STEP_SIZE 							0.5e-7
#define EVS_STEP_SIZE 						1e-5
#define DERIVATIVE_STEP_SIZE				1e-4
#define BOUNDARY_OPTIMIZATION_THRESHOLD		1e-5

#define OMEGAM_INDEX 						0
#define OMEGABH2_INDEX 						1
#define H0_INDEX 							2
#define SIGMA8_INDEX 						3
#define NS_INDEX 							4
#define W_INDEX 							5
#define WA_INDEX 							6

double **X_ch;
int X_ch_size;

int *i_ch;
int i_ch_size;

double boundary_cost_function(double *X, double *f, double *grid, double level, double **bounds, int n, int m);
double evs_cost_function(double *X, char **param_names, bool maximize);

void evs_cost_function_grad(double *X, char** param_names, double *dX, bool maximize);
void cost_function_grad(double *X, double *dX, double *f, double *grid, double level, double **bounds, int n, int m);
void cost_function_grad2(double *X, double *dX, double *f, double *grid, double level, double **bounds, int n, int m);
void optimize_evs_gradient_descent(double *X0, char *x_name, char *y_name, double *f, double *grid, double level, int n, int m, bool maximize);

void add_to_convergence_history(double *X);
void add_to_convergence_history_index(int i);
void get_convergence_results(double *X, bool get_history);

void clear_computation_cache_dat();
void init_computation_cache_dat();
void release_2d_cache_addr(double **ptr_arr, int n);
void compute_initial_step(double *X, double *gradC, char **param_names, double* contours, int contour_arr_size, bool maximize);
