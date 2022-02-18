#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>

#define NDIM 	2

// Function prototypes
void convert_xy_to_ij(double x, double y, int *i, int *j, double* grid, double **bounds, int n, int m);
int closest_index(double *arr, double x, double *bounds, int start_index, int stop_index);
void convert_ij_to_xy(int i, int j, double *x, double *y, double *grid, int n, int m);
void get_bounds(double *grid, double **bounds, int n, int m);
double __2D_interpolation(double x, double y, double* values, double *grid, double **bounds, int n, int m);
double interp2D(double x, double y, double* values, double* grid, int n, int m);
void interp_array_2D(double* _new_grid_values, double* _new_grid, double* _old_grid_values, double* _old_grid, int _new_n, int _new_m, int _old_n, int _old_m);

// Function declaration

void convert_xy_to_ij(double x, double y, int *i, int *j, double* grid, double **bounds, int n, int m)
{
	(*i) = (int) closest_index(grid, x, bounds[0], 0, m - 1);
	(*j) = (int) closest_index(grid, y, bounds[1], m, 2 * m - 1) - m;
}

int closest_index(double *arr, double x, double *bounds, int start_index, int stop_index)
{
	int i_inf, i_sup, i_new;
	
	if(x == bounds[0]) // check if x closest to bounds (min)
	{
		return start_index;
	}
	else if(x == bounds[1]) // check if x closest to bounds (max)
	{
		return start_index + stop_index;
	}
	else
	{
		i_inf = start_index;
		i_sup = start_index + stop_index;
		
		goto1:
		i_new = (int)(0.5 * (double)(i_sup + i_inf));
		
		if (i_new == i_inf) return i_inf;
		
		if (arr[i_new] < x) 
		{
			i_inf = i_new;
			goto goto1;
		}
		else
		{
			i_sup = i_new;
			goto goto1;
		}
	}
}

void convert_ij_to_xy(int i, int j, double *x, double *y, double *grid, int n, int m)
{
	if (i >= n || i < 0 || j >= m || j < 0) printf("Error : index out of range.\n");
	else
	{
		(*x) = grid[i];
		(*y) = grid[m + j];
	}
}

void get_bounds(double *grid, double **bounds, int n, int m)
{
	int i;
	
	for(i = 0; i < NDIM; i++)
	{
		bounds[i][0] = grid[i * m]; 			// min 
		bounds[i][1] = grid[(i + 1) * m - 1]; 	// max 
	}
	
	printf("xmin=%f xmax=%f ymin=%f ymax=%f\n", bounds[0][0], bounds[0][1], bounds[1][0], bounds[1][1]);
}

double __2D_interpolation(double x, double y, double* values, double *grid, double **bounds, int n, int m)
{
	int i1, j1;
	double x1, x2, y1, y2;
	double interp_result;

	convert_xy_to_ij(x, y, &i1, &j1, grid, bounds, n, m);
	
	convert_ij_to_xy(i1, j1, &x1, &y1, grid, n, m);
	convert_ij_to_xy(i1 + 1, j1 + 1, &x2, &y2, grid, n, m);
	
	interp_result = (x2 - x1) * (y2 - y1);
	interp_result = (	(x2 - x) * (y2 - y) * values[i1 * m 		+ j1] 
					+ 	(x - x1) * (y2 - y) * values[(i1 + 1) * m  	+ j1] 
					+ 	(x2 - x) * (y - y1) * values[i1  * m 		+ j1 + 1] 
					+ 	(x - x1) * (y - y1) * values[(i1 + 1) * m 	+ j1 + 1]) 
					/ interp_result;
	
	return interp_result;
}

double interp2D(double x, double y, double* values, double* grid, int n, int m)
{
	int i, j;
	
	double **bounds = (double**)malloc(NDIM * sizeof(double*));
	
	for(i = 0; i < NDIM; i++)
		bounds[i] = (double*)malloc(2 * sizeof(double));
	
	get_bounds(grid, bounds, n, m);

	return __2D_interpolation(x, y, values, grid, bounds, n, m);
}

void interp_array_2D(double* _new_grid_values, double* _new_grid, double* _old_grid_values, double* _old_grid, int _new_n, int _new_m, int _old_n, int _old_m)
{	
	int i, j;	
		
	double x, y;
	
	double **bounds = (double**)malloc(NDIM * sizeof(double*));
	
	for(i = 0; i < NDIM; i++)
		bounds[i] = (double*)malloc(2 * sizeof(double));
	
	get_bounds(_new_grid, bounds, _new_n, _new_m);
	
	// interp 2D
	for(i = 0; i < _new_n - 1; i++)
	{
		for(j = 0; j <_new_m - 1; j++) 
		{
			x = _new_grid[i];
			y = _new_grid[_new_m + j];
						
			_new_grid_values[i * _new_m + j] = __2D_interpolation(x, y, _old_grid_values, _old_grid, bounds, _old_n, _old_m);
		}
	}
}










