#include "2D_interp.h"

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
		return 0;
	}
	else if(x == bounds[1]) // check if x closest to bounds (max)
	{
		return stop_index;
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

double interp2D(double x, double y, double* values, double* grid, double **bounds, int n, int m)
{
	int i;

	return __2D_interpolation(x, y, values, grid, bounds, n, m);
}