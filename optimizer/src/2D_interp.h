#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define NDIM 2

int closest_index(double *arr, double x, double *bounds, int start_index, int stop_index);

double interp2D(double x, double y, double* values, double* grid, double **bounds, int n, int m);
double __2D_interpolation(double x, double y, double* values, double *grid, double **bounds, int n, int m);

void convert_ij_to_xy(int i, int j, double *x, double *y, double *grid, int n, int m);
void convert_xy_to_ij(double x, double y, int *i, int *j, double* grid, double **bounds, int n, int m);