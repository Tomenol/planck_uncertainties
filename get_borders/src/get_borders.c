#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>

double *x, *y;
int *x_i, *y_i;

int contour_size;

void get_point_on_square(int i, int j, double *pdf, double *grid, double level, int m);
double get_norm_factor(double *data, int n, int m);
void get_contour(double *pdf, double *grid, double level, int n, int m);
void get_point_on_vertex(int i1, int j1, int i2, int j2, double *pdf, double *grid, double level, int* flag, double *x_tmp, double *y_tmp, int m);
int get_contour_size();
void get_contour_values(double *_x_contours, double *_y_contours);
void get_contour_indices(int *_xi_contours, int *_yi_contours);

double get_norm_factor(double *data, int n, int m)
{
	int i, j;
	double max;
	
	max = 0;
	
	for(i = 0; i < n; i++)
		for(j = 0; j < m; j++)
			if (max < data[i * m + j]) max = data[i * m + j];
	
	return 1 / max;
}

/**
*	Function used to compute contours of the pdf associated with a given level based on linear interpolation
*		Returns the (x, y) points associated with the contour
**/
void get_contour(double *pdf, double *grid, double level, int n, int m)
{
	int i, j;
	
	contour_size = 0;
			
	for(i = 0; i < n-1; i++)
		for(j = 0; j < m-1; j++)
			get_point_on_square(i, j, pdf, grid, level, m);

	printf("done\n");
}

void get_point_on_square(int i, int j, double *pdf, double *grid, double level, int m)
{
	int flag;
	
	double *x_tmp, *y_tmp;
	
	x_tmp = (double*)malloc(2 * sizeof(double));
	y_tmp = (double*)malloc(2 * sizeof(double));
	
	flag = 0;
	
	get_point_on_vertex(i, j, i + 1, j, pdf, grid, level, &flag, x_tmp, y_tmp, m);
	get_point_on_vertex(i + 1, j, i + 1, j + 1, pdf, grid, level, &flag, x_tmp, y_tmp, m);
	get_point_on_vertex(i + 1, j + 1, i, j + 1, pdf, grid, level, &flag, x_tmp, y_tmp, m);
	get_point_on_vertex(i, j + 1, i, j, pdf, grid, level, &flag, x_tmp, y_tmp, m);
	
	if (flag == 2)
	{		
		contour_size++;
		
		x = realloc(x, contour_size * sizeof(double));
		y = realloc(y, contour_size * sizeof(double));
		
		x_i = realloc(x_i, contour_size * sizeof(int));
		y_i = realloc(y_i, contour_size * sizeof(int));
		
		//printf("reallocation succesful : size %d\n", contour_size);

		x[contour_size - 1] = 0.5 * (x_tmp[0] + x_tmp[1]);
		y[contour_size - 1] = 0.5 * (y_tmp[0] + y_tmp[1]);
		
		x_i[contour_size - 1] = i;
		y_i[contour_size - 1] = j;

		printf("x=%f, y=%f added (%d, %d) -> value in [%f, %f]\n", x[contour_size - 1], y[contour_size - 1], x_i[contour_size - 1], y_i[contour_size - 1], pdf[i * m + j], pdf[(i + 1) * m + j]);
	}
	
	free(x_tmp);
	free(y_tmp);
}

void get_point_on_vertex(int i1, int j1, int i2, int j2, double *pdf, double *grid, double level, int* flag, double *x_tmp, double *y_tmp, int m)
{
	double a;
	
	a = (level - pdf[i1 * m + j1]) / (pdf[i2 * m + j2] - pdf[i1 * m + j1]);
	
	if (a >= 0 && a < 1 && (*flag) < 2)
	{
		//printf("found border between (%d, %d) and (%d, %d) \n", i1, j1, i2, j2);
		
		x_tmp[(*flag)] = grid[i1];
		y_tmp[(*flag)] = grid[m + j1];
		
		//printf("x=%f, y=%f found -> value in [%f, %f]\n", x_tmp[(*flag)], y_tmp[(*flag)], pdf[i1 * m + j1], pdf[i2 * m + j2]);

		if (a != 0)
		{						
			if (i1 == i2) // vertex along x axis
			{
				y_tmp[(*flag)] = grid[m + j1] + a * (grid[m + j2] - grid[m + j1]);
			}
			else if (j1 == j2) // vertex along y axis
			{
				x_tmp[(*flag)] = grid[i1] + a * (grid[i2] - grid[i1]);
			}
			else 
			{
				printf("ERROR ! \n");
			}
		}
		
		//printf("x=%f, y=%f found -> value in [%f, %f]\n", x_tmp[(*flag)], y_tmp[(*flag)], pdf[i1 * m + j1], pdf[i2 * m + j2]);
		
		(*flag) += 1;
	}
}

int get_contour_size()
{
	if (x == NULL || y == NULL) return -1;
	else return contour_size;
}

void get_contour_values(double *_x_contours, double *_y_contours)
{
	int i;
	
	for(i = 0; i < contour_size; i++)
	{
		_x_contours[i] = x[i];
		_y_contours[i] = y[i];
		
		//printf("%d : %f %f\n", i, _x_contours[i], _y_contours[i]);
	}
}

void get_contour_indices(int *_xi_contours, int *_yi_contours)
{
	int i;
	
	for(i = 0; i < contour_size; i++)
	{
		_xi_contours[i] = x_i[i];
		_yi_contours[i] = y_i[i];
		//printf("%d : %d %d\n", i, _xi_contours[i], _yi_contours[i]);
	}
}
/**
*		This function generates contours which organization is optimized for uncertainty calculations
**/
/*void generate_organized_contours(double* _x, double grid, int n)
{
	int i, j;
	double x_min, x_max;
	
	x_min = grid[0][n - 1];
	x_max = grid[0][0];
	
	// get min and max values of x contour values ()
	for(i = 0; i < n; i++)
	{
		if ()
	}
	
	for(i = 0; i < n; i++)
}
*/
















