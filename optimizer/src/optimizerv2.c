#include "optimizerv2.h"

/* 
*		function that has to be minimize in order to maximize a given parameter x
*/
double f_minim(double x)
{
	return -x;
}

double evs_cost_function(double *X, double *X_mean, char** param_names, double *coeffs, bool maximize)
{
	int param_index, index;
	double coeff_sum, cost_fnc, cost, coeff;
	
	coeff = -1.0;
	cost_fnc = 0.0;
	
	if(maximize == true) coeff = 1.0;
	
	coeff_sum = 0;
	for (param_index = 0; param_index < 7; param_index++)
		if(!isnan(coeffs[param_index])) coeff_sum += fabs(coeffs[param_index]);

	for (param_index = 0; param_index < NDIM; param_index++)
	{
		if (strcmp(param_names[param_index], "omegam") == 0) index = OMEGAM_INDEX;
		if (strcmp(param_names[param_index], "omegabh2") == 0) index = OMEGABH2_INDEX;
		if (strcmp(param_names[param_index], "H0") == 0) index = H0_INDEX;
		if (strcmp(param_names[param_index], "sigma8") == 0) index = SIGMA8_INDEX;
		if (strcmp(param_names[param_index], "ns") == 0) index = NS_INDEX;
		if (strcmp(param_names[param_index], "w") == 0) index = W_INDEX;
		if (strcmp(param_names[param_index], "wa") == 0) index = WA_INDEX;
		
		cost = f_minim(coeff * (X[param_index] - X_mean[param_index]) * coeffs[index] / coeff_sum);
		
		printf("Cost function : parameter %s (mean %f / val %f)-> derivative %f : cost %f\n", param_names[param_index], X_mean[param_index], X[param_index], coeffs[index], cost);
		
		cost_fnc += cost;
	}
	
	return cost_fnc;
}

void evs_cost_function_grad(double *X, char** param_names, double *dX, double *coeffs, bool maximize)
{
	int param_index, index;
	double coeff_sum, coeff;
	
	if(maximize == false) coeff = 1.0;
	else coeff = -1.0;
	
	coeff_sum = 0;
	for (param_index = 0; param_index < 7; param_index++)
		if(!isnan(coeffs[param_index])) coeff_sum += fabs(coeffs[param_index]);

	for (param_index = 0; param_index < NDIM; param_index++)
	{
		if (strcmp(param_names[param_index], "omegam") == 0) index = OMEGAM_INDEX;
		if (strcmp(param_names[param_index], "omegabh2") == 0) index = OMEGABH2_INDEX;
		if (strcmp(param_names[param_index], "H0") == 0) index = H0_INDEX;
		if (strcmp(param_names[param_index], "sigma8") == 0) index = SIGMA8_INDEX;
		if (strcmp(param_names[param_index], "ns") == 0) index = NS_INDEX;
		if (strcmp(param_names[param_index], "w") == 0) index = W_INDEX;
		if (strcmp(param_names[param_index], "wa") == 0) index = WA_INDEX;
		
		dX[param_index] = coeff * coeffs[index] / coeff_sum;
	}
}

void parameters_optimisation(double *X, char **X_names, double *f, double *grid, double* contours, double* coefficients, int n, int m, int contour_arr_size, bool maximize, bool firststep,double **x_ret, double **y_ret, int *ret_size)
{
	int i, j, k;
	int index, iflag;
	int i1, i2;
	
	double min_dist;
	double dotp1, dotp2, dist;
	double *gradC;
	
	double *x_ch, *y_ch;
	int ch_size;
	
	ch_size = 0;
	add_to_convergence_history(X, &x_ch, &y_ch, &ch_size);
	
	printf("Executing EVS optimisation script\n");
	
	gradC = (double*)malloc(NDIM * sizeof(double));
		
	evs_cost_function_grad(X, X_names, gradC, coefficients, maximize);

	if (firststep == true) 
	{
		compute_initial_step(X, gradC, X_names, contours, contour_arr_size, maximize);
		add_to_convergence_history(X, &x_ch, &y_ch, &ch_size);
	}

	CONVERGENCE_LOOP:
	printf("Convergence loop\n");
	// get 2 oposite closest points to X
	min_dist = 1e16;
	
	for(j = 0; j < contour_arr_size; j++)
	{		
		dist = 0;
	
		for(k = 0; k < NDIM; k++)
		{
			dist += pow(X[k] - contours[k * contour_arr_size + j], 2);
		}
		
		dist = sqrt(dist);
		
		if(dist < min_dist)
		{
			min_dist = dist;
			i1 = j;
		}
	}
	
	min_dist = 1e16;
	
	for(j = 0; j < contour_arr_size; j++)
	{
		if(j != i1)
		{
			dist = 0;
			dotp1 = 0;
		
			for(k = 0; k < NDIM; k++)
			{
				dist += pow(X[k] - contours[k * contour_arr_size + j], 2);
				dotp1 += (X[k] - contours[k * contour_arr_size + j]) * (X[i1] - contours[k * contour_arr_size + i1]);
			}

			dist = sqrt(dist);
			
			if(dist < min_dist && dotp1 < 0)
			{
				min_dist = dist;
				i2 = j;
				
				//printf("dot product : %f (M1 : x=%f, y=%f / M2 : x=%f, y=%f)\n", dotp1, contours[i1], contours[contour_arr_size + i1], contours[i2], contours[contour_arr_size + i2]);
			}
		}
	}
	
	dotp1 = 0;
	dotp2 = 0;
	
	for(i = 0; i < NDIM; i++)
	{
		dotp1 += (X[i] - contours[i * contour_arr_size + i1]) * gradC[i];
		dotp2 += (X[i] - contours[i * contour_arr_size + i2]) * gradC[i];
	}
	
	//printf("dot products : dotp1=%f / dotp2=%f\n", dotp1, dotp2);
	
	// test 
	if (dotp1 <= 0 && dotp2 <= 0)
	{
		//printf("dead zone, exiting...\n");
		iflag = 2;
	}
	else if(dotp1 < 0 && dotp2 > 0)
	{
		//printf("going towards M2.\n");
		
		index = i2;
		iflag = 1;
	}
	else if(dotp1 > 0 && dotp2 < 0)
	{
		//printf("going towards M1.\n");
		
		index = i1;
		iflag = 1;
	}
	
	if(iflag == 1)
	{
		for(i = 0; i < NDIM; i++)
			X[i] = contours[i * contour_arr_size + index];
		
		add_to_convergence_history(X, &x_ch, &y_ch, &ch_size);
		
		goto CONVERGENCE_LOOP;
	}
	
	//printf("interpolation result : x=%f y=%f\n", X[0], X[1]);
	
	(*x_ret) = x_ch;
	(*y_ret) = y_ch;
	(*ret_size) = ch_size;
	
	free(gradC);
}

double sgn(double x)
{
	if (x >= 0.0) return 1.0;
	if (x < 0.0) return -1.0;
}

void get_contour_intersection(double *contours, int contour_arr_size, double x_value, int *contour_index, int contour_index_size, double *x, double *y)
{
	int i, i1, i2;
	double dist;
	
	for(i = 0; i < contour_index_size; i++)
	{
		dist = x_value - contours[contour_index[i]];
		
		//printf("i=%d -> dist = %f (%f - %f)\n", i, dist, x_value, contours[contour_index[i]]);
		
		i1 = i;
		i2 = i;
		
		if (dist <= 0)
		{
			if(dist < 0 && i > 0) i1 = i - 1;
			break;
		}
	}
	
	(*x) = x_value;
	
	// interpolation
	if (i1 == i2) 			(*y) = contours[contour_index[i1] + contour_arr_size];
	else
	{
		double x1, x2;
		double y1, y2;
		
		x1 = contours[contour_index[i1]];
		x2 = contours[contour_index[i2]];
		
		y1 = contours[contour_index[i1] + contour_arr_size];
		y2 = contours[contour_index[i2] + contour_arr_size];
	
		printf("interpolation between (%f, %f) and (%f, %f)\n", x1, y1, x2, y2);

		if(y1 == y2) 		(*y) = y1;
		else
		{
			if (x1 == x2) 	(*y) = 0.5 * (y1 + y2); // simple mean 
			else 			(*y) = (y2 - y1) / (x2 - x1) * (x_value - x1) + y1; // linear interpolation
		}
	}
	
	printf("result : (%f, %f)\n", (*x), (*y));
}

void parameters_optimisation_2(double *X, double x_value, char **X_names, double* contours, double* coefficients, int contour_arr_size, bool maximize, double **x_ret, double **y_ret, int *ret_size)
{
	int i;
	int index;
	int ch_size, inf_size, sup_size;
	int *contour_sup, *contour_inf;
	
	double C_min, C;
	double *Yi, *Xi;
	double *X_mean;
	double *x_ch, *y_ch;

	Xi = (double*)malloc(2 * sizeof(double));
	Yi = (double*)malloc(2 * sizeof(double));
	X_mean = (double*)malloc(2 * sizeof(double));
	
	// mean coefficient values -> used to compute the evs cost function
	for(i = 0; i < NDIM; i++)
		X_mean[i] = X[i];
	
	ch_size = 0;
	add_to_convergence_history(X, &x_ch, &y_ch, &ch_size);
	
	splitContoursAlongY(contours, contour_arr_size, &contour_inf, &contour_sup, &inf_size, &sup_size);

	printf("Executing 2nd EVS optimisation script\n");
	
	// find 2 closest points on the left of the x value 
	get_contour_intersection(contours, contour_arr_size, x_value, contour_sup, sup_size, &(Xi[0]), &(Yi[0]));
	get_contour_intersection(contours, contour_arr_size, x_value, contour_inf, inf_size, &(Xi[1]), &(Yi[1]));
	
	for(i = 0; i < 2; i++)
	{
		X[0] = Xi[i];
		X[1] = Yi[i];
		
		printf("X[%d] : (%f, %f)\n", i, Xi[i], Yi[i]);
		
		add_to_convergence_history(X, &x_ch, &y_ch, &ch_size);
	}

	// optimisation
	printf("Running optimisation script\n");
	C_min = 1e40;
	index = 0;
	
	for(i = 0; i < 2; i++)
	{
		X[0] = Xi[i];
		X[1] = Yi[i];
		
		C = evs_cost_function(X, X_mean, X_names, coefficients, maximize);
		
		printf("C=%f / C_min=%f\n", C, C_min);
		
		if(C < C_min)
		{
			index = i;
			C_min = C;
		}
	}
	
	X[0] = Xi[index];
	X[1] = Yi[index];
	
	add_to_convergence_history(X, &x_ch, &y_ch, &ch_size);
	
	(*x_ret) = x_ch;
	(*y_ret) = y_ch;
	(*ret_size) = ch_size;
	
	printf("interpolation result : x=%f y=%f\n", X[0], X[1]);
	
	free(Xi);
	free(Yi);
	free(X_mean);
	
	for(i = 0; i < ch_size; i++)
		printf("i=%d : xch=%f, ych=%f\n", i, (*x_ret)[i], (*y_ret)[i]);
	
	printf("ok\n");
}

/*
*	Functiion used to find the closest point on the contour in the specified quadrant (+-x / +- y)
* 		Arguments : 
*			-	index				: (int) index of the initial point in the contour array
*			-	contour				: (double array) array containing contour points
*			-	contour_arr_size	: (int) size of the contour array
*			-	x_sign				: +1 or -1, used to set the quadrant in which the next point has to be (along x)
*			-	y_sign				: +1 or -1, used to set the quadrant in which the next point has to be (along y)
*
*		Return value : index of the next point 
*/ 
int findClosestPointAlongDirection(int index, double *contour, int contour_arr_size, int x_sign, int y_sign)
{
	int i, new_index;
	double x, y, min_dist, dist;
	
	new_index = -1;
	min_dist = 1e50;
	for(i = 0; i < contour_arr_size; i++)
	{
		if(i != index)
		{
			x = contour[i] - contour[index];
			y = contour[contour_arr_size + i] - contour[contour_arr_size + index];
			
			dist = sqrt(pow(x, 2) + pow(y, 2));
				
			if(dist < min_dist && (sgn(x) == sgn(x_sign) || x_sign == 0) && (sgn(y) == sgn(y_sign) || y_sign == 0))
			{
				min_dist = dist;
				new_index = i;
			}
		}
	}
		
	return new_index;
}

void add_value_to_array(int **array, int value, int *new_size)
{
	int *tmp_array;
	
	if((*new_size) == 0) // allocate new memory
	{
		(*new_size) = 1;
		tmp_array = (int*)malloc(sizeof(**array));
		
		if(tmp_array == NULL)
		{
			printf("ERROR : Could not allocate memory !\n");
			free(tmp_array);
			
			exit(0);
		}
	}
	else 
	{
		(*new_size) += 1;
		
		tmp_array = (int*)realloc(*array, (*new_size) * sizeof(**array));
		
		if(tmp_array == NULL)
		{
			printf("ERROR : Could not allocate memory !\n");
			
			free(tmp_array);			
			
			exit(0);
		}
	}
	
	tmp_array[(*new_size) - 1] = value;
	(*array) = tmp_array;
}

void splitContoursAlongY(double *contour, int contour_arr_size, int **_contour_inf, int **_contour_sup, int *_size_inf, int *_size_sup)
{
	double min, max;
	int min_index, max_index, i, index, new_index;
	
	int *contour_sup, *contour_inf;
	int size_inf, size_sup;
	
	size_inf = 0;
	size_sup = 0;
	
	contour_sup = (int*)malloc(0);
	contour_inf = (int*)malloc(0);
	
	max = -1e50;
	min = -max;
	
	// get min and max x contour value index
	for(i = 0; i < contour_arr_size; i++)
	{
		if(contour[i] < min)
		{
			min = contour[i];
			min_index = i;
		}
		else if(contour[i] > max)
		{
			max = contour[i];
			max_index = i;
		}
	}
	
	printf("min : %f, max : %f\n", min, max);
	
	// Loop inf contour
	add_value_to_array(&contour_inf, min_index, &size_inf);
	
	index = min_index;// start from min index
	
	LOOP_INF_CONTOUR:
	// find closest points on the +x direction and -y direction
	if (size_inf == 1) new_index = findClosestPointAlongDirection(index, contour, contour_arr_size, 1, -1);
	else new_index = findClosestPointAlongDirection(index, contour, contour_arr_size, 1, 0);
	
	if (new_index != -1)
	{
		index = new_index;
		
		// goto point and add index to inf contour
		add_value_to_array(&contour_inf, new_index, &size_inf);
			
		if(contour[new_index] < max) goto LOOP_INF_CONTOUR;
	}

	index = min_index;// start from min index
	
	// Loop sup contour
	LOOP_SUP_CONTOUR:
	// find closest points on the +x direction and +y direction
	if (size_sup == 0) new_index = findClosestPointAlongDirection(index, contour, contour_arr_size, 1, 1);
	else new_index = findClosestPointAlongDirection(index, contour, contour_arr_size, 1, 0);
	
	if (new_index != -1)
	{
		index = new_index;
		add_value_to_array(&contour_sup, new_index, &size_sup);

		if(contour[new_index] < max) goto LOOP_SUP_CONTOUR;
	}
	
	(*_contour_sup) = contour_sup;
	(*_contour_inf) = contour_inf;
	(*_size_inf) = size_inf;
	(*_size_sup) = size_sup;
}

void compute_initial_step(double *X, double *gradC, char **param_names, double* contours, int contour_arr_size, bool maximize)
{
	int i, j;
	int index;
	
	double dotp_min;
	double s1, s2;
	double dotp;
	
	// compute optimal initial step
	dotp_min = 1e16;
	
	for(i = 0; i < contour_arr_size; i++)
	{
		dotp = 0;
		
		s1 = 0;
		s2 = 0;
	
		for(j = 0; j < NDIM; j++)
		{
			dotp += (X[j] - contours[j * contour_arr_size + i]) * gradC[j];
			
			s1 += pow(X[j] - contours[j * contour_arr_size + i], 2);
			s2 += pow(gradC[j], 2);
		}
		
		dotp = dotp / (sqrt(s1) * sqrt(s2)); // calculate the normalized dot product of the two vectors
		
		if(dotp - 1 < dotp_min - 1 && dotp > 0.0)
		{
			dotp_min = dotp;
			index = i;
		}
	}
	
	X[0] = contours[index];
	X[1] = contours[contour_arr_size + index];
}

void add_to_convergence_history(double *X, double **x_ch, double **y_ch, int *ch_size)
{	
	double *x_ch_tmp, *y_ch_tmp;

	if((*ch_size) == 0) // allocate new memory
	{
		(*ch_size) = 1;

		x_ch_tmp = (double*)malloc(sizeof(**x_ch));
		y_ch_tmp = (double*)malloc(sizeof(**y_ch));
		
		if(x_ch_tmp == NULL || y_ch_tmp == NULL)
		{
			printf("ERROR : Could not allocate memory !\n");
			
			free(x_ch_tmp);
			free(y_ch_tmp);
			
			exit(0);
		}
	}
	else 
	{
		(*ch_size) += 1;
		
		x_ch_tmp = (double*)realloc(*x_ch, (*ch_size) * sizeof(**x_ch));
		y_ch_tmp = (double*)realloc(*y_ch, (*ch_size) * sizeof(**y_ch));
		
		if(x_ch_tmp == NULL || y_ch_tmp == NULL)
		{
			printf("ERROR : Could not allocate memory !\n");
			
			free(x_ch_tmp);
			free(y_ch_tmp);
			
			exit(0);
		}
	}
	
	x_ch_tmp[(*ch_size)-1] = X[0];
	y_ch_tmp[(*ch_size)-1] = X[1];
	
	(*x_ch) = x_ch_tmp;
	(*y_ch) = y_ch_tmp;
}

void clear_optimisation_cache_dat(double *X, double *Y)
{
	free(X);
	free(Y);
}

void clear_contour_cache_dat(int *contour_sup, int *contour_inf)
{
	free(contour_sup);
	free(contour_inf);
}