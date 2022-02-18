#include "optimizer.h"

double boundary_cost_function(double *_X, double *f, double *grid, double level, double **bounds, int n, int m)
{
	return pow(interp2D(_X[0], _X[1], f, grid, bounds, n, m) - level, 2);
}

/* 
*		function that has to be minimize in order to maximize a given parameter x
*/
double f_minim(double x)
{
	return -x;
}

double evs_cost_function(double *X, char** param_names, bool maximize)
{
	int param_index;
	double cost_fnc, coeff;
	
	coeff = -1.0;
	cost_fnc = 0.0;
	
	if(maximize == false) coeff = 1.0;

	for (param_index = 0; param_index < NDIM; param_index++)
	{
		if(strcmp(param_names[param_index], "omegam") == 0 || strcmp(param_names[param_index], "omegabh2") == 0 || strcmp(param_names[param_index], "sigma8") == 0 || strcmp(param_names[param_index], "wa") == 0)
		{
			cost_fnc += f_minim(coeff * (X[param_index] - X_ch[0][param_index]));
		}
		else if (strcmp(param_names[param_index], "H0") == 0 || strcmp(param_names[param_index], "ns") == 0 || strcmp(param_names[param_index], "w") == 0)
		{
			cost_fnc += f_minim(-coeff * (X[param_index] - X_ch[0][param_index]));
		}
		else if (strcmp(param_names[param_index], "H0") == 0)
		{
			cost_fnc += f_minim(-coeff * (X[param_index] - X_ch[0][param_index]));
		}
		else
		{
			printf("Error : unknown param name : '%s' !\n", param_names[param_index]);
		}
	}
	
	return cost_fnc;
}

void evs_cost_function_grad_2(double *X, char** param_names, double *dX, bool maximize)
{
	int param_index;
	double coeff;
	
	coeff = -1.0;
	
	if(maximize == false) coeff = 1.0;

	for (param_index = 0; param_index < NDIM; param_index++)
	{
		if(strcmp(param_names[param_index], "omegam") == 0 || strcmp(param_names[param_index], "omegabh2") == 0 || strcmp(param_names[param_index], "sigma8") == 0 || strcmp(param_names[param_index], "wa") == 0)
		{
			dX[param_index] = coeff;
		}
		else if (strcmp(param_names[param_index], "ns") == 0 || strcmp(param_names[param_index], "w") == 0 || strcmp(param_names[param_index], "H0") == 0)
		{
			dX[param_index] = -coeff;
		}
		else
		{
			printf("Error : unknown param name : '%s' !\n", param_names[param_index]);
		}
	}
}

void evs_cost_function_grad_3(double *X, char** param_names, double *dX, double *coeffs, bool maximize)
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

void evs_cost_function_grad(double *X, char** param_names, double *dX, bool maximize)
{
	int param_index, param_index_tmp, i;
	
	double *X_tmp = (double*)malloc(NDIM * sizeof(double));
	int epsilon;
	
	for(param_index = 0; param_index < NDIM; param_index++)
	{
		dX[param_index] = 0;
		
		for(param_index_tmp = 0; param_index_tmp < NDIM; param_index_tmp++)
			X_tmp[param_index_tmp] = X[param_index_tmp];
		
		for(epsilon = 1; epsilon >= 0; epsilon--) // epsilon = 1 -> 0
		{
			X_tmp[param_index] = X[param_index] + (double)epsilon * DERIVATIVE_STEP_SIZE;
			dX[param_index] += (double)(2.0 * ((double)epsilon - 0.5)) * evs_cost_function(X_tmp, param_names, maximize);
			
			printf("%d : evs cost function : x=%f %f\n", param_index, X_tmp[param_index], evs_cost_function(X_tmp, param_names, maximize));
		}
		
		dX[param_index] = dX[param_index] / (DERIVATIVE_STEP_SIZE);
		
		printf("%d : grad evs cost function : %f\n", param_index, dX[param_index]);
	}
	
	free(X_tmp);
}

void cost_function_grad(double *X, double *dX, double *f, double *grid, double level, double **bounds, int n, int m)
{
	int param_index, param_index_tmp, i;
	
	double *X_tmp = (double*)malloc(NDIM * sizeof(double));
	int epsilon;
	
	for(param_index = 0; param_index < NDIM; param_index++)
	{
		dX[param_index] = 0;
		
		for(param_index_tmp = 0; param_index_tmp < NDIM; param_index_tmp++)
			X_tmp[param_index_tmp] = X[param_index_tmp];
				
		for(epsilon = 1; epsilon >= 0; epsilon--)
		{
			X_tmp[param_index] = X[param_index] + (double)epsilon * DERIVATIVE_STEP_SIZE;
			dX[param_index] += (double)(2.0 * ((double)epsilon - 0.5)) * boundary_cost_function(X_tmp, f, grid, level, bounds, n, m);
		}
		
		dX[param_index] = dX[param_index] / (DERIVATIVE_STEP_SIZE);
	}
	
	free(X_tmp);
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

void optimize_evs_gradient_descent(double *X, char *x_name, char *y_name, double *f, double *grid, double level, int n, int m, bool maximize)
{
	int evs_step, step, i;
	double cost_boundary;
	double *dX, ds;
	
	double **bounds = (double**)malloc(NDIM * sizeof(double*));
	for(i = 0; i < NDIM; i++)
		bounds[i] = (double*)malloc(2 * sizeof(double));
	
	get_bounds(grid, bounds, n, m);
	
	char *param_names[2] = {x_name, y_name};
	printf("X start : %f %f!\n", X[0], X[1]);
	
	dX = (double*)malloc(NDIM * sizeof(double));
	
	init_computation_cache_dat();
	add_to_convergence_history(X); // add initial position to history
	
	for(evs_step = 0; evs_step < MAX_EVS_STEPS; evs_step++)
	{
		evs_cost_function_grad_2(X, param_names, dX, maximize);

		for(i = 0; i < NDIM; i++)
			X[i] = X[i] - EVS_STEP_SIZE * dX[i];
		
		for(step = 0; step < MAX_STEPS; step++)
		{
			cost_function_grad(X, dX, f, grid, level, bounds, n, m);
			
			/*printf("X step : %f %f : (%f, %f) -> (%f, %f) !\n", 
			 -STEP_SIZE * dX[0], -STEP_SIZE * dX[1], X[0], X[1],
			X[0] - STEP_SIZE * dX[0], X[1] - STEP_SIZE * dX[1]);*/
			
			for(i = 0; i < NDIM; i++)
				X[i] = X[i] - STEP_SIZE * dX[i];

			cost_boundary = boundary_cost_function(X, f, grid, level, bounds, n, m);
			
			//printf("cost : %.10f -> pdf = %.10f  (x=%.5f, y=%.5f)\n", cost_boundary, interp2D(X[0], X[1], f, grid, bounds, n, m), X[0], X[1]);
			
			if (cost_boundary <= BOUNDARY_OPTIMIZATION_THRESHOLD)
			{
				printf("x = %f, y = %f\n", X[0], X[1]);
                break;
			}

			if (step == MAX_STEPS - 1)
			{
				printf("too much steps in optimisation script : %d\n", step);
                return;
			}
		}
		
		add_to_convergence_history(X);
		
		ds = 0;
		
		for(i = 0; i < NDIM; i++)
		{
			ds += pow(X_ch[i][X_ch_size - 1] - X_ch[i][X_ch_size - 3], 2);
		}
		
		if (sqrt(ds) <= BOUNDARY_OPTIMIZATION_THRESHOLD)
		{
			printf("x = %f, y = %f\n", X[0], X[1]);
			break;
		}
	}
	
	if (evs_step == MAX_EVS_STEPS - 1)
	{
		printf("too much steps in evs optimisation script : %d\n", evs_step);
		return;
	}
	
	printf("convergence : x=%f y=%f\n", X[0], X[1]);
	
	free(dX);
	
	release_2d_cache_addr(bounds, NDIM);
}

void optimize_evs_gradient_descent2(double *X, char *x_name, char *y_name, double *f, double *grid, double* contours, int n, int m, int contour_arr_size, bool maximize)
{
	int evs_step, i, j;
	
	double cost_boundary, ds, s, s_min;	
	double *dX, *X_tmp;
	
	char *param_names[2] = {x_name, y_name};
	
	double **bounds = (double**)malloc(NDIM * sizeof(double*));	
	for(i = 0; i < NDIM; i++)
		bounds[i] = (double*)malloc(2 * sizeof(double));
	
	get_bounds(grid, bounds, n, m);

	printf("X start : %f %f!\n", X[0], X[1]);
	
	dX = (double*)malloc(NDIM * sizeof(double));
	X_tmp = (double*)malloc(NDIM * sizeof(double));
	
	init_computation_cache_dat();
	
	add_to_convergence_history(X); // add initial position to history
	
	for(evs_step = 0; evs_step < MAX_EVS_STEPS; evs_step++)
	{
		evs_cost_function_grad(X, param_names, dX, maximize);

		for(i = 0; i < NDIM; i++)
			X[i] = X[i] - dX[i];
		
		printf("step %d : X = (%f, %f)\n", evs_step, X[0], X[1]);

		// get closest point on contour
		s_min = 1e16;
		
		for(j = 0; j < contour_arr_size; j++)// array of the form [[x0, x1], [y0, y1]]
		{
			s = 0;
			
			for(i = 0; i < NDIM; i++)
				s += pow(X[i] - contours[i * contour_arr_size + j], 2);
			
			s = sqrt(s);
			
			//printf("s=%f and s_min=%f\n", s, s_min);
			
			if(s < s_min)
			{
				s_min = s;
				
				//printf("s < s_min -> replacing point\n", s, s_min);
				
				for(i = 0; i < NDIM; i++)
					X_tmp[i] = contours[i * contour_arr_size + j];
			}
		}
		
		//printf("end -> X = (%f, %f)\n", X_tmp[0], X_tmp[1]);
		
		for(i = 0; i < NDIM; i++)
			X[i] = X_tmp[i];
		
		//printf("end -> X = (%f, %f)\n", X[0], X[1]);
		
		add_to_convergence_history(X);
		
		ds = 0;
		
		for(i = 0; i < NDIM; i++)
		{
			ds += pow(X_ch[i][X_ch_size - 1] - X_ch[i][X_ch_size - 2], 2);
		}
		
		if (sqrt(ds) <= 1e-3)
		{
			printf("x = %f, y = %f\n", X[0], X[1]);
			break;
		}
	}
	
	if (evs_step == MAX_EVS_STEPS - 1)
	{
		printf("too much steps in evs optimisation script : %d\n", evs_step);
		return;
	}
	
	printf("convergence : x=%f y=%f\n", X[0], X[1]);
	
	free(dX);
	free(X_tmp);
	
	release_2d_cache_addr(bounds, NDIM);
}

void optimize_evs_gradient_descent3(double *X, char *x_name, char *y_name, double *f, double *grid, double* contours, double* coefficients, int n, int m, int contour_arr_size, bool maximize, bool firststep)
{
	int i, j, k;
	int index, iflag;
	int i1, i2;
	
	double min_dist;
	double dotp1, dotp2, dist;
	double *gradC;
	
	char *param_names[2] = {x_name, y_name};
	
	printf("Executing 3rd EVS script\n");
	
	gradC = (double*)malloc(NDIM * sizeof(double));
	
	init_computation_cache_dat();
	
	evs_cost_function_grad_3(X, param_names, gradC, coefficients, maximize);

	if (firststep == true) 
	{
		compute_initial_step(X, gradC, param_names, contours, contour_arr_size, maximize);
		add_to_convergence_history(X);
	}

	CONVERGENCE_LOOP:
	
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
	
	printf("dot products : dotp1=%f / dotp2=%f\n", dotp1, dotp2);
	
	// test 
	if (dotp1 <= 0 && dotp2 <= 0)
	{
		printf("dead zone, exiting...\n");
		iflag = 2;
	}
	else if(dotp1 < 0 && dotp2 > 0)
	{
		printf("going towards M2.\n");
		
		index = i2;
		iflag = 1;
	}
	else if(dotp1 > 0 && dotp2 < 0)
	{
		printf("going towards M1.\n");
		
		index = i1;
		iflag = 1;
	}
	
	if(iflag == 1)
	{
		for(i = 0; i < NDIM; i++)
			X[i] = contours[i * contour_arr_size + index];
		
		add_to_convergence_history(X);
		
		goto CONVERGENCE_LOOP;
	}
	
	printf("interpolation result : x=%f y=%f\n", X[0], X[1]);
	
	free(gradC);
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

void add_to_convergence_history_index(int i)
{
	int *i_ch_tmp_i;
	
	i_ch_size += 1;
		
	i_ch_tmp_i = (int*)realloc(i_ch, i_ch_size * sizeof(*i_ch));
	
	if (i_ch == NULL) printf("Error : could not allocate memory for tmp buffer.\n");
	else
	{
		i_ch = i_ch_tmp_i;
		i_ch[i_ch_size - 1] = i;
	}
}

void add_to_convergence_history(double *X)
{
	int i;
	double *X_ch_tmp_i;
	
	X_ch_size += 1;
		
	for(i = 0; i < NDIM; i++)
	{		
		X_ch_tmp_i = (double*)realloc(X_ch[i], X_ch_size * sizeof(*X_ch[i]));
		
		if (X_ch[i] == NULL) printf("Error : could not allocate memory for tmp buffer.\n");
		else
		{
			X_ch[i] = X_ch_tmp_i;
			X_ch[i][X_ch_size - 1] = X[i];
		}
	}
}

void get_convergence_results(double *X, bool get_history)
{
	int i, j;
	
	if (get_history == false)
	{
		for(i = 0; i < NDIM; i++)
			X[i] = X_ch[i][X_ch_size - 1];
	}
	else if (get_history == true)
	{
		for(i = 0; i < NDIM; i++)
			for(j = 0; j < X_ch_size; j++)
				X[i * X_ch_size + j] = X_ch[i][j];
	}
}

int get_convergence_results_size()
{
	return X_ch_size;
}

void clear_computation_cache_dat()
{
	release_2d_cache_addr(X_ch, NDIM);
	free(i_ch);
	
	X_ch_size = 0;
	i_ch_size = 0;
}

void init_computation_cache_dat()
{
	int i;
	
	X_ch_size = 0;
	i_ch_size = 0;

	X_ch = (double**)malloc(NDIM * sizeof(double*));
	i_ch = (int*)malloc(0 * sizeof(int));

	for(i = 0; i < NDIM; i++)
		X_ch[i] = (double*)malloc(0 * sizeof(double));
}

void release_2d_cache_addr(double **ptr_arr, int n)
{
	int i;
	
	for(i = 0; i < n; i++)
		free(ptr_arr[i]);
	
	free(ptr_arr);
}










