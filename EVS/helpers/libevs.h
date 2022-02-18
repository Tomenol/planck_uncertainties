#pragma once

#define _USE_MATH_DEFINES

#include <cmath>
#include <string>
#include <iostream>
#include <sstream> 

#define PI M_PI
#define c 299790.0

/**
*		FUNCTIONS RELATIVE TO MATH COMPUTATIONS
**/
template <class T> T max_fnc(T a, T b)
{
	return (a < b) ? b : a;     // or: return comp(a,b)?b:a; for version (2)
}

template <class T> std::string debugArray(T* arr, const int& len)
{
	std::stringstream ss;
	ss.precision(17);

	for(int i = 0; i < len; i++)
	{		
		ss << std::fixed << arr[i] << ' ';
	}

	return ss.str();
}

template <class T> std::string debugArray(T** arr, const int& len1, const int& len2)
{
	std::stringstream ss;

	for(int i = 0; i < len1; i++)
		for(int j = 0; j < len2; j++)
		{
			std::stringstream sstmp;
			sstmp.precision(17);
		
			sstmp << arr[i][j] << " ";
			ss << sstmp.str() << " ";
		}
	
	return ss.str();
}

template <typename T> int sgn(T val)
{
	return (int)((T(0) < val) - (val < T(0)));
}

#ifdef __cplusplus
extern "C"
{
#endif

	/**
	*		COSMOLOGY COMMONS
	*			Parameters and functions relative to the cosmology
	**/

	// LCDM parameters
	double OMEGA_M;
	double OMEGA_B;
	double SIGMA_8;
	double h;
	double ANS;
	double w_0;
	double w_a;

	//other cosmological parameters
	double RHO_M;
	double OMEGABH2;

	double OMEGA_R;

	double R_to_M;

	// observationnal parameters 
	double ZMAX, ZMIN, DZ;
	double DELTA_C;
	double MDELTA_MIN;
	double fsky;

	std::string select_mf_option;
	std::string select_delta_units;

	// cosmo functions
	void setCosmology(double Om, double Obh2, double h, double ns, double sigma8, double w_0, double w_a);
    void setObsParams(double _zmin, double _zmax, double _dz, double _delta_c, double _fsky, double _Mlim, char* _select_mf_option, char* _select_delta_units);

	double f_DE(double z);
	double inv_hubblez(double z);
	double dVdz(double z);
	double DELTAC(double OMEGAM, double zred);
	void COSMIC_VOLUME();
	double D_PLUS(double z);
	void GROWTH_FACTOR();
	double SELECTION(double M, double Z);

	/**
	*		POWER SPECTRUM DEFINITIONS
	**/
	const int NK = 1000;
	const int NM = 1000;
	const int NMRED = NM - 2;

	double M[NM], R[NM];
	double K_WAV[NK], DELTA2K[NK], D2DELTA2K[NK], TK_MATTER[NK];
	double SIGMA[NMRED], DSIGMADR[NMRED], RADIUS[NMRED], D2SIGMA[NMRED], D2DSIGMADR[NMRED];

	double hubble;

	// Power spectrum
	void MATTER_VARIANCE();
	double SIGMA_R(double R_Filter);
	double DSIGMADR_R(double R_Filter);

	// global variables power spectrum
	double theta_cmb;
	double z_equality;
	double k_equality;
	double z_drag;
	double R_drag;
	double R_equality;
	double sound_horizon;
	double k_silk;
	double alpha_c;
	double beta_c;
	double alpha_b;
	double beta_b;
	double beta_node;

	// other functions
	double DELTA2OFK(double K);
	void VARIANCE(double KMIN, double KMAX, int NRAD, double* RAD, double SIGMA8, double* SIG, double* DSIGDR);
	double WINDOW(double X);
	void TFfit(int NK, double OMEGA_M, double OMEGA_B, double h, double* K_WAV, double* TK_MATTER);
	void TFset_parameters(double omhh, double f_baryon, double Tcmb);
	void TFtransfer_function(double k, double omhh, double f_baryon, double* tf_full, double* tf_baryon, double* tf_cdm);
	double TF_pressureless(double q, double a, double b);
	double TF_zerobaryon(double q);
	double TF_nowiggles(double k, double omhh, double f_baryon, double Tcmb);
	double sound_horizon_fit(double omhh, double f_baryon);
	double k_peak(double omhh, double f_baryon);
	double alpha_gamma(double omhh, double f_baryon);

	/**
	*		MASS FUNCTIONS DEFINITIONS
	**/
	const int N_TINKER = 9;

	double deltab[N_TINKER], acap[N_TINKER], asmal[N_TINKER], b[N_TINKER], c_[N_TINKER];
	double d2acap[N_TINKER], d2asmal[N_TINKER], d2b[N_TINKER], d2c[N_TINKER];

	void TINKER_PARAMETER_DATA();
	double acap_zero_tink(double d);
	double asmal_zero_tink(double d);
	double b_zero_tink(double d);
	double c_zero_tink(double d);
	double F_TINKER(double Z, double SIGMAR, double DELTAB);
	double dNdM(double MASS, double REDSHIFT, double DELTAB);


	/**
	*		EVS INTEGRATION
	*			Parameters and functions relative to the integration of the evs formula and to processing
	**/
	const int NZ = 1000;
	double ZVOL[NZ], dVOL[NZ], d2dVOL[NZ], AFIN[NZ], DPLUS[NZ], D2DPLUS[NZ];
	
	void setup_evs_computation();
	
	void integration_over_zbins(int IZ_BIN, double &Z_BIN, double &NZ_CL, double LNM_SUP, double LNM_INF, double MASS_SUP, double MASS_INF, int NM_INT, double DELTA_MINT);
	void integration_over_Mmax(double& PDFHALOS_1, double& CUMULATIVEHALOS_1, double N_CL, double Z_BIN, double MASS_INF, double MASS_SUP, int NM_INT, double LNM_INF, double LNM_SUP, double DELTA_MINT);

	void findMassPeak(int NM_MAX, double* pdf, double* M, double &Mmax, double& pdfMax);
	void findMassConfIntervals(int NM_MAX, double* cumulDistrib, double* M, double &s1_p, double& s1_m, double &s2_p, double& s2_m);

	/**		DRIVER EVS
	*			Parameters and functions relative to math calculations
	**/
	double clever(double a, double b, double(*f)(double));
	double nsplint(double* xa, double* ya, double* y2a, int n, double x);
	void spline(double* x, double* y, int n, double yp1, double ypn, double* y2);
	void splie2(double* x1a, double* x2a, double** ya, int m, int n, double** y2a);
	double splin2(double* x1a, double* x2a, double** ya, double** y2a, int m, int n, double x1, double x2);
	void odeint(double* ystart, int nvar, double alfa1, double alfa2, double EPS, double hh1, double hmin, void(*derivs)(int, double, double*, double*), void(*rkqs)(double*, double*, int, double&, double, double, double*, double&, double&, void(*)(int, double, double*, double*)));
	void rkqs(double* y, double* dydx, int n, double &x, double htry, double EPS, double* yscal, double& hdid, double& hnext, void(*derivs)(int, double, double*, double*));
	void rkck(double* y, double* dydx, int n, double x, double h, double* yout, double* yerr, void(*derivs)(int, double, double*, double*));
	void fderivs(int n, double x, double* y, double* dydx);
#ifdef __cplusplus
}
#endif