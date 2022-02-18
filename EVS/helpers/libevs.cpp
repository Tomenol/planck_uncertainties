#include "libevs.h"

#ifdef __cplusplus
extern "C"
{
#endif
    /**     
    *   EVS INTEGRATION
    *			Parameters and functions relative to the integration of the evs formula
    **/
	void setup_evs_computation()
	{
		std::cout.precision(17);
		//std::cout << "starting setup for EVS computations" 	<< std::endl;
		
		TINKER_PARAMETER_DATA();
		COSMIC_VOLUME();
		GROWTH_FACTOR();
		MATTER_VARIANCE();
	}
	
	void findMassPeak(int NM_MAX, double* pdf, double* M, double &Mmax, double& pdfMax)
	{
		int i;
		double MPeak_tmp, max_PDF_tmp;
		
		Mmax = M[0];
		pdfMax = pdf[0];
		
		for (i = 1; i < NM_MAX; i++)
		{
			if (pdfMax < pdf[i])
			{
				pdfMax = pdf[i];
				Mmax = M[i];
			}
		}
	}
	
	void findMassConfIntervals(int NM_MAX, double* cumulDistrib, double* M, double &s1_p, double& s1_m, double &s2_p, double& s2_m)
	{
		int i;

		for (i = 0; i < NM_MAX; i++)
		{
          if((cumulDistrib[i] - 0.16) 	<= 1e-5) 	s1_m = M[i];
          if((cumulDistrib[i] - 0.84) 	<= 1e-5) 	s1_p = M[i];
          if((cumulDistrib[i] - 0.025) 	<= 1e-5) 	s2_m = M[i];
          if((cumulDistrib[i] - 0.975) 	<= 1e-5) 	s2_p = M[i];
		}
	}
	
    void setCosmology(double _Om, double _Obh2, double _h, double _ns, double _sigma8, double _w_0, double _w_a)
    {
	    OMEGA_M = _Om;
	    OMEGABH2 = _Obh2;
	    OMEGA_B = _Obh2 / pow(_h, 2);
	    OMEGA_R = 2.48 * 1e-5 / pow(_h, 2);

	    RHO_M = OMEGA_M * 2.78 * 1e11;

	    h = _h;
	    ANS = _ns;
	    SIGMA_8 = _sigma8;
		
		w_0 = _w_0;
		w_a = _w_a;

	    R_to_M = 3.0 / 4.0 / PI / RHO_M;
		
		/*std::cout.precision(17);
		
		std::cout << "STATUS : New cosmological parameters have been set : " << std::endl;
		std::cout << "OMEGA_M : " << OMEGA_M << std::endl;
		std::cout << "OMEGABH2 : " << OMEGABH2 << std::endl;
		std::cout << "OMEGA_B : " << OMEGA_B << std::endl;
		std::cout << "OMEGA_R : " << OMEGA_R << std::endl;
		std::cout << "RHO_M : " << RHO_M << std::endl;
		std::cout << "h : " << h << std::endl;
		std::cout << "NS : " << ANS << std::endl;
		std::cout << "SIGMA_8 : " << SIGMA_8 << std::endl;
		std::cout << "w_0 : " << w_0 << std::endl;
		std::cout << "w_a : " << w_a << std::endl;
		std::cout << "R_to_M : " << R_to_M << std::endl;*/
    }

    void setObsParams(double _zmin, double _zmax, double _dz, double _delta_c, double _fsky, double _Mlim, char* _select_mf_option, char* _select_delta_units)
    {
		std::stringstream ss;
		
	    ZMAX = (double) _zmax;
	    ZMIN = (double) _zmin;
	    DZ = (double) _dz;
	    DELTA_C = (double) _delta_c;
	    fsky = (double) _fsky;
		MDELTA_MIN = (double) _Mlim;
		
		ss << _select_mf_option << ' ' << _select_delta_units;
	    ss >> select_mf_option >> select_delta_units;
		
		/*std::cout.precision(17);
		
		std::cout << "STATUS : New observational and computational parameters have been set : " << std::endl;
		std::cout << "z_max : " << ZMAX << std::endl;
		std::cout << "z_min : " << ZMIN << std::endl;
		std::cout << "dZ : " << DZ << std::endl;
		std::cout << "Delta_C : " << DELTA_C << std::endl;
		std::cout << "f_sky : " << fsky << std::endl;
		std::cout << "Mlim : " << MDELTA_MIN << std::endl;
		std::cout << "Mass function : " << select_mf_option << std::endl;
		std::cout << "delta units : " << select_delta_units << std::endl;*/
    }

    double f_DE(double z)
    {
	    double a = 1.0 / (1.0 + z);
	    return exp(3.0 * w_a * (a - 1.0)) * pow((1.0 + z), (3.0 * (1.0 + w_0 + w_a)));
    }


    double inv_hubblez(double z)
    {
	    double hz = OMEGA_M * pow(1.0 + z, 3) + OMEGA_R * pow(1.0 + z, 4) + (1.0 - OMEGA_M - OMEGA_R) * f_DE(z);
	    return 1.0 / sqrt(hz);
    }


    double dVdz(double z)
    {
	    if (z == 0.0) return 0.0;
	    else return nsplint(ZVOL, dVOL, d2dVOL, NZ, z);
    }

    double DELTAC(double OMEGAM, double zred)
    {
	    double OMEGAMZ = OMEGAM * pow(1.0 + zred, 3) * pow(inv_hubblez(zred), 2);
	    return 3.0 / 20.0 * pow(12.0 * PI, (2.0 / 3.0)) * (1.0 + 0.0123 * log10(OMEGAMZ));
    }


    void COSMIC_VOLUME()
    {
	    int IZ;
	    double DIN, DEND;
	    double RDIST;
		
	    for (IZ = 0; IZ < NZ; IZ++)
	    {
		    ZVOL[IZ] = ZMIN + (ZMAX - ZMIN) * (double)(IZ) / (double)(NZ - 1);
		    RDIST = clever(0.0, ZVOL[IZ], inv_hubblez);
		    dVOL[IZ] = pow(RDIST, 2) * inv_hubblez(ZVOL[IZ]) * pow((c / (h * 100.)), 3);
	    }

	    DIN = 1e30;
	    DEND = 1e30;

	    spline(ZVOL, dVOL, NZ, DIN, DEND, d2dVOL);
    }

    double D_PLUS(double z)
    {
	    return nsplint(AFIN, DPLUS, D2DPLUS, NZ, 1.0 / (1.0 + z));
    }

    //-- - COMPUTE D + (z) for given cosmology------------------------------------

    void GROWTH_FACTOR()
    {
	    double EPS = 1e-7;
	    int nvar;
	    int IA;
	    double h1, hmin;
	    double XINI, AINI, XFIN, ZINI;
	    double y[2];
	    double ALFAFIN, ALFASTART, DIN, DEND;

	    double DPLUS_TEMP[NZ];

	    nvar = 2;
	    h1 = 1e-7;
	    hmin = 0.0;
	    ZINI = 100.0;
	    AINI = 1.0 / (ZINI + 1.0);
	    XINI = std::log(1.0 / (1.0 + ZINI));
	    XFIN = 0.0;

	    y[0] = std::log(3.0 / 5.0 * AINI);
	    y[1] = 0.6 * AINI / exp(y[0]);

	    for (IA = 0; IA < NZ; IA++)
	    {
		    ALFAFIN = XINI + (double)(IA + 1.0) * (XFIN - XINI) / (double)(NZ);
		    ALFASTART = XINI + (double)(IA) * (XFIN - XINI) / (double)(NZ);

		    odeint(y, nvar, ALFASTART, ALFAFIN, EPS, h1, hmin, fderivs, rkqs);
			
		    DPLUS_TEMP[IA] = exp(y[0]);
		    AFIN[IA] = exp(ALFAFIN);
	    }

	    for (IA = 0; IA < NZ; IA++)
		    DPLUS[IA] = DPLUS_TEMP[IA] / DPLUS_TEMP[NZ - 1];

	    DIN = 1e30;
	    DEND = 1e30;

	    spline(AFIN, DPLUS, NZ, DIN, DEND, D2DPLUS);
    }

    //-- - COMPUTE W(M, z) for given cosmologyand F_x---------------------------- -

    double SELECTION(double M, double Z)
    {
	    if (M >= MDELTA_MIN) return 1.0;
	    else return 0.0;
    }

    /**     EVS INTEGRATION
    *			Parameters and functions relative to the integration of the evs formula
    **/
    void MATTER_VARIANCE()
    {
        int IK, IM, JM;
        double DIN, DEND;
        double PKINI, PK;
        double LG10M_MIN, LG10M_MAX;

        TFfit(NK, OMEGA_M, OMEGA_B, h, K_WAV, TK_MATTER);

        for (IK = 0; IK < NK; IK++)
        {		
            PKINI = pow(K_WAV[IK] * h, (ANS - 1.0));
            PK = 2 * pow(PI, 2) * PKINI * (K_WAV[IK] * h) * pow(TK_MATTER[IK], 2) * pow(h, 3);
            DELTA2K[IK] = 1.0 / 2.0 / pow(PI, 2) * pow(K_WAV[IK], 3) * PK;
        }

        DIN = 1e30;
        DEND = 1e30;
		
        spline(K_WAV, DELTA2K, NK, DIN, DEND, D2DELTA2K);

        LG10M_MIN = 8.0;       //M_MIN = 10 ^ 8 M_sun / h;
        LG10M_MAX = 19.0;      //M_MAX = 10 ^ 19 M_sun / h;

        for (IM = 0; IM < NM; IM++)
        {
            M[IM] = pow(10, (LG10M_MIN + (LG10M_MAX - LG10M_MIN) * ((double)IM) / ((double)NM - 1.0)));
            R[IM] = pow((3.0 / 4.0 / PI * M[IM] / RHO_M), (1.0 / 3.0));   //!h ^ -1 Mpc;
        }

        VARIANCE(K_WAV[0], K_WAV[NK - 1], NM, R, SIGMA_8, SIGMA, DSIGMADR);

        for (JM = 1; JM < NMRED; JM++)
            RADIUS[JM - 1] = R[JM];

        spline(RADIUS, SIGMA, NMRED, DIN, DEND, D2SIGMA);
        spline(RADIUS, DSIGMADR, NMRED, DIN, DEND, D2DSIGMADR);
    }

    double SIGMA_R(double R_Filter)
    {
        return nsplint(RADIUS, SIGMA, D2SIGMA, NMRED, R_Filter);
    }


    double DSIGMADR_R(double R_Filter)
    {
        return nsplint(RADIUS, DSIGMADR, D2DSIGMADR, NMRED, R_Filter);
    }


    double DELTA2OFK(double K)
    {
		return nsplint(K_WAV, DELTA2K, D2DELTA2K, NK, K);
    }

    void VARIANCE(double KMIN, double KMAX, int NRAD, double* RAD, double SIGMA8, double* SIG, double* DSIGDR)
    {
        const int NKSTEP = 1000;

        double ANORM, ANORM_MIN, ANORM_MAX;
        double LN_KMIN, LN_KMAX, XMIN, XMAX, RADIUS, X, K;
        double DSIG2_MIN, DSIG2_MAX, DSIG2, DELTALNK;

        int IK, IR;

        double* SIG2;
        SIG2 = (double*)malloc(NRAD * sizeof(double));

        if (SIG2 == NULL)
        {
            //std::cout << "ERROR : memory for SIG2 could not be allocated !" << std::endl;
            exit(EXIT_FAILURE);
        }

        LN_KMIN = log(KMIN);
        LN_KMAX = log(KMAX);
        DELTALNK = (LN_KMAX - LN_KMIN) / (double)(NKSTEP);

        RADIUS = 8.0;
        XMIN = RADIUS * KMIN;
        XMAX = RADIUS * KMAX;

        ANORM_MIN = pow(WINDOW(XMIN), 2) * DELTA2OFK(KMIN);
        ANORM_MAX = pow(WINDOW(XMAX), 2) * DELTA2OFK(KMAX);
		
        ANORM = 0.0;

        for (IK = 0; IK < NKSTEP; IK++)
        {
            K = exp(LN_KMIN + (IK + 1) * DELTALNK);
            X = RADIUS * K;

            ANORM = ANORM + pow(WINDOW(X), 2) * DELTA2OFK(K);
        }

        ANORM = DELTALNK * (ANORM + (ANORM_MIN + ANORM_MAX) / 2.0);
        ANORM = pow(SIGMA8, 2) / ANORM;

        for (IR = 0; IR < NRAD; IR++)
        {
            RADIUS = RAD[IR];
            XMIN = RADIUS * KMIN;
            XMAX = RADIUS * KMAX;

            DSIG2_MIN = pow(WINDOW(XMIN), 2) * ANORM * DELTA2OFK(KMIN);
            DSIG2_MAX = pow(WINDOW(XMAX), 2) * ANORM * DELTA2OFK(KMAX);
            DSIG2 = 0.0;

            for (IK = 0; IK < NKSTEP; IK++)
            {
                K = exp(LN_KMIN + (IK + 1) * DELTALNK);
                X = RADIUS * K;
                DSIG2 = DSIG2 + pow(WINDOW(X), 2) * ANORM * DELTA2OFK(K);
            }

            SIG2[IR] = DELTALNK * (DSIG2 + (DSIG2_MIN + DSIG2_MAX) / 2.0);
        }

        for (IR = 1; IR < NRAD - 1; IR++)
        {
            SIG[IR - 1] = sqrt(SIG2[IR]);
            DSIGDR[IR - 1] = (sqrt(SIG2[IR + 1]) - sqrt(SIG2[IR - 1])) / 2.0 / (RAD[IR + 1] - RAD[IR - 1]);
        }

        free(SIG2);
    }

    double WINDOW(double X)
    {
        return 3.0 * (sin(X) - X * cos(X)) / pow(X, 3);
    }


    void TFfit(int NK, double OMEGA_M, double OMEGA_B, double _h, double* K_WAV, double* TK_MATTER)
    {
        double omega0, f_baryon, Tcmb, kmax, kmin;
        double k;
        double tf_full, tf_baryon, tf_cdm;
        double omhh;
        int i, numk;

        omega0 = (double)(OMEGA_M);
        f_baryon = (double)(OMEGA_B) / (double)(OMEGA_M);
        hubble = (double)(_h);
        Tcmb = 2.7255;
        omhh = omega0 * hubble * hubble;

        TFset_parameters(omhh, f_baryon, Tcmb);

        kmax = 10.0;
        numk = NK;

        kmin = 0.0001;

        for (i = 0; i < numk; i++)
        {
            k = pow(10.0, ((i + 1) * (log10(kmax / kmin) / numk))) * kmin;

            TFtransfer_function(k * hubble, omhh, f_baryon, &tf_full, &tf_baryon, &tf_cdm);

            K_WAV[i] = (double)(k);
            TK_MATTER[i] = tf_full;
        }
    }

    void TFset_parameters(double omhh, double f_baryon, double Tcmb)
    {
        double obhh, y;

        hubble = h;

        if (f_baryon <= 0.0) f_baryon = 1e-5;
        if (Tcmb <= 0.0) Tcmb = 2.728;
        if (hubble > 10.0)
            std::cout << "TFset_parameters(): WARNING, Hubble constant in 100km/s/Mpc desired" << std::endl;

        obhh = omhh * f_baryon;
        theta_cmb = Tcmb / 2.7;

        z_equality = 2.50e4 * omhh * pow(theta_cmb, -4.) - 1.0;
        k_equality = 0.0746 * omhh * pow(theta_cmb, -2.);

        z_drag = 0.313 * pow(omhh, -0.419) * (1.0 + 0.607 * pow(omhh, 0.674));
        z_drag = 1e0 + z_drag * pow(obhh, (0.238 * pow(omhh, 0.223)));
        z_drag = 1291.0 * pow(omhh, 0.251) / (1e0 + 0.659 * pow(omhh, 0.828)) * z_drag;

        R_drag = 31.5 * obhh * pow(theta_cmb, -4.) * 1000e0 / (1e0 + z_drag);
        R_equality = 31.5 * obhh * pow(theta_cmb, -4.) * 1000e0 / (1e0 + z_equality);

        sound_horizon = 2. / 3. / k_equality * sqrt(6. / R_equality) * log((sqrt(1. + R_drag) + sqrt(R_drag + R_equality)) / (1. + sqrt(R_equality)));

        k_silk = 1.6 * pow(obhh, 0.52) * pow(omhh, 0.73) * (1e0 + pow(10.4 * omhh, -0.95));

        alpha_c = pow((46.9 * omhh), 0.670) * (1e0 + pow((32.1 * omhh), -0.532));
        alpha_c = pow(alpha_c, -f_baryon);
        alpha_c = alpha_c * pow((pow(12.0 * omhh, 0.424) * (1.0 + pow(45.0 * omhh, -0.582))), -pow(f_baryon, 3.0));

        beta_c = 0.944 / (1 + pow(458. * omhh, -0.708));
        beta_c = 1.0 + beta_c * (pow(1.0 - f_baryon, pow(0.395 * omhh, -0.0266)) - 1.0);
        beta_c = 1.0 / beta_c;

        y = (1e0 + z_equality) / (1e0 + z_drag);
        alpha_b = y * (-6. * sqrt(1. + y) + (2. + 3. * y) * log((sqrt(1. + y) + 1.) / (sqrt(1. + y) - 1.)));
        alpha_b = 2.07 * k_equality * sound_horizon * pow(1. + R_drag, -0.75) * alpha_b;


        beta_b = 0.5 + f_baryon + (3. - 2. * f_baryon) * sqrt(pow(17.2 * omhh, 2.0) + 1e0);

        beta_node = 8.41 * pow(omhh, 0.435);
    }

    void TFtransfer_function(double k, double omhh, double f_baryon, double* tf_full, double* tf_baryon, double* tf_cdm)
    {
        double q, ks;
        double s_tilde;

        q = k / 13.41 / k_equality;
        ks = k * sound_horizon;

        (*tf_cdm) = 1.0 / (1.0 + pow(ks / 5.4, 4.0));
        (*tf_cdm) = (*tf_cdm) * TF_pressureless(q, 1.0, beta_c) + (1.0 - (*tf_cdm)) * TF_pressureless(q, alpha_c, beta_c);


        s_tilde = sound_horizon / pow(1. + pow(beta_node / ks, 3.0), 1.0 / 3.0);
		
        (*tf_baryon) = TF_pressureless(q, 1.0, 1.0) / (1.0 + pow(ks / 5.2, 2.0));		
        (*tf_baryon) = (*tf_baryon) + alpha_b / (1. + pow(beta_b / ks, 3)) * exp(-pow(k / k_silk, 1.4));
        (*tf_baryon) = (*tf_baryon) * (sin(k * s_tilde) / (k * s_tilde));		
        (*tf_full) = f_baryon * (*tf_baryon) + (1 - f_baryon) * (*tf_cdm);
    }

    double TF_pressureless(double q, double a, double b)
    {
        double _TF_pressureless;

        _TF_pressureless = log(exp(1.0) + 1.8 * b * q);
        return _TF_pressureless / (_TF_pressureless + (14.2 / a + 386 / (1. + 69.9 * pow(q, 1.08))) * pow(q, 2));		
    }

    double TF_zerobaryon(double q)
    {
        double _TF_zerobaryon;

        _TF_zerobaryon = log(2.0 * exp(1.) + 1.8 * q);
        return _TF_zerobaryon / (_TF_zerobaryon + (14.2 + 731.0 / (1 + 62.5 * q)) * pow(q, 2));
    }

    double TF_nowiggles(double k, double omhh, double f_baryon, double Tcmb)
    {
        double q_eff, a;

        if (Tcmb <= 0.0) Tcmb = 2.728;
        a = alpha_gamma(omhh, f_baryon);
        q_eff = k / omhh * pow(Tcmb / 2.7, 2);
        q_eff = q_eff / (a + (1. - a) / (1. + pow(0.43 * k * sound_horizon_fit(omhh, f_baryon), 4)));

        return TF_zerobaryon(q_eff);
    }

    double sound_horizon_fit(double omhh, double f_baryon)
    {
        return 44.5 * log(9.83 / omhh) / sqrt(1.0 + 10.0 * pow(f_baryon * omhh, 0.75));
    }

    double k_peak(double omhh, double f_baryon)
    {
        return 5. * 3.14159 / 2. * (1. + 0.217 * omhh) / sound_horizon_fit(omhh, f_baryon);
    }

    double alpha_gamma(double omhh, double f_baryon)
    {
        return 1. - 0.328 * log(431.0 * omhh) * f_baryon + 0.38 * log(22.3 * omhh) * pow(f_baryon, 2);
    }

    /**
    *		MASS FUNCTIONS DEFINITIONS
    **/
    void TINKER_PARAMETER_DATA()
    {
        double din, den;

        deltab[0] = 200;
        deltab[1] = 300;
        deltab[2] = 400;
        deltab[3] = 600;
        deltab[4] = 800;
        deltab[5] = 1200;
        deltab[6] = 1600;
        deltab[7] = 2400;
        deltab[8] = 3200;

        acap[0] = 0.186;
        acap[1] = 0.200;
        acap[2] = 0.212;
        acap[3] = 0.218;
        acap[4] = 0.248;
        acap[5] = 0.255;
        acap[6] = 0.260;
        acap[7] = 0.260;
        acap[8] = 0.260;

        din = 1e30;
        den = 1e30;
		
        spline(deltab, acap, N_TINKER, din, den, d2acap);

        asmal[0] = 1.47;
        asmal[1] = 1.52;
        asmal[2] = 1.56;
        asmal[3] = 1.61;
        asmal[4] = 1.87;
        asmal[5] = 2.13;
        asmal[6] = 2.30;
        asmal[7] = 2.53;
        asmal[8] = 2.66;

        spline(deltab, asmal, N_TINKER, din, den, d2asmal);

        b[0] = 2.57;
        b[1] = 2.25;
        b[2] = 2.05;
        b[3] = 1.87;
        b[4] = 1.59;
        b[5] = 1.51;
        b[6] = 1.46;
        b[7] = 1.44;
        b[8] = 1.41;

        spline(deltab, b, N_TINKER, din, den, d2b);

        c_[0] = 1.19;
        c_[1] = 1.27;
        c_[2] = 1.34;
        c_[3] = 1.45;
        c_[4] = 1.58;
        c_[5] = 1.80;
        c_[6] = 1.97;
        c_[7] = 2.23;
        c_[8] = 2.44;

        spline(deltab, c_, N_TINKER, din, den, d2c);
    }

    double acap_zero_tink(double d)
    {
        return nsplint(deltab, acap, d2acap, N_TINKER, d);
    }

    double asmal_zero_tink(double d)
    {
        return nsplint(deltab, asmal, d2asmal, N_TINKER, d);
    }

    double b_zero_tink(double d)
    {
        return nsplint(deltab, b, d2b, N_TINKER, d);
    }

    double c_zero_tink(double d)
    {
        return nsplint(deltab, c_, d2c, N_TINKER, d);
    }

    double F_TINKER(double Z, double SIGMAR, double DELTAB)
    {
        double ACAP, ASMAL, BSMAL, CSMAL, ALPHA;
        double AAZ, AZ, BZ;

        ACAP = acap_zero_tink(DELTAB);
        ASMAL = asmal_zero_tink(DELTAB);
        BSMAL = b_zero_tink(DELTAB);
        CSMAL = c_zero_tink(DELTAB);

        ALPHA = pow(10.0, -pow(0.75 / log10(DELTAB / 75.), 1.2));

        AAZ = ACAP * pow(1.0 + Z, -0.14);
        AZ = ASMAL * pow(1.0 + Z, -0.06);
        BZ = BSMAL * pow(1.0 + Z, -ALPHA);

        return AAZ * (pow(SIGMAR / BZ, -AZ) + 1.0) * exp(-CSMAL / pow(SIGMAR, 2));
    }

    double dNdM(double MASS, double REDSHIFT, double DELTAB)
    {
        double DELTA, DELTA_VIR, x, nu;
        double OMEGAMZ, RSIZE, DSIGMADM;
        double ACAP, ASMALL, ALPHA, F_DESPALI;

        RSIZE = pow(R_to_M * MASS, 1.0 / 3.0);
        DSIGMADM = DSIGMADR_R(RSIZE) * R_to_M / 3. / pow(R_to_M * MASS, 2.0 / 3.0);
        nu = pow(DELTAC(OMEGA_M, REDSHIFT) / SIGMA_R(RSIZE) / D_PLUS(REDSHIFT), 2);

        OMEGAMZ = OMEGA_M * pow(1.0 + REDSHIFT, 3) * pow(inv_hubblez(REDSHIFT), 2);

        if (select_delta_units == "rho_c")
        {
            DELTA = DELTAB / OMEGAMZ;
        }
        else
        {
            DELTA = DELTAB;
        }

        if (select_mf_option == "Tinker")
        {
            return F_TINKER(REDSHIFT, SIGMA_R(RSIZE) * D_PLUS(REDSHIFT), DELTA) * RHO_M / MASS * (-1.0 / SIGMA_R(RSIZE) * DSIGMADM);
        }
        else
        {
            DELTA_VIR = 18.0 * pow(PI, 2) + 82.0 * (OMEGAMZ - 1.0) - 39.0 * pow(OMEGAMZ - 1.0, 2);

            x = log10(DELTA / DELTA_VIR);

            ACAP = -0.1362 * x + 0.3292;
            ASMALL = 0.4332 * pow(x, 2) + 0.2263 * x + 0.7665;
            ALPHA = -0.1151 * pow(x, 2) + 0.2554 * x + 0.2488;

            F_DESPALI = 2.0 * ACAP * sqrt(2.0 * ASMALL / PI * nu) * (1.0 + pow(1.0 / (nu * ASMALL), ALPHA)) * exp(-ASMALL * nu / 2.0);

            return F_DESPALI * RHO_M / MASS * (-1.0 / SIGMA_R(RSIZE) * DSIGMADM);
        }
    }
	
	void integration_over_zbins(int IZ_BIN, double &Z_BIN, double &NZ_CL, double LNM_SUP, double LNM_INF, double MASS_SUP, double MASS_INF, int NM_INT, double DELTA_MINT)
	{
		double DELTAB_ZBIN, DNDMDZ_INF, DNDMDZ_SUP, DNDMDZ, MASS_TEMP, DNDZ;
		int IM_INT;
		
		Z_BIN = ZMIN + (IZ_BIN + 0.50) * DZ;
		
		DELTAB_ZBIN = DELTA_C/(OMEGA_M * pow(1.0 + Z_BIN, 3) * pow(inv_hubblez(Z_BIN), 2));

		DNDMDZ_INF = dNdM(MASS_INF, Z_BIN, DELTAB_ZBIN) * MASS_INF * SELECTION(MASS_INF, Z_BIN);
		DNDMDZ_SUP = dNdM(MASS_SUP, Z_BIN, DELTAB_ZBIN) * MASS_SUP * SELECTION(MASS_SUP, Z_BIN);
		DNDMDZ = 0.0;

		for(IM_INT = 1; IM_INT < NM_INT - 1; IM_INT++)
		{
			MASS_TEMP = exp(LNM_INF + (LNM_SUP - LNM_INF) * (double)(IM_INT)/(double)(NM_INT - 1));
			DNDMDZ = DNDMDZ + dNdM(MASS_TEMP, Z_BIN, DELTAB_ZBIN) * MASS_TEMP * SELECTION(MASS_TEMP, Z_BIN);
		}

		DNDZ = (DNDMDZ + (DNDMDZ_INF + DNDMDZ_SUP)/2.0) * DELTA_MINT;

		NZ_CL = DNDZ * fsky * dVdz(Z_BIN) * (double)DZ;
		
		//std::cout << "integration over bin " << IZ_BIN << " done ! result : " << NZ_CL << std::endl;
	}
	
	void integration_over_Mmax(double& PDFHALOS_1, double& CUMULATIVEHALOS_1, double N_CL, double Z_BIN, double MASS_INF, double MASS_SUP, int NM_INT, double LNM_INF, double LNM_SUP, double DELTA_MINT)
	{
		double DELTAB_ZBIN, DFDMDZ_INF, DFDMDZ_SUP, DFDMDZ, MASS_TEMP, DFDZ, DFDM, FMZ;
		int IM_INT;
		
		DELTAB_ZBIN = DELTA_C/(OMEGA_M*pow(1.0 + Z_BIN, 3)*pow(inv_hubblez(Z_BIN), 2.0));

		DFDMDZ_INF = dNdM(MASS_INF, Z_BIN, DELTAB_ZBIN) * MASS_INF * SELECTION(MASS_INF, Z_BIN);
		DFDMDZ_SUP = dNdM(MASS_SUP, Z_BIN, DELTAB_ZBIN) * MASS_SUP * SELECTION(MASS_SUP, Z_BIN);
		DFDMDZ = 0.0;

		for(IM_INT = 1; IM_INT < NM_INT-1; IM_INT++)
		{
			MASS_TEMP = exp(LNM_INF + (LNM_SUP - LNM_INF) * (double)(IM_INT)/(double)(NM_INT - 1));
			DFDMDZ = DFDMDZ + dNdM(MASS_TEMP, Z_BIN, DELTAB_ZBIN) * MASS_TEMP * SELECTION(MASS_TEMP, Z_BIN);
		}

		DFDZ = (DFDMDZ + (DFDMDZ_INF + DFDMDZ_SUP)/2.0)*DELTA_MINT;
		FMZ =  DFDZ * fsky * dVdz(Z_BIN) * DZ/N_CL;

		DFDM = fsky * dVdz(Z_BIN) * DZ/N_CL * dNdM(MASS_SUP, Z_BIN, DELTAB_ZBIN) * SELECTION(MASS_SUP, Z_BIN);
 
		PDFHALOS_1 = MASS_SUP * (N_CL * pow(FMZ, N_CL - 1.0) * DFDM);
		CUMULATIVEHALOS_1 = pow(FMZ, N_CL);
	}

    /**
    *		DRIVER EVS
    *			Parameters and functions relative to math calculations
    **/
    double clever(double a, double b, double (*f)(double))
    {
        int LORR[200], IFLAG, LVL;

        double FV[5], FIT[200], F2T[200], F3T[200], DAT[200];
        double ARESTT[200], ESTT[200], EPST[200], PSUM[200];

        double U, ACC, FOURU, EPS, _ERROR, _ALPHA_param, _DA, AREA, AREST, KOUNT;
        double WT, EST, DX, ESTL, ESTR, SUM;
        double ARESTL, ARESTR, DIFF;

        U = 1.0E-7;
        ACC = 1.0E-4;
        FOURU = 4.0 * U;
        IFLAG = 1;
        EPS = ACC;
        _ERROR = 0.0;

        LVL = 0;
        LORR[LVL] = 1;
        PSUM[LVL] = 0.0;

        _ALPHA_param = a;
        _DA = b - a;
        AREA = 0.0;
        AREST = 0.0;

        FV[0] = f(_ALPHA_param);
        FV[2] = f(_ALPHA_param + 0.5 * _DA);
        FV[4] = f(_ALPHA_param + _DA);
        KOUNT = 3;
        WT = _DA / 6.0;
        EST = WT * (FV[0] + 4.0 * FV[2] + FV[4]);

    goto1:
        DX = 0.5 * _DA;
        FV[1] = f(_ALPHA_param + 0.5 * DX);
        FV[3] = f(_ALPHA_param + 1.5 * DX);

        KOUNT = KOUNT + 2;

        WT = DX / 6.0;
        ESTL = WT * (FV[0] + 4.0 * FV[1] + FV[2]);
        ESTR = WT * (FV[2] + 4.0 * FV[3] + FV[4]);
        SUM = ESTL + ESTR;

        //Added by Max :
        if (!((SUM >= 0) || (SUM <= 0))) return SUM; //Sum = NaN, so might as well bail out right now to save time

        ARESTL = WT * (std::abs(FV[0]) + std::abs(4. * FV[1]) + std::abs(FV[2]));
        ARESTR = WT * (std::abs(FV[2]) + std::abs(4. * FV[3]) + std::abs(FV[4]));
        AREA = AREA + ((ARESTL + ARESTR) - AREST);
        DIFF = EST - SUM;

        if (std::abs(DIFF) <= EPS * std::abs(AREA)) goto goto2;
        if (std::abs(DX) <= FOURU * std::abs(_ALPHA_param)) goto goto51;
        if (LVL >= 199) goto goto5;
        if (KOUNT >= 1999) goto goto6;

        LVL = LVL + 1;
        LORR[LVL] = 0;
        FIT[LVL] = FV[2];
        F2T[LVL] = FV[3];
        F3T[LVL] = FV[4];
        _DA = DX;
        DAT[LVL] = DX;
        AREST = ARESTL;
        ARESTT[LVL] = ARESTR;
        EST = ESTL;
        ESTT[LVL] = ESTR;
        EPS = EPS / 1.4;
        EPST[LVL] = EPS;
        FV[4] = FV[2];
        FV[2] = FV[1];
        goto goto1;

    goto2:
        _ERROR = _ERROR + DIFF / 15.;

    goto3:
        if (LORR[LVL] == 0) goto goto4;
        SUM = PSUM[LVL] + SUM;
        LVL = LVL - 1;

        if (LVL > 0) goto goto3;
        return SUM;

    goto4:
        PSUM[LVL] = SUM;
        LORR[LVL] = 1;
        _ALPHA_param = _ALPHA_param + _DA;
        _DA = DAT[LVL];
        FV[0] = FIT[LVL];
        FV[2] = F2T[LVL];
        FV[4] = F3T[LVL];
        AREST = ARESTT[LVL];
        EST = ESTT[LVL];
        EPS = EPST[LVL];
        goto goto1;

    goto5:
        IFLAG = 2;
        goto goto2;
    goto6:
        IFLAG = 3;
        goto goto2;
    goto51:
        IFLAG = 4;
        goto goto2;
    }

    double nsplint(double* xa, double* ya, double* y2a, int n, double x)
    {
        /**
        *   USE nrtype
        *       Given the arrays xa(1:n) and ya(1:n) of length n, which tabulate a function
        *       (with the xa(i) in order), and given the array y2a(1:n), which is the output
        *       from the subroutine spline, and given a value of x, this routine returns a
        *       cubic spline interpolated value y.
        *       (adopted from Numerical Recipes in FORTRAN 77)
        **/

        int k, khi, klo;
        double a, b, _h;

        klo = 0;
        khi = n-1;

    goto1:
        if (khi - klo > 1)
        {
            k = (khi + klo) / 2;

            if (xa[k] > x)
                khi = k;
            else
                klo = k;
            goto goto1;
        }

        _h = xa[khi] - xa[klo];

        if (_h == 0.0) return -1; //'bad xa input in splint'

        a = (xa[khi] - x) / _h;
        b = (x - xa[klo]) / _h;
		
        return a * ya[klo] + b * ya[khi] + ((pow(a, 3) - a) * y2a[klo] + (pow(b, 3) - b) * y2a[khi]) * pow(_h, 2) / 6.0;
    }

    void spline(double* x, double* y, int n, double yp1, double ypn, double* y2)
    {
        /** Given arrays x(1:n) and y(1:n) containing a tabulated function, i.e.
        *       y(i) = f(x(i)), with x(1) < x(2) < ... < x(n), and given values yp1and ypn for
        *       the first derivative of the interpolating function at points 1 and n,
        *       respectively, this routine returns an array y2(1:n) of length n which
        *       contains the second derivatives of the interpolating function at the
        *       tabulated points x(i).If yp1and /or ypn are equal to 1.e30 or larger,
        *       the routine is signaled to set the corresponding boundary condition for a
        *       natural spline with zero second derivative on that boundary.
        *       Parameter: nmax is the largest anticipiated value of n
        *       (adopted from Numerical Recipes in FORTRAN 77)
        **/

        int i, k;
        double p, qn, sig, un;

        double* u;

        u = (double*)malloc(n * sizeof(double));
		if (u == NULL) std::cout << "could not allocate memory for u (size " << n << ")" << std::endl;
		
        if (yp1 >= 1e30)
        {
            y2[0] = 0.0;
            u[0] = 0.0;
        }
        else
        {
            y2[0] = -0.5;
            u[0] = (3. / (x[1] - x[0])) * ((y[1] - y[0]) / (x[1] - x[0]) - yp1);
        }
        for (i = 1; i < n - 1; i++)
        {
            sig = (x[i] - x[i - 1]) / (x[i + 1] - x[i - 1]);
            p = sig * y2[i - 1] + 2.0;
            y2[i] = (sig - 1.) / p;
            u[i] = (6. * ((y[i + 1] - y[i]) / (x[i + 1] - x[i]) - (y[i] - y[i - 1]) / (x[i] - x[i - 1])) / (x[i + 1] - x[i - 1]) - sig * u[i - 1]) / p;
        }

        if (ypn >= 1e30)
        {
            qn = 0.0;
            un = 0.0;
        }
        else
        {
            qn = 0.5;
			un = (3.0 / (x[n-1] - x[n-2])) * (ypn - (y[n-1] - y[n - 2]) / (x[n-1] - x[n - 2]));
        }

        y2[n - 1] = (un - qn * u[n - 2]) / (qn * y2[n - 2] + 1.0);

        for (k = n - 2; k >= 0; k--)
            y2[k] = y2[k] * y2[k + 1] + u[k];

        free(u);
    }

    /** sline for a 2d function
        x1a(m), x2a(n), y2a(m, n), ya(m, n)
    **/
    void splie2(double* x1a, double* x2a, double** ya, int m, int n, double** y2a)
    {
        int j, k;
        double* y2tmp, * ytmp;

        y2tmp = (double*)malloc(n * m * sizeof(double));
        ytmp = (double*)malloc(n * m * sizeof(double));

        for (j = 0; j < m; j++)
        {
            for (k = 0; k < n; k++)
                ytmp[k] = ya[j][k];

            spline(x2a, ytmp, n, 1e30, 1e30, y2tmp);

            for (k = 0; k < n; k++)
                y2a[j][k] = y2tmp[k];
        }

        free(y2tmp);
        free(ytmp);
    }

    /*
        int m, n
        REAL(8) ::x1, x2, y, x1a(m), x2a(n), y2a(m, n), ya(m, n)
    */

    double splin2(double* x1a, double* x2a, double** ya, double** y2a, int m, int n, double x1, double x2)
    {
        int j, k;
        double* y2tmp, * ytmp, * yytmp;
        double res_splin2;

        y2tmp = (double*)malloc(n * m * sizeof(double));
        ytmp = (double*)malloc(n * m * sizeof(double));
        yytmp = (double*)malloc(n * sizeof(double));

        for (j = 0; j < m; j++)
        {
            for (k = 0; k < n; k++)
            {
                ytmp[k] = ya[j][k];
                y2tmp[k] = y2a[j][k];
            }

            yytmp[j] = nsplint(x2a, ytmp, y2tmp, n, x2);
        }

        spline(x1a, yytmp, m, 1e30, 1e30, y2tmp);

        res_splin2 = nsplint(x1a, yytmp, y2tmp, m, x1);

        free(y2tmp);
        free(yytmp);
        free(ytmp);

        return res_splin2;
    }

    /*
        integer nbad, nok, nvar, KMAXX, MAXSTP, NMAX, kmax, kount
        real * 8 EPS, hh1, hmin, ystart(nvar), TINY, alfa1, alfa2
        external derivs, rkqs
    */
    void odeint(double* ystart, int nvar, double alfa1, double alfa2, double EPS, double hh1, double hmin, void(*derivs)(int, double, double*, double*), void(*rkqs)(double*, double*, int, double&, double, double, double*, double&, double&, void(*)(int, double, double*, double*)))
    {
        const double TINY = 1e-30;
        const int MAXSTP = 1000000;

        int i, nstp, nok, nbad;
        double h, hdid, hnext, x, * dydx, * y, * yscal;

        dydx = (double*)malloc(nvar * sizeof(double));
        y = (double*)malloc(nvar * sizeof(double));
        yscal = (double*)malloc(nvar * sizeof(double));

        x = alfa1;
		
        h = std::abs((double) hh1) * sgn(alfa2 - alfa1);
		
        nok = 0;
        nbad = 0;

        hdid = 0;
        hnext = 0;

        for (i = 0; i < nvar; i++)
            y[i] = ystart[i];
		
        for (nstp = 0; nstp < MAXSTP; nstp++)
        {
            derivs(nvar, x, y, dydx);

            for (i = 0; i < nvar; i++)
                yscal[i] = std::abs(y[i]) + std::abs(h * dydx[i]) + TINY;

            if ((x + h - alfa2) * (x + h - alfa1) > 0.0) 
				h = alfa2 - x;

            rkqs(y, dydx, nvar, x, h, EPS, yscal, hdid, hnext, derivs);

            if (hdid == h)
                nok = nok + 1;
            else
                nbad = nbad + 1;

            if ((x - alfa2) * (alfa2 - alfa1) >= 0.0)
            {
                for (i = 0; i < nvar; i++)
                    ystart[i] = y[i];
				
                goto goto_return;
            }

            h = hnext;

            if (std::abs(hnext) < hmin)
                h = hmin;
        }

        std::cout << "too many steps in odeint" << std::endl;
	
	goto_return:
		free(dydx);
		free(y);
		free(yscal);
    }

    void rkqs(double* y, double* dydx, int n, double &x, double htry, double EPS, double* yscal, double& hdid, double& hnext, void(*derivs)(int, double, double*, double*))
    {
        int i;

        double errmax, h, xnew;
        double* yerr, * ytemp;

        yerr = (double*)malloc(n * sizeof(double));
        ytemp = (double*)malloc(n * sizeof(double));

        const double SAFETY = 0.9;
        const double PGROW = -0.2;
        const double PSHRNK = -0.25;
        const double ERRCON = 1.89e-4;

        h = htry;

    goto1:
        rkck(y, dydx, n, x, h, ytemp, yerr, derivs);
        errmax = 0.0;
        for (i = 0; i < n; i++)
            errmax = max_fnc((double) errmax, (double) std::abs(yerr[i] / yscal[i]));

        errmax = errmax / EPS;
        if (errmax > 1.0)
        {
            h = SAFETY * h * (pow(errmax, PSHRNK));

            if (h < 0.1 * h) h = .1 * h;

            xnew = x + h;
            goto goto1;
        }
        else
        {
            if (errmax > ERRCON) hnext = SAFETY * h * (pow(errmax, PGROW));
            else hnext = 5. * h;

            hdid = h;
            x = x + h;

            for (i = 0; i < n; i++)
                y[i] = ytemp[i];
        }

        free(yerr);
        free(ytemp);
    }


    void rkck(double* y, double* dydx, int n, double x, double h, double* yout, double* yerr, void(*derivs)(int, double, double*, double*))
    {
        int i;
        const double A2 = .2, A3 = .3, A4 = .6, A5 = 1., A6 = .875, B21 = .2, B31 = 3. / 40.,
            B32 = 9. / 40., B41 = .3, B42 = -.9, B43 = 1.2, B51 = -11. / 54., B52 = 2.5,
            B53 = -70. / 27., B54 = 35. / 27., B61 = 1631. / 55296., B62 = 175. / 512.,
            B63 = 575. / 13824., B64 = 44275. / 110592., B65 = 253. / 4096., C1 = 37. / 378.,
            C3 = 250. / 621., C4 = 125. / 594., C6 = 512. / 1771., DC1 = C1 - 2825. / 27648.,
            DC3 = C3 - 18575. / 48384., DC4 = C4 - 13525. / 55296., DC5 = -277. / 14336.,
            DC6 = C6 - .25;

        double* ak2, * ak3, * ak4, * ak5, * ak6, * ytemp;

        ak2 = (double*)malloc(n * sizeof(double));
        ak3 = (double*)malloc(n * sizeof(double));
        ak4 = (double*)malloc(n * sizeof(double));
        ak5 = (double*)malloc(n * sizeof(double));
        ak6 = (double*)malloc(n * sizeof(double));
        ytemp = (double*)malloc(n * sizeof(double));

        for (i = 0; i < n; i++)
            ytemp[i] = y[i] + B21 * h * dydx[i];

        derivs(n, x + A2 * h, ytemp, ak2);

        for (i = 0; i < n; i++)
            ytemp[i] = y[i] + h * (B31 * dydx[i] + B32 * ak2[i]);

        derivs(n, x + A3 * h, ytemp, ak3);

        for (i = 0; i < n; i++)
            ytemp[i] = y[i] + h * (B41 * dydx[i] + B42 * ak2[i] + B43 * ak3[i]);

        derivs(n, x + A4 * h, ytemp, ak4);

        for (i = 0; i < n; i++)
            ytemp[i] = y[i] + h * (B51 * dydx[i] + B52 * ak2[i] + B53 * ak3[i] + B54 * ak4[i]);

        derivs(n, x + A5 * h, ytemp, ak5);

        for (i = 0; i < n; i++)
            ytemp[i] = y[i] + h * (B61 * dydx[i] + B62 * ak2[i] + B63 * ak3[i] + B64 * ak4[i] + B65 * ak5[i]);

        derivs(n, x + A6 * h, ytemp, ak6);

        for (i = 0; i < n; i++)
            yout[i] = y[i] + h * (C1 * dydx[i] + C3 * ak3[i] + C4 * ak4[i] + C6 * ak6[i]);

        for (i = 0; i < n; i++)
            yerr[i] = h * (DC1 * dydx[i] + DC3 * ak3[i] + DC4 * ak4[i] + DC5 * ak5[i] + DC6 * ak6[i]);

        free(ak2);
        free(ak3);
        free(ak4);
        free(ak5);
        free(ak6);
        free(ytemp);
    }



    void fderivs(int n, double x, double* y, double* dydx)
    {
        double z;
        double OMEGA_DE, f_OM, f_OR, f_Q, f;

        z = exp(-x) - 1.0;

        OMEGA_DE = 1.0 - OMEGA_R - OMEGA_M;
        f_OM = OMEGA_M * exp(-3.0 * x) * pow(inv_hubblez(z), 2);
        f_OR = OMEGA_R * exp(-4.0 * x) * pow(inv_hubblez(z), 2);
        f_Q = OMEGA_DE * f_DE(z) * pow(inv_hubblez(z), 2);
        f = y[1];

        dydx[0] = f;
        dydx[1] = -pow(f, 2) - f / 2.0 * (1.0 - f_OR - 3.0 * (w_0 - z * w_a) * f_Q) + 3.0 / 2.0 * f_OM;
    }

#ifdef __cplusplus
}
#endif