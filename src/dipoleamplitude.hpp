/*
 * Solve DGLAP equation with given parametrization 
 * and evaluate dipole amplitude
 *
 * Uses exactly the same DGLAP solver (LO_evolution_routine.f and
 * alphaS.f) as is used in the fit code
 *
 * H. Mäntysaari and P. Zurita, 2017
 */

// LO DGLAP solver and alphas
// gluon = alphas(x,Q^2) * gluon!
// For documentation of the varaibles see alphaS.f and LO_evolution_routine.f
// DipoleAmplitude class takes care of calling these, so the user doesn't have to care!
extern "C"
{
    void lo_evol_(double *x, double* Q2, double *gluon, int* coupling,
                  double *mc, double *mb, double* mu0 , double *asmur,
                  double* Ag, double* lambdag, double *As, double* lambdas
                  );
    void init_();
    double alphas_(double *mu);
    void initalphas_(int *iord, double *fr2, double *mur, double* asmur, double *mc, double *mb, double* mt);
};

class DipoleAmplitude
{
    public:
        DipoleAmplitude(double C_, double mu0_, double lambda_g_, double A_g_, double mc_, double mb_=4.75, double mt_=175 );
    
        /*
         * Compute alphas(mu^2) * xg(x, mu^2), musqr in GeV^2
         */
        double Alphas_xg(double x, double musqr);
    
        // Strong coupling, Q in GeV
        double Alphas(double Q);
    
        // Only xg, calculated internally as Alphas_xg / Alphas, so not as effective!
        double xg(double x, double musqr);
    
        /*
         * Dipole amplitude
         * [r] = GeV^-1, dipole size
         * [b]=Gev^-1, impact parameter
         */
        double N(double r, double xbj, double b);
    
        /*
         * Proton profile, normalized to unity
         * \int d^2 b T_p = 1
         *
         * [b] = GeV^-1
         */
        double Tp(double b);
    
        // Other methods
        double GetMu0() { return mu0; }
        double GetLambdaG() { return lambda_g; }
        double GetAg() { return A_g; }
        double GetMc() { return mc; }
        double GetMb() { return mb; }
        double GetMt() { return mt; }
        bool GetSaturation() { return saturation; }
        void SetSaturation(bool s) { saturation = s; }
    
    private:
        bool saturation;   // True for ipsat, false for ipnonsat
        double InitAlphas();
        double C;
        double mu0;
        double lambda_g;
        double A_g;
        double mc;
        double mb;
        double mt;
        double B_p; // Proton size, GeV^-2
        double alphas_mur;      // alpha_s at mu0
        int Nc;
};
