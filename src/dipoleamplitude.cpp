/*
 * Solve DGLAP equation with given parametrization
 * and evaluate dipole amplitude
 *
 * Uses exactly the same DGLAP solver (LO_evolution_routine.f and
 * alphaS.f) as is used in the fit code
 *
 * H. Mäntysaari and P. Zurita, 2017
 */

#include "dipoleamplitude.hpp"
#include <string>
#include <cmath>
#include <iostream>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>

using namespace std;

// Example
int main()
{
    // Parameters are: DipoleAmplitude(C, mu0 [GeV], lambda_g, A_g, m_c [GeV]
    DipoleAmplitude amplitude(4, 1.35179, 0.09708, 2.29831, 1.37497);
    
    cout << "mu^2=5, as*xg (x=0.001) = " << amplitude.Alphas_xg(0.001, 5) << endl;
    
    
    
    return 0;
}


///////////////////////////////
// DipoleAmplitude class, wrapper to all DGLAP fortran etc


DipoleAmplitude::DipoleAmplitude(double C_, double mu0_, double lambda_g_, double A_g_, double mc_, double mb_, double mt_)
{
    C=C_; mu0=mu0_; lambda_g = lambda_g_; A_g = A_g_; mc = mc_; mb = mb_; mt = mt_;
    saturation = true;
    B_p=4.0;    // All fits are done with fixed B_p=4
    Nc=3;
    
    alphas_mur = InitAlphas();
    init_(); // Init Mellin momenta
}

double DipoleAmplitude::Alphas_xg(double x, double musqr)
{
    double as_xg=0;
    int coupling = 0;
    double As=0;
    double lambdas=0; //singlet
    
    lo_evol_(&x, &musqr, &as_xg, &coupling, &mc, &mb,
             &mu0, &alphas_mur,
             &A_g, &lambda_g, &As, &lambdas );
}

double DipoleAmplitude::Alphas(double Q)
{
    return alphas_(&Q);
}

double DipoleAmplitude::xg(double x, double musqr)
{
    return Alphas_xg(x, musqr) / Alphas(std::sqrt(musqr));
}

double DipoleAmplitude::N(double r, double xbj, double b)
{
    double musqr = mu0 + C / r*r;
    double exponent = M_PI*M_PI / (2.0 * Nc) * r*r * Alphas_xg(xg, musqr) * Tp(b);
    
    if (!saturation)     // IPnonsat
        return exponent;
    else
        return 1.0 - std::exp(-exponent);

}

double DipoleAmplitude::Tp(double b)
{
    return 1.0/(2.0*M_PI*B_p) * std::exp(-b*b / (2.0*B_p));
}

/////
// Alpha_s initialization

double alphas_helper(double asmur, void* p);
double DipoleAmplitude::InitAlphas()
{
    // Init Alpha_s
    // Find alphas at the initial scale
    gsl_function f; f.function = &alphas_helper;
    f.params = this;
    
    const gsl_root_fsolver_type *T = gsl_root_fsolver_brent;
    gsl_root_fsolver *s = gsl_root_fsolver_alloc (T);
    double min = 0.1; double max = 0.5;
    gsl_root_fsolver_set (s, &f, min, max);
    
    
    int iter=0; int maxiter = 100;
    double asmur = 0.4; int status;
    do
    {
        iter++;
        gsl_root_fsolver_iterate (s);
        asmur = gsl_root_fsolver_root (s);
        min = gsl_root_fsolver_x_lower (s);
        max = gsl_root_fsolver_x_upper (s);
        status = gsl_root_test_interval (min, max,   0, 0.001);
        
    }
    while (status == GSL_CONTINUE && iter < maxiter);
    
    gsl_root_fsolver_free(s);
    
    if (iter >= maxiter)
    {
        cerr << "Initializing alphas failed! " << endl;
        exit(1);
    }
    
    return asmur;
    
}
// Helper function to solve asmur s.t. alphas(M_z)=0.11884
double alphas_helper(double asmur, void* p)
{
    //cout << "Trying " << asmur << endl;
    // Init alphas
    DipoleAmplitude *dipole = (DipoleAmplitude*)p;
    int iord=0;
    double fr2 = 1.0;
    double mur = dipole->GetMu0();
    double mc = dipole->GetMc();
    double mb = dipole->GetMb();
    double mt = dipole->GetMt();
    if (mur > mc)
    {
        cerr << "mu_r > m_c at alphas_helperf(), not possible! " << endl;
        exit(1);
    }
    
    initalphas_(&iord, &fr2, &mur, &asmur, &mc, &mb, &mt);
    double mz=91.1876;
    return alphas_(&mz) - 0.1183; // From HERA II 1506.06042
}
