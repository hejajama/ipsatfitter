
/*
 * IPsat class for fitting
 * Heikki Mäntysaari <mantysaari@bnl.gov>, 2017
 */

#include <cmath>
#include <iostream>
#include "ipsat.hpp"
#include "wave_function.hpp"

#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_expint.h>
#include <Minuit2/MnUserParameterState.h>

using namespace std;

const double MINB = 1e-7;
const double MAXB = 2e2;
const double BINTACCURACY = 0.001;
const int BINTEGRATIONDEPTH = 54;

using namespace std;


// IPsat 2012 - to test (extrat xg from this)
/*
extern "C" {
    double dipole_amplitude_(double* xBj, double* r, double* b, int* param);
};
int IPSAT12_PAR = 1;    // 1: m_c=1.27 GeV,   2: m_c=1.4GeV
*/

// LO DGLAP solver
// Modified by H.M. such that
// gluon = alphas(x,Q^2) * gluon!
extern "C"
{
    void lo_evol_(double *x, double* Q2, double *gluon, int* coupling, double* Ag, double* lambdag, double* mu0 , double *asmur );
    //void init_(); // Init Mellin momenta
};


/*
 * Evaluate dipole amplitude
 *
 * config can be used to specify a given eventy-by-event configuration, default is -1
 * which refers to non-fluctuating case
 */
double IPsat::DipoleAmplitude(double r, double b, double x, FitParameters parameters,  int config) const
{
    // Use Amir's dipole
    //return dipole_amplitude_(&x, &r, &b, &IPSAT12_PAR)/2.0;
    
    double mu_0 = parameters.values->at( parameters.parameter->Index("mu_0"));
    double C = 4;
    
    double musqr = mu_0*mu_0 + C / SQR(r);
    
    double exponent = SQR(M_PI*r)/(2.0*NC) * Alphas(musqr, parameters) * xg(x, musqr, parameters) * Tp(b, parameters, config);
    if (saturation)
        return 1.0 - exp(-exponent);
    else
        return exponent;
    
    
}


/*
 * Dipole amplitude integrated over d^2 b
 * \int d^2b IPsat::DipoleAmplitude(double r, double b, double x, int config)
 */
struct inthelper_bint
{
    double r,x;
    int config;
    const IPsat* ipsat;
    FitParameters parameters;
};
double inthelperf_bint(double b, void* p)
{
    inthelper_bint* par = (inthelper_bint*) p;
    return b*par->ipsat->DipoleAmplitude(par->r, b, par->x, par->parameters, par->config);
}

double IPsat::DipoleAmplitude_bint(double r, double x, FitParameters parameters, int  config) const
{
    double analres=0;
    
    if (config == -1 and saturation )
    {
        // Assume Gaussian profile exp(-b^2/(2B)) in the IPsat, can calculate
        // b integral analytically, as
        // \int d^2 b (1-Exp[ -a * exp(-b^2/(2B) ] )
        // = 2\pi B (gamma_E - CoshIntegral[a] + Log[a] + SinhIntegral[a]
        // GSL can evaluate CoshIntegral and SinhIntgral, so this is most likely the
        // most effective way to calculate b integral
        //
        // Note that in this notation a = pi^2 r^2 / (2NC) alphas * xg / (2 pi B_p)
        double B = parameters.values->at( parameters.parameter->Index("B_G"));
        double mu_0 = parameters.values->at( parameters.parameter->Index("mu_0"));
        double C = 4;
        double musqr = mu_0*mu_0 + C / SQR(r);
        double a = SQR(M_PI*r)/(2.0*NC) * Alphas(musqr, parameters) * xg(x, musqr, parameters) / (2.0 * M_PI * B);
        
        
        
        gsl_sf_result sinres;
        gsl_sf_result cosres;
        int sinint = gsl_sf_Shi_e(a, &sinres);
        
        if (!sinint)
        {
            // No overflows
            int cosint = gsl_sf_Chi_e(a, &cosres);
            
            if (!cosint)
            {
                // No overflows, use analytical result, otherwise we fall back to numerics
                // in the region where the contribution anyway is small
                analres = 2.0*M_PI*B * ( M_EULER - cosres.val + log(a) + sinres.val);
                return analres;
            }
            else
                return 0;
        }
        else
        {
            return 0;
            // Overflow, very large r most likely. Fall back to numerical integral.
            //cerr << "Error with argument " << a << " r " << r << endl;
        }

        
        
        
        
        
        
    }
    
    
    // b integral numerically
    
    gsl_function fun; fun.function=inthelperf_bint;
    inthelper_bint par;
    par.r=r; par.x=x; par.config=config;
    par.parameters = parameters;
    par.ipsat = this;
    fun.params=&par;
    
    double acc = BINTACCURACY;
    if (r < 1e-5)
        acc*=10; // At very small dipoles this is numerically dificult, but also
        // contribution is very small, so let's not demand too much from GSL
    
    double result,abserr;
    gsl_integration_workspace* ws = gsl_integration_workspace_alloc(BINTEGRATIONDEPTH);
    int status = gsl_integration_qag(&fun, MINB, MAXB, 0, acc,
                                     BINTEGRATIONDEPTH, GSL_INTEG_GAUSS51, ws, &result, &abserr);
    if (status)
        cerr << "bintegral failed in IPsat::DipoleAmplitude_bit with r=" << r <<", result " << result << " relerror " << abserr/result << endl;
    gsl_integration_workspace_free(ws);
    
    return 2.0*M_PI*result; //2pi from angular integral
}

double IPsat::xg(double x, double musqr, FitParameters parameters) const
{
    if (x < minx or x>maxx or musqr < minQ2 or musqr > maxQ2)
    {
        cerr << "xg evaluated outside the validity range, x=" << x << ", mu^2=" << musqr << endl;
        return 0;
    }
    // For testing: initial condition
    //return A_g * std::pow(x, -lambda_g) * pow((1.0 - x), 5.6);
    
    double lambdag =parameters.values->at( parameters.parameter->Index("lambda_g"));
    double Ag =parameters.values->at( parameters.parameter->Index("A_g"));
    double mu0 =parameters.values->at( parameters.parameter->Index("mu_0"));
    
    double gluon=0;
    int coupling = 0;
    lo_evol_(&x, &musqr, &gluon, &coupling, &Ag, &lambdag, &mu0, &parameters.alphas_mur);
    
    return  gluon;
    
    
}

/*
 * Transverse profile
 * Here, the standard parametrization is T = 1/(2\pi B_G) exp(-b^2/(2B_G))
 *
 * In order to support event-by-event fluctuations, there is parameter config 
 * by default, config = -1, and in that case we have the round non-fluctuating profile
 */
double IPsat::Tp(double b, FitParameters parameters, int config) const
{
    if (config != -1)
    {
        cerr << "Event-by-event fluctuations for the proton are not supproted at this point!" << endl;
        return 0;
    }
    double B_G = parameters.values->at( parameters.parameter->Index("B_G"));
    return 1.0 / (2.0 * M_PI * B_G) * exp(-b*b / (2.0 * B_G));
}

/*
 * Strong coupling
 */
double IPsat::Alphas(double musqr, FitParameters parameters) const
{
    // Currently, xg returns actually alphas*xg, because I modified the
    // Fortran code...
    
    return 1.0;
    
}

/*
 * Constructor, set general settings
 */
IPsat::IPsat()
{
    minx=1e-9;
    maxx=0.02;
    minQ2=0;
    maxQ2=1e99;
    saturation = true;
    maxalphas = 0.5;
    
    
    
    
}
