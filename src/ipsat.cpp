
/*
 * IPsat class for fitting
 * Heikki Mäntysaari <mantysaari@bnl.gov>, 2017
 */

#include <cmath>
#include <iostream>
#include "ipsat.hpp"
#include "wave_function.hpp"

#include <gsl/gsl_integration.h>
#include <Minuit2/MnUserParameterState.h>

using namespace std;

const double MINB = 1e-7;
const double MAXB = 1e2;
const double BINTACCURACY = 0.01;
const int BINTEGRATIONDEPTH = 20;

using namespace std;


// IPsat 2012 - to test (extrat xg from this)
extern "C" {
    double dipole_amplitude_(double* xBj, double* r, double* b, int* param);
};

int IPSAT12_PAR = 1;    // m_c=1.27 GeV


/*
 * Evaluate dipole amplitude
 *
 * config can be used to specify a given eventy-by-event configuration, default is -1
 * which refers to non-fluctuating case
 */
double IPsat::DipoleAmplitude(double r, double b, double x, FitParameters parameters,  int config) const
{
    double mu_0 = parameters.values->at( parameters.parameter->Index("mu_0"));
    double C = 4;
    
    return dipole_amplitude_(&x,&r,&b,&IPSAT12_PAR)/2.0;
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
    gsl_function fun; fun.function=inthelperf_bint;
    inthelper_bint par;
    par.r=r; par.x=x; par.config=config;
    par.parameters = parameters;
    par.ipsat = this;
    fun.params=&par;
    
    double result,abserr;
    gsl_integration_workspace* ws = gsl_integration_workspace_alloc(BINTEGRATIONDEPTH);
    int status = gsl_integration_qag(&fun, MINB, MAXB, 0, BINTACCURACY,
                                     BINTEGRATIONDEPTH, GSL_INTEG_GAUSS51, ws, &result, &abserr);
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
    
    // Extrat xg from fit
    // mu2 = mu_0^2 + C/r^2
    // r^2 = C/(mu^2 - mu_0^2)
    double mu_0 = parameters.values->at( parameters.parameter->Index("mu_0"));
    double r=sqrt( C/ ( musqr - mu_0*mu_0));
    double b=0;
    
    double n =1 - dipole_amplitude_(&x, &r, &b, &IPSAT12_PAR)/2.0;
    if (n<=0) return 0;
    double tp = 1.0/(2.0*M_PI*4.0)*std::exp(- b*b / (2.0*4.0));
    
    
    double c = -log(n);
    c /= SQR(M_PI*r)/(2.0*NC) * Alphas(musqr, parameters) * tp;
    
    return c;
    
    
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
    // Todo: flavor scheme?
    double NF=4;
    double b0 = 11.0 - 2.0/3.0 * NF;
    
    double lqcd = parameters.values->at( parameters.parameter->Index("lqcd"));
    
    if (musqr / (lqcd*lqcd) < 0)
        return maxalphas;
    
    double as = 4.0*M_PI / (b0 * log(musqr / (lqcd*lqcd)));
    if (as > maxalphas)
        return maxalphas;
    return as;
    
}

/*
 * Constructor, set general settings
 */
IPsat::IPsat()
{
    minx=1e-9;
    maxx=0.02;
    minQ2=1;
    maxQ2=1e9;
    saturation = true;
    maxalphas = 0.5;
    
    
    
    
}
