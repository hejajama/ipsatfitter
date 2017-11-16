
/*
 * IPsat class for fitting
 * Heikki Mäntysaari <mantysaari@bnl.gov>, 2017
 */

#include <cmath>
#include <iostream>
#include <sstream>
#include <vector>
#include "ipsat.hpp"
#include "wave_function.hpp"
#include "woodsaxon.hpp"
#include "interpolation.hpp"

//#include "dglap_sartre/Dglap.h"

#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_expint.h>
#include <Minuit2/MnUserParameterState.h>

using namespace std;

const double MINB = 1e-7;
const double MAXB = 200;
const double BINTACCURACY = 0.001;
const int BINTEGRATIONDEPTH = 540;

using namespace std;

#ifdef INCLUDE_SARTRE_DGLAP
DglapEvolution sartre_dglap;
#endif

// IPsat 2012 - to test (extrat xg from this)

extern "C" {
    double dipole_amplitude_(double* xBj, double* r, double* b, int* param);
};
int IPSAT12_PAR = 2;    // 1: m_c=1.27 GeV,   2: m_c=1.4GeV

const bool USE_AMIR_FIT = false;

// LO DGLAP solver
// Modified by H.M. such that
// gluon = alphas(x,Q^2) * gluon!
extern "C"
{
    void lo_evol_(double *x, double* Q2, double *gluon, int* coupling,
                  double *mc, double *mb, double* mu0 , double *asmur,
                  double* Ag, double* lambdag, double *As, double* lambdas
                   );
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
    if (USE_AMIR_FIT)
        return dipole_amplitude_(&x, &r, &b, &IPSAT12_PAR)/2.0;
    
    double mu_0 = parameters.values->at( parameters.parameter->Index("mu_0"));
    double C = parameters.values->at( parameters.parameter->Index("C"));
    
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
    double B = parameters.values->at( parameters.parameter->Index("B_G"));
    if (config == -1 and saturation and !USE_AMIR_FIT and A==1 )
    {
        // Assume Gaussian profile exp(-b^2/(2B)) in the IPsat, can calculate
        // b integral analytically, as
        // \int d^2 b (1-Exp[ -a * exp(-b^2/(2B) ] )
        // = 2\pi B (gamma_E - CoshIntegral[a] + Log[a] + SinhIntegral[a]
        // GSL can evaluate CoshIntegral and SinhIntgral, so this is most likely the
        // most effective way to calculate b integral
        //
        // Note that in this notation a = pi^2 r^2 / (2NC) alphas * xg / (2 pi B_p)
        
        double mu_0 = parameters.values->at( parameters.parameter->Index("mu_0"));
        double C = parameters.values->at( parameters.parameter->Index("C"));
        double musqr = mu_0*mu_0 + C / SQR(r);
        double a = SQR(M_PI*r)/(2.0*NC) * Alphas(musqr, parameters) * xg(x, musqr, parameters) / (2.0 * M_PI * B);
       if (a==0) // Basically so small r that xg =0 as we are outside the dglap evolution grid
			   return 0; 
        
        gsl_sf_result sinres;
        gsl_sf_result cosres;
        int sinint = gsl_sf_Shi_e(a, &sinres);
        int cosint = gsl_sf_Chi_e(a, &cosres);
        
        if (!sinint and !cosint and cosres.val < 1e4 and sinres.val < 1e4)
        {
            // Require sinint and cosint to be so small that we can reliably calculate their difference
            // Otherwise fall back to numerical integration
  
            analres = 2.0*M_PI*B * ( M_EULER - cosres.val + log(a) + sinres.val);
            return analres;
        }
        
    }
    else if (config == -1 and !saturation and !USE_AMIR_FIT and A==1 )
    {
        // b integral analytically, now this is trivial as \int d^2 T_b = 1
        // so actually the result is 2\pi B N(r, b=0)
        return 2.0 * M_PI * B * DipoleAmplitude(r, 0, x, parameters, -1);
    }
    
    
    // b integral numerically
    
    gsl_function fun; fun.function=inthelperf_bint;
    inthelper_bint par;
    par.r=r; par.x=x; par.config=config;
    par.parameters = parameters;
    par.ipsat = this;
    fun.params=&par;
    
    double acc = BINTACCURACY;
    
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
    if (x < minx or x>maxx or musqr < minQ2 or musqr > maxQ2 )
    {
        cerr << "xg evaluated outside the validity range, x=" << x << ", mu^2=" << musqr << endl;
        return 0;
    }
    
    double lambdag =parameters.values->at( parameters.parameter->Index("lambda_g"));
    double Ag =parameters.values->at( parameters.parameter->Index("A_g"));
    double lambdas =parameters.values->at( parameters.parameter->Index("lambda_s"));
    double As =parameters.values->at( parameters.parameter->Index("A_s"));
    double mu0 =parameters.values->at( parameters.parameter->Index("mu_0"));
    double mc = parameters.values->at( parameters.parameter->Index("charm_mass"));
    double mb = parameters.values->at( parameters.parameter->Index("bottom_mass"));
    
    double gluon = 0;
    
    int coupling = 0;
    if (GetSinglet())
        coupling = 1;
    
    if (dglapsolver == PIA)
    {
        lo_evol_(&x, &musqr, &gluon, &coupling, &mc, &mb,
             &mu0, &parameters.alphas_mur,
             &Ag, &lambdag, &As, &lambdas );
        return gluon;
    }
#ifdef INCLUDE_SARTRE_DGLAP
    else if (dglapsolver == SARTRE)
    {
        gluon = x*sartre_dglap.G(x, musqr);
    }
#endif
    else if (dglapsolver == CPPPIA)
    {
        return parameters.cppdglap->alphasxG(x, musqr, mu0, coupling, Ag,
                   lambdag,  As,  lambdas);
    }
    else
    {
        cerr << "Unknown DGLAP solver!" << endl;
        exit(1);
    }
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
    double B_G = parameters.values->at( parameters.parameter->Index("B_G"));
    if (A>1)
    {
        return density_interpolator->Evaluate(b);
        
    }
    if (config != -1)
    {
        cerr << "Event-by-event fluctuations for the proton are not supproted at this point!" << endl;
        return 0;
    }
    
    return 1.0 / (2.0 * M_PI * B_G) * exp(-b*b / (2.0 * B_G));
}

/*
 * Strong coupling
 */
double IPsat::Alphas(double musqr, FitParameters parameters) const
{
    if (dglapsolver == PIA)
    {
        // Currently, xg returns actually alphas*xg, because I modified the
        // Fortran code...
        
        return 1.0;
    }
    else if (dglapsolver == CPPPIA)
    {
        return 1.0;
    }
    else
    {
        cerr << "Unknown DGLAP solver" << endl;
        return -1;
    }
    
}

/*
 * Init SARTRE DGLAP solver
 */
void IPsat::InitializeDGLAP(FitParameters par) const
{
    
}

/*
 * Constructor, set general settings
 */
IPsat::IPsat()
{
    minx=1e-10;
    maxx=0.1;
    minQ2=0;
    maxQ2=1e99;
    saturation = true;
    enable_singlet = false;
    maxalphas = 0.5;
    
    dglapsolver = CPPPIA;
    density_interpolator=NULL;
    A=1;
}

/*
 * Change to nuclear mode
 */
void IPsat::InitNucleus(int A_)
{

    A=A_;
    
    // Init for nucleus
    WoodsSaxon nuke(197);
    vector<double> bvals;
    vector<double> tavals;
    for (double b=0; b<100; b+=0.1)
    {
        bvals.push_back(b);
        tavals.push_back(A*nuke.T_A(b));
    }
    density_interpolator = new Interpolator(bvals, tavals);
    density_interpolator->SetOverflow(0);
    density_interpolator->SetUnderflow(0);
    density_interpolator->SetFreeze(true);

}

IPsat::~IPsat()
{
    if (A>1 and density_interpolator != NULL)
        delete density_interpolator;
}


std::ostream& operator<<(std::ostream& os, IPsat& ipsat)
{
    
    if (ipsat.GetSaturation())
        os << "IPsat";
    else
        os << "IPnonsat";

    os << ", singlet contribution: ";
    if (ipsat.GetSinglet())
        os << "enabled";
    else
        os << "disabled";
    return os;
    
}
