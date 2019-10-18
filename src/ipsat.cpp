
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

using namespace ROOT::Minuit2;



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
    cout << "Why I'm here?" << endl;
    
    
    
    
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
    double total_gammap;    // Used when we have "lumpy nucleus"
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
    double qs02 = parameters.values->at( parameters.parameter->Index("qs02"));
    double lambda = parameters.values->at( parameters.parameter->Index("lambda"));
    double gamma = parameters.values->at( parameters.parameter->Index("gamma"));
    double sigma02 = parameters.values->at( parameters.parameter->Index("sigma02"));
    double lnx0 = parameters.values->at( parameters.parameter->Index("lnx0"));
    double x0 = std::exp(lnx0);
    double qs2 = qs02 * std::pow(x0/x, lambda);
    
    return 2.0*M_PI*sigma02*(1.0 - std::exp(-std::pow(r*r*qs2, gamma)/4.0));
    /*
    double B = parameters.values->at( parameters.parameter->Index("B_G"));
    int A =parameters.values->at( parameters.parameter->Index("A"));
    if (config == -1 and saturation and A==1 )
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
        
        
        gsl_sf_result sinres;
        gsl_sf_result cosres;
        int sinint = gsl_sf_Shi_e(a, &sinres);
        int cosint = gsl_sf_Chi_e(a, &cosres);
        
        if (!sinint and !cosint and cosres.val < 1e4 and sinres.val < 1e4)
        {
            // Require sinint and cosint to be so small that we can reliably calculate their difference
            // Otherwise fall back to numerical integration
  
            return 2.0*M_PI*B * ( M_EULER - cosres.val + log(a) + sinres.val);
        }
        
    }
    else if (config == -1 and !saturation and A==1 )
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
     */
}

/*
 * b integral for nucleus, assuming "lumpy nucleus" from KT hep-ph/0304189 e.q. 47
 */
double inthelperf_bint_lumpyA(double b, void* p)
{
    inthelper_bint* par = (inthelper_bint*) p;
    int A = par->parameters.values->at( par->parameters.parameter->Index("A"));
    double TA = par->ipsat->GetDensityInterpolator()->Evaluate(b)/A ;
    
    // For ligh nulcei do not take large-A limit
    if (A < 100)
    {
        //cout << 1.0 - TA*par->total_gammap/2.0  << " " << 1.0 - std::pow(1.0 - TA*par->total_gammap/2.0 , A) << endl;
        double power  = std::pow(1.0 - TA*par->total_gammap/2.0 , A);
        
        // Sometimes we get that power > 1, which happens if total gammap xs is too large
        // NOTE: This happens always at crazy large r, because the total dipole-proton xs
        // goes like log(r). Of course, the photon wave function kills these crazy contributions.
        // But here we have to return 0, otherwise we may end up returning a large (negative) number which is not
        // damped by the virtual photon wave function
        if (power > 1) return 0;
        return b * (1.0 - power );
    }
    else
    {
        
        return b * (1.0 - std::exp( -A*TA*par->total_gammap/2.0));
    }
}
double IPsat::DipoleAmplitude_bint_lumpyA(double r, double x, FitParameters parameters, int config) const
{
    
    // First calculate total gamma-p
    // Need to create a local copy of MnUserParameters
    MnUserParameters tmp_parameters;
    // Constants
    tmp_parameters.Add("B_G", parameters.values->at( parameters.parameter->Index("B_G")) );
    tmp_parameters.Add("light_mass", parameters.values->at( parameters.parameter->Index("light_mass")));
    tmp_parameters.Add("charm_mass", parameters.values->at( parameters.parameter->Index("charm_mass")));
    tmp_parameters.Add("bottom_mass", parameters.values->at( parameters.parameter->Index("bottom_mass")));
    tmp_parameters.Add("C", parameters.values->at( parameters.parameter->Index("C")));
    tmp_parameters.Add("mu_0", parameters.values->at( parameters.parameter->Index("mu_0")));
    tmp_parameters.Add("lambda_g", parameters.values->at( parameters.parameter->Index("lambda_g")));
    tmp_parameters.Add("A_g", parameters.values->at( parameters.parameter->Index("A_g")));
    tmp_parameters.Add("lambda_s", 0);
    tmp_parameters.Add("A_s", 0);
    tmp_parameters.Add("A", 1);
   
    
    FitParameters tmp_proton_parameters;
    tmp_proton_parameters.parameter = &tmp_parameters;
    tmp_proton_parameters.alphas_mur = parameters.alphas_mur;
    tmp_proton_parameters.cppdglap = parameters.cppdglap;
    tmp_proton_parameters.alpha_strong = parameters.alpha_strong;
    
    vector<double> tmp;
    for (int i=0; i < tmp_parameters.Params().size(); i++)
    {
        tmp.push_back(tmp_parameters.Params()[i]);
    }
    tmp_proton_parameters.values = &tmp;
    double total_gammp = 2.0*DipoleAmplitude_bint(r, x, tmp_proton_parameters);
   
    
        /// gamma_p done
    
    gsl_function fun; fun.function=inthelperf_bint_lumpyA;
    inthelper_bint par;
    par.r=r; par.x=x; par.config=config;
    par.parameters = parameters;
    par.total_gammap = total_gammp;
    par.ipsat = this;
    fun.params=&par;
    
    double result,abserr;
    gsl_integration_workspace* ws = gsl_integration_workspace_alloc(BINTEGRATIONDEPTH);
    int status = gsl_integration_qag(&fun, MINB, MAXB, 0, BINTACCURACY,
                                     BINTEGRATIONDEPTH, GSL_INTEG_GAUSS51, ws, &result, &abserr);
    if (status)
        cerr << "bintegral failed in IPsat::DipoleAmplitude_bint_lumpyA with r=" << r <<", result " << result << " relerror " << abserr/result << endl;
    gsl_integration_workspace_free(ws);
    
    return 2.0*M_PI*result;
    
    
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
    
    if (A>1 )
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
    maxx=0.6;
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
    WoodsSaxon nuke(A);
    vector<double> bvals;
    vector<double> tavals;
    for (double b=0; b<100; b+=0.1)
    {
        bvals.push_back(b);
        tavals.push_back(A * nuke.T_A(b));
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
