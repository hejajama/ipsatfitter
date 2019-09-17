
/*
 * Solve DGLAP equation with given parametrization
 * and evaluate dipole amplitude
 *
 * Uses exactly the same DGLAP solver (LO_evolution_routine.f and
 * alphaS.f) as is used in the fit code
 *
 * H. Mäntysaari, 2018
 *
 * Read in HERA reduced cross section data, and subtract at all points
 * heavy quark contribution, and output "light quark sigma_r"
 */

#include "../src/data.hpp"
#include "../src/dglap_cpp/AlphaStrong.h"
#include "../src/dglap_cpp/EvolutionLO_nocoupling.h"
#include "../src/virtual_photon.hpp"
#include "../src/wave_function.hpp"
#include "../src/dis.hpp"
#include <string>
#include <cmath>
#include <iostream>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>
#include <sstream>
#include <iomanip>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_errno.h>
#include <Minuit2/FCNBase.h>
#include <Minuit2/MnUserParameterState.h>
#include <Minuit2/MnPrint.h>

using namespace std;
void ErrHandler(const char * reason,
                const char * file,
                int line,
                int gsl_errno);


const double AS_MZ = 0.1183;
using namespace ROOT::Minuit2;

double StrToReal(std::string str)
{
    std::stringstream buff(str);
    double tmp;
    buff >> tmp;
    return tmp;
}

// Example
int main(int argc, char* argv[])
{
    gsl_set_error_handler(&ErrHandler);
    double mu0 = std::sqrt(1.1);
    
    // We will use DISFitter class, so need to set up MINUIT
    MnUserParameters parameters;
    double mc =1.35277437092;
    double mb = 4.75;
    double ml=0.03;
    // Constants
    parameters.Add("B_G", 4.0);
    
    parameters.Add("light_mass",  ml);
    parameters.Add("charm_mass", mc);
    parameters.Add("bottom_mass", mb);
    parameters.Add("C", 2.289363553168);
    parameters.Add("mu_0", std::sqrt(1.1));
    parameters.Add("lambda_g", 0.08289088639946);
    parameters.Add("A_g", 2.195310911936);
    
    parameters.Add("lambda_s", 0);
    parameters.Add("A_s", 0);
    parameters.Add("A", 1);
    
    DISFitter fitter(parameters);
    vector<double> parvec = parameters.Params();
    
    cout << "# IPsat parameters" << endl;
    cout << parameters << endl;
    
    // Initialize alpha_s(M_Z=91.1876)=0.1183
    
    AlphaStrong *alphas = new AlphaStrong(0, 1.0, 91.1876, 0.1183, mc, mb, 175);
    // DGLAP_Solver will take care of deleting alphas when it is deleted
    EvolutionLO_gluon *cppdglap = new EvolutionLO_gluon(alphas);
    cppdglap->generateLookupTable(parvec[parameters.Index("mu_0")], 0, parvec[parameters.Index("A_g")], parvec[parameters.Index("lambda_g")], 0, 0);
    cppdglap->useLookupTable(false);

    
    
    FitParameters params;
    params.parameter = &parameters;
    params.values = &parvec;
    
    params.alpha_strong = alphas;
    params.cppdglap = cppdglap;
    
    // Wave functions
    VirtualPhoton wf_lightquark;
    VirtualPhoton wf_charm;
    VirtualPhoton wf_bottom;
    wf_lightquark.SetQuark(LIGHT, ml);
    wf_charm.SetQuark(C, mc);
    wf_bottom.SetQuark(B, mb);
    
    // Add datafiles, if 2nd parameter=CHARM, then this is only charmdata
    Data data;
    data.SetMaxQsqr(1000.1);
    data.SetMinQsqr(0.74);
    //data.LoadData("./data/hera_combined_sigmar.txt", TOTAL);
    data.LoadData("./data/hera_combined_sigmar_eminusp.txt", TOTAL);
    
    cout << "# Q2 [GeV^2] " << setw(10) << " x " << setw(10) << " y " << setw(10) << " sigma_r " << setw(10) << " abs.err" << endl;
    
    double chisqr=0; int points=0;
    for (int i=0; i<data.NumOfPoints(); i++)
    {
        double x = data.xbj(i);
        double y = data.y(i);
        double Q2 = data.Qsqr(i);
        double sigmar = data.ReducedCrossSection(i);
        double sigmar_err = data.ReducedCrossSectionError(i);
        DataType point_type =data.DataPointType(i);
        double sqrts = sqrt( Q2 / (x * y) );
        
        double charmx = x * (1.0 + 4.0*mc*mc / Q2);
        double bottomx = x*(1.0 + 4.0*mb*mb / Q2);
        
        double relative_uncert = sigmar_err/sigmar;
        
        // charm contribution
        double th_sigmar_c = fitter.ReducedCrossSection(Q2, charmx, sqrts, &wf_charm, params);
        //double th_sigmar =fitter.ReducedCrossSection(Q2, x, sqrts, &wf_lightquark, params);
        
        double th_sigmar_b=0;
        if (bottomx < 0.1)
            th_sigmar_b = fitter.ReducedCrossSection(Q2, bottomx, sqrts, &wf_bottom, params);
        
        double new_sigmar = sigmar - th_sigmar_c - th_sigmar_b;
        
        cout << setw(10) << Q2 << " " << setw(10)  << x << " " << setw(10) << y << " " << setw(10) << new_sigmar << " "  << setw(10)  << sigmar_err << endl;
        
        
    }
    
    //cout << "chi^2/N = " << chisqr/points<< " (N=" << points << ")" << endl;
        
    
}

int errors;
void ErrHandler(const char * reason,
                const char * file,
                int line,
                int gsl_errno)
{
    
    // 14 = failed to reach tolerance
    // 18 = roundoff error prevents tolerance from being achieved
    // 11 = maximum number of subdivisions reached
    // 15: underflows
    
    if (gsl_errno == 15 or gsl_errno == 16) return;
    // Ugly hack, comes from the edges of the z integral in virtual_photon.cpp
    // Overflows come from IPsat::bint when it is done analytically
    // Hope is that these errors are handled correctly everywhere
    
    errors++;
    std::cerr << file << ":"<< line <<": Error " << errors << ": " <<reason
    << " (code " << gsl_errno << ")." << std::endl;
}

    
    
    
