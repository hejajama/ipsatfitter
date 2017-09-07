#include <iostream>
#include <string>
#include <gsl/gsl_errno.h>
#include <vector>

#include "wave_function.hpp"
#include "virtual_photon.hpp"
#include "ipsat.hpp"
#include "dis.hpp"
#include "data.hpp"

#ifdef INCLUDE_SARTRE_DGLAP
    #include "dglap_sartre/AlphaStrong.h"
#endif

#include <Minuit2/MnUserParameterState.h>
#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/MnApplication.h>
#include <Minuit2/MnPrint.h>
#include <Minuit2/MnScan.h>
#include <Minuit2/MnMinimize.h>
#include <Minuit2/MnSimplex.h>

using namespace std;
using namespace ROOT::Minuit2;
void ErrHandler(const char * reason,
                const char * file,
                int line,
                int gsl_errno);

int main()
{
    gsl_set_error_handler(&ErrHandler);
    
    Data data;
    data.SetMinQsqr(1.99);
    data.SetMaxQsqr(501);
    
    // Add datafiles, if 2nd parameter=CHARM, then this is only charmdata
    data.LoadData("./data/hera_combined_sigmar.txt", TOTAL);
    data.LoadData("./data/hera_combined_sigmar_eminusp.txt", TOTAL);
    data.LoadData("data/hera_combined_sigmar_cc.txt", CHARM, 1.0); // charm data

    
    MnUserParameters parameters;
    // Constants
    parameters.Add("B_G", 4.0);
    parameters.Add("light_mass", 0.05); // Having very small mass is numerically difficult
    parameters.Add("charm_mass", 1.3408, 0.1); // 1.27
    parameters.Add("bottom_mass", 4.75);
    parameters.Add("C", 4.0);
    
    
    // Start using some reasonable parameters
    
    // IPsat, mc fitted to 1.381
     parameters.Add("mu_0", 1.28802, 0.2 );
    parameters.Add("lambda_g", 0.0976, 0.02);
    parameters.Add("A_g", 2.1955, 0.4);
    parameters.Add("lambda_s", 0);
    parameters.Add("A_s", 0);
    


	parameters.Add("lambda_s", 0);
	parameters.Add("A_s", 0);

    
    // Set limits
    
    parameters.SetLowerLimit("A_g", 0);
    parameters.SetLowerLimit("A_s", 0);
    parameters.SetLowerLimit("mu_0", 0.5); // In priciple can go to anything >0 (right?)
    parameters.SetUpperLimit("mu_0", sqrt(2.0)); // Upper limit is min Q^2
    DISFitter fitter(parameters);
    fitter.AddDataset(data);
    
    fitter.SetSaturation(true);
    fitter.SetSinglet(false);
    
    fitter.SetDGLAPSolver(PIA);
    
    FitParameters p;
    p.parameter = &parameters;
    vector<double> parvec = parameters.Params();
    p.values = &parvec;
    
    
    /*
     fitter(parvec);
    cout << "Q^2 [GeV^2]    alphas*xg(x=0.01, Q^2)    alphas*xg(x=0.001 Q^2)" << endl;
    for (double q2 = 1.5; q2 < 1e7; q2*=1.1)
    {
        cout << q2 << " " << fitter.GetDipole().xg(0.01, q2, p) << " " << fitter.GetDipole().xg(0.001, q2, p) << endl;
    }
    exit(1);
   */
	/*
    double xbj=5e-3;
    cout << "#x = " << xbj <<  ", all quarks included" << endl;
    cout << "# Q^2   F_2    F_L" << endl;
    for (double q2=1; q2<10000; q2*=1.2)
    {
        cout << q2 << " " << fitter.F2(q2, xbj, p) << " " << fitter.FL(q2, xbj, p) << endl;
    }
    exit(1);
    */
    
     //parameters.SetPrecision(0.001);
    
    cout << "=== Initial parameters ===" << endl;
    
    cout << parameters << endl;
    cout << "DGLAP solver: ";
    if (fitter.GetDGLAPSolver() == PIA) cout << "PIA" << endl;
    else if (fitter.GetDGLAPSolver() == SARTRE) cout << "SARTRE" << endl;
    
    
    cout << "=== Dipole amplitude ===" << endl;
    cout << endl;
    cout << fitter.GetDipole() << endl<<endl;
    
    cout << "=== Starting fit ===" << endl;
    
    
    
    
    // MnMinimize: use MIGRAD, if it fails, fall back to SIMPLEX
    //MnSimplex fit(fitter,parameters);
    MnMinimize fit(fitter, parameters);
    //MnMigrad fit(fitter, parameters);
    //MnScan fit(fitter, parameters);
    
    
    // minimize
    FunctionMinimum min = fit();
    // output
    std::cout<<"minimum: "<<min<<std::endl;
    
    
    return 0;
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
