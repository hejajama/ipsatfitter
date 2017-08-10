#include <iostream>
#include <string>
#include <gsl/gsl_errno.h>
#include <vector>

#include "wave_function.hpp"
#include "virtual_photon.hpp"
#include "ipsat.hpp"
#include "dis.hpp"
#include "data.hpp"


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
    data.SetMinQsqr(0.75);
    data.SetMaxQsqr(650);
    
    // Add datafiles, if 2nd parameter=CHARM, then this is only charmdata
    data.LoadData("./data/hera_combined_sigmar.txt", TOTAL);
    data.LoadData("data/hera_combined_sigmar_cc.txt", CHARM); // charm data

    
    MnUserParameters parameters;
    // Constants
    parameters.Add("B_G", 4.0);
    parameters.Add("light_mass", 0.05); // Having very small mass is numerically difficult
    parameters.Add("charm_mass", 1.270); // 1.27
    parameters.Add("bottom_mass", 4.75);
    parameters.Add("C", 4.0);
    
    
    // Start using some reasonable parameters
    
    // IPsat, mc fitted to 1.381
    /* parameters.Add("mu_0", 1.339, 0.2 );
    parameters.Add("lambda_g", 0.0955, 0.02);
    parameters.Add("A_g", 2.309, 0.4);
    parameters.Add("lambda_s", 0);
    parameters.Add("A_s", 0);
    */

	// Amir
	parameters.Add("mu_0", 1.188, 0.01);
    parameters.Add("A_g", 2.373, 0.1);
    parameters.Add("lambda_g", 0.043, 0.01);

    //parameters.Add("mu_0", sqrt(1.428), 0.01);
    //parameters.Add("A_g", 2.373, 0.1);
    //parameters.Add("lambda_g", 0.052, 0.01);
	parameters.Add("lambda_s", 0);
	parameters.Add("A_s", 0);
    // IPsat with singlet
    /*
    parameters.Add("mu_0", 1.538, 0.2 );
    parameters.Add("lambda_g", 0.1129, 0.02);
    parameters.Add("A_g", 2.233, 0.4);
    parameters.Add("lambda_s", 0);
    parameters.Add("A_s", 0);
    */
    
    // IPsat, mc=1.27
    /*
    parameters.Add("mu_0", 1.352, 0.1);
    parameters.Add("lambda_g", 0.0921, 0.05);
    parameters.Add("A_g", 2.309, 0.1);
    parameters.Add("lambda_s", 0);
    parameters.Add("A_s", 0);
    */
    
    // IPsat without singlet, free charm mass, fitted to all data
    /*
    parameters.Add("mu_0", 1.338937588098, 0.2 );
    parameters.Add("lambda_g", 0.09548245405113, 0.02);
    parameters.Add("A_g", 2.308901027904, 0.4);
    parameters.Add("lambda_s", 0);
    parameters.Add("A_s", 0);
    */
    
    // IPsat fitted to charm data
    /*
    parameters.Add("mu_0", 0.8381738043555, 0.2 );
    parameters.Add("lambda_g", 0.06520164038353, 0.02);
    parameters.Add("A_g", 2.002630013683, 0.4);
    parameters.Add("lambda_s", 0);
    parameters.Add("A_s", 0);
     */
    
    
    // IPnonsat w/o singlet
    /*
    parameters.Add("mu_0", 0.706, 0.1);
    parameters.Add("lambda_g", -0.199, 0.05);
    parameters.Add("A_g", 4.105, 0.1);
    parameters.Add("lambda_s", 0);
    parameters.Add("A_s", 0);
    */
    
    // IPnonsat with singlet
    /*
    parameters.Add("mu_0", 0.7223, 0.1);
    parameters.Add("lambda_g", -0.1759, 0.05);
    parameters.Add("A_g", 3.781, 0.1);
    parameters.Add("lambda_s", 0);
    parameters.Add("A_s", 0);
    */
    
    // Set limits
    
    parameters.SetLowerLimit("A_g", 0);
    parameters.SetLowerLimit("A_s", 0);
    parameters.SetLowerLimit("mu_0", 0.4); // In priciple can go to anything >0 (right?)
    //parameters.SetUpperLimit("mu_0", 1.43); // For some reason alphas code does not work with larger mu_0
    
    DISFitter fitter(parameters);
    fitter.AddDataset(data);
    
    fitter.SetSaturation(true);
    fitter.SetSinglet(false);
    
    fitter.SetDGLAPSolver(SARTRE);
    
    FitParameters p;
    p.parameter = &parameters;
    vector<double> parvec = parameters.Params();
    p.values = &parvec;
   
	/*
    for (double q2=1; q2<10000; q2*=1.2)
    {
        cout << q2 << " " << fitter.F2(q2, 5e-3, p) << " " << fitter.FL(q2, 5e-3, p) << endl;
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
