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

int main(int argc, char* argv[])
{
    gsl_set_error_handler(&ErrHandler);
    
    Data data;
    data.SetMinQsqr(0.74);
    data.SetMaxQsqr(650.1);
    
    // Add datafiles, if 2nd parameter=CHARM, then this is only charmdata
        data.LoadData("./data/hera_combined_sigmar.txt", TOTAL);
    data.LoadData("./data/hera_combined_sigmar_eminusp.txt", TOTAL);
    //data.LoadData("data/hera_combined_sigmar_cc.txt", CHARM, 1.0); // charm data

    
    MnUserParameters parameters;
    // Constants
    parameters.Add("B_G", 4.0);
    
    //parameters.Add("charm_mass", 1.3231);
    
    parameters.Add("light_mass", 0.1388639702255); // Having very small mass is numerically difficult
    parameters.Add("charm_mass", 1.342035015621,  0.1 ); // 1.27 // Ipsat 1.361410284911 // Nonsat 1.350324669808,
    //parameters.Add("light_mass", 0.03);
    //parameters.Add("charm_mass", 1.354062489611);
    parameters.Add("bottom_mass", 4.75);  // 4.75
    parameters.Add("C", 4.939286653112, 1.0);
    //parameters.Add("C", 2.321526423259);
    //parameters.Add("C", 2.41367);
    // Start using some reasonable parameters
    
    
    
    parameters.Add("mu_0", 1.1 );
    
    //parameters.Add("lambda_g", 0.09106887412584, 0.02);
    //parameters.Add("A_g", 2.155203998342, 0.4);
    
    // maxQ2 500
    //parameters.Add("lambda_g", 0.09661, 0.02);
    //parameters.Add("A_g", 2.0667, 0.4);
    
    parameters.Add("lambda_g", -0.009631194037871, 0.02);
    parameters.Add("A_g", 3.058791613883, 0.4);
    
    parameters.Add("lambda_s", 0);
    parameters.Add("A_s", 0);
    


	parameters.Add("lambda_s", 0);
	parameters.Add("A_s", 0);

    
    // Set limits
    
    parameters.SetLowerLimit("A_g", 0);
    parameters.SetLowerLimit("A_s", 0);
    parameters.SetLowerLimit("mu_0", 0.5); // In priciple can go to anything >0 (right?)
    parameters.SetUpperLimit("mu_0", sqrt(2.7)); // Upper limit is min Q^2
    parameters.SetUpperLimit("charm_mass", 3);
    DISFitter fitter(parameters);
    fitter.AddDataset(data);
    
    fitter.SetSaturation(false);
    fitter.SetSinglet(false);
    
    fitter.SetDGLAPSolver(CPPPIA);
    
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
    cout << "# Q^2   F_2(x=1e-2)    F_L(x=1e-2)   F_2(x=5e-3)   F_L(x=5e-3)    F_2(x=1e-3)   F_L(x=1e-3)    F_L(x=1e-4)    F_2(x=1e-4)   F_L(x=1e-5)    F_2(x=1e-5) " << endl;
    for (double q2=1; q2<10000; q2*=1.2)
    {
        cout << q2 << " " << fitter.F2(q2, 1e-2, p) << " " << fitter.FL(q2, 1e-2, p) << " " << fitter.F2(q2, 5e-3, p) << " " << fitter.FL(q2, 5e-3, p)  << " " << fitter.F2(q2, 1e-3, p) << " " << fitter.FL(q2, 1e-3, p) << " " << fitter.F2(q2, 1e-4, p) << " " << fitter.FL(q2, 1e-4, p) << " " << fitter.F2(q2, 1e-5, p) << " " << fitter.FL(q2, 1e-5, p)  << endl;
    }
    exit(1);
    */
    
     //parameters.SetPrecision(0.001);
    
    cout << "=== Initial parameters ===" << endl;
    
    cout << parameters << endl;
    cout << "DGLAP solver: ";
    if (fitter.GetDGLAPSolver() == PIA) cout << "PIA" << endl;
    else if (fitter.GetDGLAPSolver() == SARTRE) cout << "SARTRE" << endl;
    else if (fitter.GetDGLAPSolver() == CPPPIA) cout << "C++ version of PIA's" << endl;
    
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
