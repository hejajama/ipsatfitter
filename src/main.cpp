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
    data.SetMaxQsqr(50);
    
    // Add datafiles, if 2nd parameter=CHARM, then this is only charmdata
    data.LoadData("./data/hera_combined_sigmar.txt", TOTAL);
    //data.LoadData("data/hera_combined_sigmar_cc.txt", CHARM); // charm data

    
    MnUserParameters parameters;
    // Constants
    parameters.Add("B_G", 4.0);
    parameters.Add("light_mass", 0.05); // Having very small mass is numerically difficult
    parameters.Add("charm_mass", 1.27);
    parameters.Add("bottom_mass", 4.75);
    parameters.Add("C", 4.0);
    
    
    // Start using some reasonable parameters
    
    // Singlet fit
    /*
    parameters.Add("mu_0", 1.44628, 0.2 );
    parameters.Add("lambda_g", 0.1111, 0.02);
    parameters.Add("A_g", 2.16795, 0.4);
    parameters.Add("lambda_s", -0.26047, 0.1);
    parameters.Add("A_s", 0.0961, 0.2);
    */
    
    // IPsat
    // Gives chi^2/N~1.46 in range 0.75-50
    parameters.Add("mu_0", 1.31777, 0.1);
    parameters.Add("lambda_g", 0.0878, 0.05);
    parameters.Add("A_g", 2.3369, 0.1);
    parameters.Add("lambda_s", 0);
    parameters.Add("A_s", 0);
    
    
    // Set limits
    parameters.SetLowerLimit("A_g", 0);
    parameters.SetLowerLimit("A_s", 0);
    parameters.SetLowerLimit("mu_0", 0.4); // In priciple can go to anything >0 (right?)
    //parameters.SetUpperLimit("mu_0", 1.43); // For some reason alphas code does not work with larger mu_0
    
    DISFitter fitter(parameters);
    fitter.AddDataset(data);
    
    fitter.SetSaturation(true);
    fitter.SetSinglet(false);
    
    FitParameters p;
    p.parameter = &parameters;
    vector<double> parvec = parameters.Params();
    p.values = &parvec;
    
    /*
    for (double q2=2; q2<10000; q2*=1.2)
    {
        cout << q2 << " " << fitter.F2(q2, 0.01, p) << endl;
    }*/
    
    //parameters.SetPrecision(0.001);
    
    cout << "=== Initial parameters ===" << endl;
    
    cout << parameters << endl;
    
    
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
