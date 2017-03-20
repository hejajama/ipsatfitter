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
    
    // Add datafiles, if 2nd parameter=true, then this is only charmdata
    data.LoadData("./data/hera_combined_sigmar.txt", TOTAL);
    //data.LoadData("data/hera_combined_sigmar_cc.txt", CHARM); // charm data

    
    MnUserParameters parameters;
    parameters.Add("B_G", 4.0);
    parameters.Add("heavy_mass", 1.27);
    parameters.Add("light_mass", 0.05); // Having very small mass is numerically difficult
    parameters.Add("mu_0", sqrt(1.51) );  // From Amir&Raju
    // mu_0 can't be free currently, as it is hardcoded in DGLAP evolution
    // When mu_0 is changed, one should also change alpha_s s.t. one still keeps
    // e.g. as(M_z) = 0.1184
    
    // Start using parameters by Amir&Raju
    parameters.Add("lambda_g", 0.072, 0.01);
    parameters.Add("A_g", 2.41, 0.4);
    
    parameters.SetLowerLimit("A_g", 0);
    
    //parameters.SetPrecision(0.001);
    
    cout << "=== Initial parameters ===" << endl;
    
    cout << parameters << endl;
    
    cout << "=== Starting fit ===" << endl;
    
    DISFitter fitter(parameters);
    fitter.AddDataset(data);
    
    
    
    // MnMinimize: use MIGRAD, if it fails, fall back to SIMPLEX
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
