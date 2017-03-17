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


using namespace std;
using namespace ROOT::Minuit2;
void ErrHandler(const char * reason,
                const char * file,
                int line,
                int gsl_errno);

int main()
{
    gsl_set_error_handler(&ErrHandler);
    
        
    //DISFitter fitter;
    
    
    
    Data data;
    data.SetMaxQsqr(50);
    data.LoadData("data/hera_combined_sigmar.txt", false);
    data.LoadData("data/hera_combined_sigmar_cc.txt", true);

    
    MnUserParameters parameters;
    parameters.Add("B_G", 4.0);
    parameters.Add("heavy_mass", 1.27);
    parameters.Add("light_mass", 0.05); // Having very small mass is numerically difficult
    parameters.Add("mu_0", sqrt(1.51) );  // From Amir&Raju
    parameters.Add("lambda_g", 0.06, 0.5);  // Some reasonable parameters to start
    parameters.Add("A_g", 2.3, 0.5);
    
    parameters.SetLowerLimit("A_g", 0);
    parameters.SetLowerLimit("lambda_g", 0); // Defined to be positive
    //parameters.SetLowerLimit("mu_0", 1);
    
    cout << "=== Initial parameters ===" << endl;
    
    cout << parameters << endl;
    
    cout << "=== Starting fit ===" << endl;
    
    DISFitter fitter(parameters);
    fitter.AddDataset(data);
    
    
    
    
    
    MnMigrad migrad(fitter, parameters);
    
    
    // minimize
    FunctionMinimum min = migrad();
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
