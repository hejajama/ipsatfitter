#include <iostream>
#include <string>
#include <gsl/gsl_errno.h>
#include <vector>

#include "wave_function.hpp"
#include "virtual_photon.hpp"
#include "ipsat.hpp"
#include "dis.hpp"
#include "data.hpp"
#include "interpolation2d.hpp"
#include <gsl/gsl_interp2d.h>

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
    data.LoadData("data/hera_combined_sigmar.txt", false);
    Data charmdata;
    charmdata.LoadData("data/hera_combined_sigmar_cc.txt", true);
    
    MnUserParameters parameters;
    parameters.Add("B_G", 4.0);
    parameters.Add("heavy_mass", 1.1, 0.1);
    parameters.Add("light_mass", 0.01);
    parameters.Add("mu_0",sqrt(1.51) );
    parameters.Add("lqcd",0.156);
    
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
    
    if (gsl_errno == 15) return; // Ugly hack, comes from the edges of the z integral in virtual_photon.cpp
    
    errors++;
    std::cerr << file << ":"<< line <<": Error " << errors << ": " <<reason
    << " (code " << gsl_errno << ")." << std::endl;
}
