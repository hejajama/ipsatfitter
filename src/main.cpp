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
    data.LoadData("data/hera_combined_sigmar.txt", false);
    data.SetMaxQsqr(50);
    Data charmdata;
    charmdata.LoadData("data/hera_combined_sigmar_cc.txt", true);
    
    MnUserParameters parameters;
    parameters.Add("B_G", 4.0);
    parameters.Add("heavy_mass", 1.27);
    parameters.Add("light_mass", 0.05);
    parameters.Add("mu_0",sqrt(1.51), 0.5 );
    parameters.Add("lambda_g", 0.05, 0.1);
    parameters.Add("A_g", 2.37, 0.1);
    
    //parameters.Add("lqcd",0.156);
    
    //parameters.SetLowerLimit("light_mass", 0);
    //parameters.SetLowerLimit("heavy_mass", 0);
    parameters.SetLowerLimit("A_g", 0);
    parameters.SetLowerLimit("lambda_g", 0);
    parameters.SetUpperLimit("mu_0", 100);
    
    DISFitter fitter(parameters);
    fitter.AddDataset(data);
    fitter.AddDataset(charmdata);
    
    
    
    
    
    MnMigrad migrad(fitter, parameters);
    
    
    // minimize
    FunctionMinimum min = migrad();
    // output
    std::cout<<"minimum: "<<min<<std::endl;
    
    
    
    // Test code to calculate F_2, just to compare with HERA data and test that
    // everything works
    /*
     double q2=10;
     FitParameters param;
     param.parameter = &parameters;
     vector<double> vals; vals.push_back(4); vals.push_back(1.4); vals.push_back(0.05); vals.push_back(sqrt(1.51)); vals.push_back(0.156);
     param.values = &vals;
     
     VirtualPhoton light;
     light.SetQuark(LIGHT, 0.01);
     double x=5e-3;
     //ProtonPhotonCrossSection(const double Qsqr, const double xbj, const Polarization pol,const VirtualPhoton* wf , FitParameters fitparams) const
     double light_l = fitter.ProtonPhotonCrossSection(q2, x, LONGITUDINAL, &light, param );
     double light_t = fitter.ProtonPhotonCrossSection(q2, x, TRANSVERSE, &light, param );
     cout << "Light " << 10.0/(4.0*SQR(M_PI)*ALPHA_e)*(light_t + light_l) << endl;
     
     VirtualPhoton charm; charm.SetQuark(C, 1.4);
     x = x*(1.0+4.0*1.4*1.4/10.0);
     double c_l = fitter.ProtonPhotonCrossSection(10, x, LONGITUDINAL, &charm, param );
     double c_t = fitter.ProtonPhotonCrossSection(10, x, TRANSVERSE, &charm, param );
     cout << "Charm " << 10.0/(4.0*SQR(M_PI)*ALPHA_e)*(c_t + c_l) << endl;
     cout << "Total F2 " <<10.0/(4.0*SQR(M_PI)*ALPHA_e)*(light_t + light_l + c_t + c_l) << " exp 0.633 pm 0.0105" << endl;
     
     exit(1);
     */
    
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
    
    if (gsl_errno == 15 or gsl_errno == 16) return; // Ugly hack, comes from the edges of the z integral in virtual_photon.cpp
    // Overflows come from IPsat::bint when it is done analytically
    
    errors++;
    std::cerr << file << ":"<< line <<": Error " << errors << ": " <<reason
    << " (code " << gsl_errno << ")." << std::endl;
}
