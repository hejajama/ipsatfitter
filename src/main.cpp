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
    data.SetMinQsqr(1.49);
    data.SetMaxQsqr(50.1);
    
    // Add datafiles, if 2nd parameter=CHARM, then this is only charmdata
     data.LoadData("./data/hera_combined_sigmar.txt", TOTAL);
     data.LoadData("./data/hera_combined_sigmar_eminusp.txt", TOTAL);
     data.LoadData("data/hera_combined_sigmar_cc.txt", CHARM, 1.0); // charm data
     //data.LoadData("data/lhec_ncepp.txt", TOTAL);
    
    MnUserParameters parameters;
    // Constants
    parameters.Add("B_G", 4.0);
    
    //parameters.Add("light_mass", 0.0005);
    //parameters.Add("charm_mass", 1.4);
    
    //parameters.Add("light_mass", 0.1388639702255); // Having very small mass is numerically difficult
    //parameters.Add("charm_mass", 1.342035015621,  0.1 ); // 1.27 // Ipsat 1.361410284911 // Nonsat 1.350324669808,
    
    parameters.Add("light_mass", 0.03);
    parameters.Add("charm_mass", 1.35165);
    //parameters.Add("charm_mass", 1.4);
    parameters.Add("bottom_mass", 4.75);  // 4.75
    //parameters.Add("C", 4.939286653112, 1.0);
    parameters.Add("C", 2.146034445992);
    //parameters.Add("C", 2.41367);
    // Start using some reasonable parameters
    
    
    
    parameters.Add("mu_0", 1.1 );
    
    parameters.Add("lambda_g", 0.09665075464199, 0.02);
    parameters.Add("A_g", 2.103826220003, 0.4);
    
    // maxQ2 500
    //parameters.Add("lambda_g", 0.09661, 0.02);
    //parameters.Add("A_g", 2.0667, 0.4);
    
    //parameters.Add("lambda_g", -0.009631194037871, 0.02);
    //parameters.Add("A_g", 3.058791613883, 0.4);
    
    parameters.Add("lambda_s", 0);
    parameters.Add("A_s", 0);
    


	parameters.Add("lambda_s", 0);
	parameters.Add("A_s", 0);
    parameters.Add("A", 1);
    
    // Set limits
    
    parameters.SetLowerLimit("A_g", 0);
    parameters.SetLowerLimit("A_s", 0);
    parameters.SetLowerLimit("mu_0", 0.5); // In priciple can go to anything >0 (right?)
    parameters.SetUpperLimit("mu_0", sqrt(2.7)); // Upper limit is min Q^2
    parameters.SetUpperLimit("charm_mass", 3);
    DISFitter fitter(parameters);
    fitter.AddDataset(data);
    
    fitter.SetSaturation(true);
    fitter.SetSinglet(false);
    
    fitter.SetDGLAPSolver(CPPPIA);
    
    FitParameters p;
    p.parameter = &parameters;
    vector<double> parvec = parameters.Params();
    p.values = &parvec;
    
    
    // Initialize alpha_s(M_Z=91.1876)=0.1183
    AlphaStrong *alphas = new AlphaStrong(0, 1.0, 91.1876, 0.1183, parvec[parameters.Index("charm_mass")], 4.75, 175);
    // DGLAP_Solver will take care of deleting alphas when it is deleted
    EvolutionLO *cppdglap = new EvolutionLO(alphas);
    
    p.cppdglap = cppdglap;
    p.alpha_strong = alphas;
    
    

    /*
    //double xvals[7] = {0.338E-5, 0.501E-5, 0.696E-5, 0.103E-4, 0.149E-4, 0.238E-4, 0.476E-4}; //Q^2=2
    double xvals[7]= { 0.422E-4, 0.627E-4, 0.870E-4, 0.128E-3, 0.186E-3, 0.298E-3, 0.595E-3};
    //double xvals[7] = {0.169E-3, 0.251E-3, 0.348E-3, 0.513E-3, 0.744E-3, 0.119E-2, 0.238E-2}; // Q^2=100
    for (int i=0; i<7; i++)
    {
        double x = xvals[i];
        cout << x << " " << fitter.F2(2,x,p) << " " << fitter.FL(2,x,p) << endl;
    }
    exit(1);
    /*
     cout << "#x  F_2(Q^2=2)    F_2(Q^2=3.5)    F_2(Q^2=6.5)   F_2(Q^2=12)  F_2(Q^2=25)  F_2(Q^2=100)  F_L(Q^2=2)    F_L(Q^2=3.5)    F_L(Q^2=6.5)   F_L(Q^2=12)  F_L(Q^2=25)  F_L(Q^2=100) " << endl;
     for (double x=1e-8; x<0.02; x*=1.5)
    //for (int i=0; i<7; i++)
    {
        //double x = xvals[i];
        cout << x << " " << fitter.F2(2, x, p) << " " << fitter.F2(3.5, x, p) << " " << fitter.F2(6.5, x, p) << " " << fitter.F2(12,x,p) << " " << fitter.F2(25,x,p) << " " << fitter.F2(100,x,p) << " " <<  fitter.FL(2, x, p) << " " << fitter.FL(3.5, x, p) << " " << fitter.FL(6.5, x, p) << " " << fitter.FL(12,x,p) << " " << fitter.FL(25,x,p) << " " << fitter.FL(100,x,p) << endl;
    }
    exit(1);
    

        /*
    cout << "#x   F_L(Q^2=2)   F_L(Q^2=5)     F_L(Q^2=6.5)    F_L(12)     F_L(25)     F_L(100)    F_L(400)    F_(2000) " << endl;
    for (double x=1e-8; x<1e-2; x*=1.2)
    {
        cout << x << " " << fitter.FL(2, x, p) << " " << fitter.FL(3.5, x, p) << " " << fitter.FL(6.5, x, p) << " " << fitter.FL(12,x,p) << " " << fitter.FL(25, x, p) << " " << fitter.FL(100, x, p) << " " << fitter.FL(400,x,p) << " " <<    fitter.FL(500, x, p) << endl;
    }
    exit(1);
    */
    
    
    cout << "# Q^2   F_2 F_L" << endl;
    for (double q2=1; q2<1000000; q2*=1.4)
    //for (double q2 =10000*5; q2<10000000; q2*=5 )
    {
        cout << q2 << " " << fitter.F2(q2, 0.1, p) << " " << fitter.FL(q2, 0.1, p) << endl;
    }
    
    
    exit(1);
    
    
    //cout << fitter.FL(2, 1e-2, p) << " " << fitter.FL(15, 1e-2, p) << " " <<fitter.FL(50, 1e-2, p) << " " << fitter.FL(150, 1e-2, p)<< " " << fitter.FL(500, 1e-2, p) << endl;
    //exit(1);

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
