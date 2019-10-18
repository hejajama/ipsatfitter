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

enum INITIAL_PARAMETERS
{
    IPSAT_MAXQ2_50,
	IPSAT_MAXQ2_500,
    IPNONSAT_MAXQ2_50,
	IPNONSAT_MAXQ2_500
};


int main(int argc, char* argv[])
{
    gsl_set_error_handler(&ErrHandler);
    
    Data data;
    data.SetMinQsqr(1.49);
    data.SetMaxQsqr(50.1);
    
    // Add datafiles, if 2nd parameter=CHARM, then this is only charmdata
    data.LoadData("./data/hera_combined_sigmar.txt", TOTAL);
    data.LoadData("./data/hera_combined_sigmar_eminusp.txt", TOTAL);
    //data.LoadData("data/hera_combined_sigmar_cc.txt", CHARM, 1.0); // charm data

    //data.LoadData("data/light_quark_f2/hera_I_combined_eplus_lightq", UDS);
    
    
    //data.LoadData("data/hera_II_combined_sigmar_b.txt", BOTTOM, 1.0);
    
    INITIAL_PARAMETERS start_params = IPSAT_MAXQ2_50;
    
    MnUserParameters parameters;
    // Constants
    
    parameters.Add("light_mass", 0.14);
    parameters.Add("charm_mass", 1.4);
    parameters.Add("qs02", 1);
    parameters.Add("lnx0", std::log(0.0004), 0.1);
    parameters.Add("lambda", 0.2, 0.1);     // Qs^2 = qs2*(0.01/x)^lambda
    parameters.Add("gamma", 0.6);            // Anomalous dimension in GBW
    parameters.Add("sigma02", 30, 5);
    
    parameters.Add("A", 1);
    
    parameters.SetLowerLimit("lambda", 0);
    parameters.SetLowerLimit("gamma", 0);
    parameters.SetLowerLimit("lnx0", std::log(1e-6));
    parameters.SetUpperLimit("lnx0", std::log(0.5));

    DISFitter fitter(parameters);
	fitter.MAXR = 99;
    /*fitter.MAXR = StrToReal(argv[1])*5.068;
	cout << "MAXR set to " << fitter.MAXR/5.068 << " fm" << endl;
     */
    fitter.AddDataset(data);
    
    if (start_params == IPSAT_MAXQ2_50 or start_params == IPSAT_MAXQ2_500)
        fitter.SetSaturation(true);
    else
        fitter.SetSaturation(false);
    fitter.SetSinglet(false);
    
    fitter.SetDGLAPSolver(CPPPIA);
    
    FitParameters p;
    p.parameter = &parameters;
    vector<double> parvec = parameters.Params();
    p.values = &parvec;
    
    
    // Initialize alpha_s(M_Z=91.1876)=0.1183
    /*
    AlphaStrong *alphas = new AlphaStrong(0, 1.0, 91.1876, 0.1183, parvec[parameters.Index("charm_mass")], 4.75, 175);
    // DGLAP_Solver will take care of deleting alphas when it is deleted
    EvolutionLO_gluon *cppdglap = new EvolutionLO_gluon(alphas);
    cppdglap->generateLookupTable(parvec[parameters.Index("mu_0")], 0, parvec[parameters.Index("A_g")], parvec[parameters.Index("lambda_g")], 0, 0);
    cppdglap->useLookupTable(true);
    
    p.cppdglap = cppdglap;
    p.alpha_strong = alphas;
    */
    /*double q2vals[4] = {2,15, 50, 500};
    double totvals_f2[4];
    double totvals_fl[4];
    double x= 0.0001;
    for (int q2ind=0; q2ind < 4; q2ind++)
    {
        fitter.MAXR=99;
        totvals_f2[q2ind] = fitter.F2(q2vals[q2ind], x, p );
        //cout << "f2c at " << q2vals[q2ind] << " : " <<totvals_f2[q2ind] << endl;
        totvals_fl[q2ind]=fitter.FL(q2vals[q2ind], x, p );
    }
    for (double maxr=0.02; maxr < 5; maxr += 0.02)
    {
        fitter.MAXR = maxr * 5.068;
        cout << maxr << " ";
        for (int q2ind=0; q2ind < 4; q2ind++)
        {
            double q2 = q2vals[q2ind];
            double f2 = fitter.F2(q2, x, p );
            double fl =fitter.FL(q2, x, p );
            cout << f2/totvals_f2[q2ind] << " " << fl/totvals_fl[q2ind] << " ";
        }
        cout << endl;
    }
    exit(1);*/
    
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
