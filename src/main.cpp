#include <iostream>
#include <string>
#include <stdexcept>
#include <gsl/gsl_errno.h>
#include <vector>
#include <algorithm>
#include <cctype>

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

namespace
{
    struct DatasetSpec
    {
        std::string filename;
        DataType type;
        double weight = 1.0;
        bool f2csv = false;     // F2 values in csv format, see Data::LoadF2CSV
    };

    std::string ToUpper(std::string s)
    {
        std::transform(s.begin(), s.end(), s.begin(),
                       [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
        return s;
    }

    bool TryParseDataType(const std::string& s, DataType& out)
    {
        const std::string t = ToUpper(s);
        if (t == "TOTAL") { out = TOTAL; return true; }
        if (t == "CHARM") { out = CHARM; return true; }
        if (t == "UDS")   { out = UDS;   return true; }
        return false;
    }

    bool IsKnownOption(const std::string& s)
    {
        return s == "--data" || s == "--f2data" || s == "--heraerrors" ||
               s == "--help" || s == "-h" ||
               s == "-minQ2" || s == "--minQ2" ||
               s == "-maxQ2" || s == "--maxQ2";
    }

    void PrintHelp(const char* prog)
    {
        std::cout
            << "Usage:\n"
            << "  " << prog << " [options]\n\n"
            << "Options:\n"
            << "  -h, --help                      Show this help and exit\n"
            << "  -minQ2, --minQ2 <value>         Minimum Q^2 (default: 1.49)\n"
            << "  -maxQ2, --maxQ2 <value>         Maximum Q^2 (default: 50.1)\n"
            << "  --data <file> <type> [weight]   Add dataset (repeatable)\n"
            << "                                  type: TOTAL | CHARM | UDS\n"
            << "                                  weight defaults to 1.0\n"
            << "  --f2data <file> <type> [weight] Add F2 dataset in csv format with columns\n"
            << "                                  \"Q2\",\"x\",\"F2\" (repeatable). Theory computes F2,\n"
            << "                                  uncertainty is the relative stat+uncorr. HERA\n"
            << "                                  uncertainty at the same (Q2,x)\n"
            << "  --heraerrors <file>             HERA data used for the --f2data uncertainties\n"
            << "                                  (default: data/hera_II_combined_sigmar.txt)\n\n"
            << "If no --data is given, defaults are loaded:\n"
            << "  ./data/hera_combined_sigmar_eminusp.txt TOTAL 1\n"
            << "  data/hera_combined_sigmar_cc.txt CHARM 1\n\n"
            << "Examples:\n"
            << "  " << prog << " --data ./data/hera_combined_sigmar_eminusp.txt TOTAL\n"
            << "  " << prog << " --data data/hera_combined_sigmar_cc.txt CHARM 0.5\n"
            << "  " << prog << " -minQ2 2.0 -maxQ2 100 --data mydata.txt UDS\n";
    }
}

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

    double minQ2 = 1.49;
    double maxQ2 = 50.1;
    std::vector<DatasetSpec> datasets;
    std::string heraerrorfile = "data/hera_II_combined_sigmar.txt";

    for (int i = 1; i < argc; ++i)
    {
        std::string arg = argv[i];

        if (arg == "-h" || arg == "--help")
        {
            PrintHelp(argv[0]);
            return 0;
        }
        else if (arg == "-minQ2" || arg == "--minQ2" || arg == "-maxQ2" || arg == "--maxQ2")
        {
            if (i + 1 >= argc)
            {
                std::cerr << "Missing value for " << arg << std::endl;
                return 1;
            }

            double value = 0.0;
            try
            {
                value = std::stod(argv[++i]);
            }
            catch (const std::exception&)
            {
                std::cerr << "Invalid numeric value for " << arg << ": " << argv[i] << std::endl;
                return 1;
            }

            if (arg == "-minQ2" || arg == "--minQ2") minQ2 = value;
            else maxQ2 = value;
        }
        else if (arg == "--heraerrors")
        {
            if (i + 1 >= argc)
            {
                std::cerr << "Missing value for " << arg << std::endl;
                return 1;
            }
            heraerrorfile = argv[++i];
        }
        else if (arg == "--data" || arg == "--f2data")
        {
            if (i + 2 >= argc)
            {
                std::cerr << "Usage error: " << arg << " requires FILENAME TYPE [WEIGHT]" << std::endl;
                return 1;
            }

            DatasetSpec spec;
            spec.f2csv = (arg == "--f2data");
            spec.filename = argv[++i];
            spec.weight = 1.0;

            const std::string typeArg = argv[++i];
            if (!TryParseDataType(typeArg, spec.type))
            {
                std::cerr << "Invalid dataset TYPE: " << typeArg
                          << " (expected TOTAL, CHARM, or UDS)" << std::endl;
                return 1;
            }

            if (i + 1 < argc)
            {
                std::string next = argv[i + 1];
                if (!IsKnownOption(next))
                {
                    try
                    {
                        spec.weight = std::stod(next);
                        ++i;
                    }
                    catch (const std::exception&)
                    {
                        std::cerr << "Invalid WEIGHT for --data " << spec.filename
                                  << ": " << next << std::endl;
                        return 1;
                    }
                }
            }

            datasets.push_back(spec);
        }
        else
        {
            std::cerr << "Unknown argument: " << arg << std::endl;
            std::cerr << "Use --help to see available options." << std::endl;
            return 1;
        }
    }

    if (minQ2 > maxQ2)
    {
        std::cerr << "Invalid Q2 range: minQ2 (" << minQ2 << ") is greater than maxQ2 (" << maxQ2 << ")." << std::endl;
        return 1;
    }
    
    Data data;
    data.SetMinQsqr(minQ2);
    data.SetMaxQsqr(maxQ2);
    
    if (datasets.empty())
    {
        data.LoadData("data/hera_combined_sigmar.txt", TOTAL, 1.0);
        data.LoadData("data/hera_combined_sigmar_cc.txt", CHARM, 1.0);
    }
    else
    {
        for (const auto& ds : datasets)
        {
            int status = ds.f2csv ? data.LoadF2CSV(ds.filename, heraerrorfile, ds.type, ds.weight)
                                  : data.LoadData(ds.filename, ds.type, ds.weight);
            if (status != 0)
                return 1;
        }
    }

    INITIAL_PARAMETERS start_params = IPSAT_MAXQ2_50;
    
    MnUserParameters parameters;
    // Constants
    parameters.Add("B_G", 4.0);
    
    // IPsat
    if (start_params == IPSAT_MAXQ2_50)
    {
        parameters.Add("light_mass", 0.03);
        parameters.Add("charm_mass", 1.35277437092, 0.01);
        parameters.Add("bottom_mass", 4.75, 0.01);
        parameters.Add("C", 2.289363553168, 0.01);
        parameters.Add("mu_0", std::sqrt(1.1));
        parameters.Add("lambda_g", 0.08289088639946, 0.01);
        parameters.Add("A_g", 2.195310911936, 0.01);
        /*
        parameters.Add("light_mass", 0.03);
        parameters.Add("charm_mass", 1.3210, 0.01);
        parameters.Add("bottom_mass", 4.75, 0.01);
        parameters.Add("C", 1.8178, 0.01);
        parameters.Add("mu_0", std::sqrt(1.1));
        parameters.Add("lambda_g", 0.09575, 0.01);
        parameters.Add("A_g", 2.1394, 0.01);
         */
        
    }
	if (start_params == IPSAT_MAXQ2_500)
    {
        parameters.Add("light_mass", 0.03);
        parameters.Add("charm_mass", 1.3296, 0.01);
        parameters.Add("bottom_mass", 4.75, 0.01);
        parameters.Add("C", 2.6477, 0.01);
        parameters.Add("mu_0", std::sqrt(1.1));
        parameters.Add("lambda_g", 0.0779, 0.01);
        parameters.Add("A_g", 2.2097, 0.01);
    }

    else if (start_params == IPNONSAT_MAXQ2_50)
    {
        // IPnonsat
        parameters.Add("light_mass", 0.1515769997484);
        parameters.Add("charm_mass", 1.350367375905);
        parameters.Add("bottom_mass", 4.75);
        parameters.Add("C", 4.297444629517);
        parameters.Add("mu_0", std::sqrt(1.1));
        parameters.Add("lambda_g", -0.006657294973805);
        parameters.Add("A_g", 3.039134356321);
    }
    else if (start_params == IPNONSAT_MAXQ2_500)
    {
        // IPnonsat
        parameters.Add("light_mass", 0.1332);
        parameters.Add("charm_mass", 1.3187);
        parameters.Add("bottom_mass", 4.75);
        parameters.Add("C", 5.6510);
        parameters.Add("mu_0", std::sqrt(1.1));
        parameters.Add("lambda_g", -0.03460);
        parameters.Add("A_g", 3.2820);
	}
 
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
    
    AlphaStrong *alphas = new AlphaStrong(0, 1.0, 91.1876, 0.1183, parvec[parameters.Index("charm_mass")], 4.75, 175);
    // DGLAP_Solver will take care of deleting alphas when it is deleted
    EvolutionLO_gluon *cppdglap = new EvolutionLO_gluon(alphas);
    cppdglap->generateLookupTable(parvec[parameters.Index("mu_0")], 0, parvec[parameters.Index("A_g")], parvec[parameters.Index("lambda_g")], 0, 0);
    cppdglap->useLookupTable(true);
    
    p.cppdglap = cppdglap;
    p.alpha_strong = alphas;
    
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
    /*
    double x = 5.6234e-05;
    double q2=10;
    parvec[parameters.Index("charm_mass")] = 0.03;
    parvec[parameters.Index("bottom_mass")] = 0.03;
    double mf0_F2 =fitter.F2(q2, x, p );
    double mf0_FL =fitter.FL(q2, x, p );
    
    
    for (double mc=0.03; mc < 4; mc += 0.1)
    {
        parvec[parameters.Index("charm_mass")] = mc;
        parvec[parameters.Index("bottom_mass")] = mc;
        parvec[parameters.Index("light_mass")] = mc;
        double F2 =fitter.F2(q2, x, p );
        double FL =fitter.FL(q2, x, p );
        
        cout << mc << " " << F2/mf0_F2 << " " << FL/mf0_FL << endl;
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
