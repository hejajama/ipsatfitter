#ifndef _IPSAT_H
#define _IPSAT_H

#include <vector>
#include <string>
#include <iostream>
#include <Minuit2/MnUserParameters.h>

#ifdef INCLUDE_SARTRE_DGLAP
    #include "dglap_sartre/DglapEvolution.h"
#endif

#include "dglap_cpp/AlphaStrong.h"
#include "dglap_cpp/EvolutionLO_coupling.h"
#include "dglap_cpp/EvolutionLO_nocoupling.h"
#include "wave_function.hpp"
#include "interpolation.hpp"

using namespace std;
using namespace ROOT::Minuit2;

enum DGLAP_Solver
{
    SARTRE,     // Sartre version
    PIA,         // Pia's LO_evolution_routine.f
    CPPPIA      // Pia's code translated into C++
};


/*
 * IPsat class that supports external DGLAP evolved xg(x,mu^2)
 * Designed to be easily extended to support multiple parameters
 */

// Struct to pass list of current parameters and Minuit parameter info object
struct FitParameters
{
    const vector<double> *values;
    const MnUserParameters *parameter;
    double alphas_mur;  // alphas at initial scale, solved by DISFitter
    
    // Needed for C++ DGLAP code
    AlphaStrong* alpha_strong;
    EvolutionLO_gluon *cppdglap;
};

// Note: In general this is designed to support event-by-event fluctuatoins
// Thus there is an additional parameter "config" floating around
// Default is config=-1 and no fluctuations
// When we implement fluctuations, config=-1 will refer to numerically
// averaged profile (notice that F_2 is linear in dipole amplitude, so
// average can be calculated at that level, but the situation is different
// for diffractive observables!



class IPsat
{
public:
    IPsat();
    void InitNucleus(int A_);
    ~IPsat();
    
    // Evaluate amplitude
    //config allows to specify configuration index in case of eventy-by-event
    // fluctuations (-1 = no fluctuations or average)
    // r and b are in GeV^(-1), x is Bjorken-x
    double DipoleAmplitude(double r, double b, double x, FitParameters parameter, int config=-1) const;
    
    // Dipole amplitude integrated over d^2 b
    double DipoleAmplitude_bint(double r, double x, FitParameters parameters, int config=-1) const;
    
    // Integrate over d^2b assuming "lumpy nucleus" ala KT
    double DipoleAmplitude_bint_lumpyA(double r, double x, FitParameters parameters, int config=-1) const;
    
    // Proton density profile in transverse plane
    // b is length of 2d vector, [b]=GeV^(-1)
    // config allows to specify configuration index in case of eventy-by-event
    // fluctuations
    double Tp(double b, FitParameters parameters, int config=-1) const;
    
    
    double Alphas(double musqr, FitParameters parameters) const;
    
    // Evaluate DGLAP evovled xg(x, mu^2), will call an external DGLAP solver
    double xg(double x, double musqr, FitParameters parameters) const;
    
    bool GetSaturation() const { return saturation; }
    void SetSaturation(bool s) { saturation = s; }
    
    bool GetSinglet() const { return enable_singlet; }
    void SetSinglet(bool s) { enable_singlet = s; }
    
    void SetDGLAPSolver(DGLAP_Solver s) {  dglapsolver = s; }
    
    void InitializeDGLAP(FitParameters par) const;
    Interpolator* GetDensityInterpolator() const { return density_interpolator; }
    
    
private:

    
    // General constraints, in which the model can be used
    double minx;
    double maxx;
    double minQ2;
    double maxQ2;
    // Flag to turn on/off saturation
    // If false, dipole amplitude ~ r^2
    bool saturation;
    bool enable_singlet;
    
    double maxalphas;
    
    DGLAP_Solver dglapsolver;
    
    Interpolator *density_interpolator;  // Used for nuclei
    int A;  // If >1, use WoodsSaxon nucleus
    
    
    
};

std::ostream& operator<<(std::ostream& os, IPsat& ipsat);


#endif
