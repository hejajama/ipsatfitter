#ifndef _IPSAT_H
#define _IPSAT_H

#include <vector>
#include <string>
#include <Minuit2/MnUserParameters.h>

#include "wave_function.hpp"

using namespace std;
using namespace ROOT::Minuit2;

/*
 * IPsat class that supports external DGLAP evolved xg(x,mu^2)
 * Designed to be easily extended to support multiple parameters
 */

struct FitParameters
{
    const vector<double> *values;
    const MnUserParameters *parameter;
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
    
    // Evaluate amplitude
    //config allows to specify configuration index in case of eventy-by-event
    // fluctuations (-1 = no fluctuations or average)
    // r and b are in GeV^(-1), x is Bjorken-x
    double DipoleAmplitude(double r, double b, double x, FitParameters parameter, int config=-1) const;
    
    // Dipole amplitude integrated over d^2 b
    double DipoleAmplitude_bint(double r, double x, FitParameters parameters, int config=-1) const;
    
    // Proton density profile in transverse plane
    // b is length of 2d vector, [b]=GeV^(-1)
    // config allows to specify configuration index in case of eventy-by-event
    // fluctuations
    double Tp(double b, FitParameters parameters, int config=-1) const;
    
    
    double Alphas(double musqr, FitParameters parameters) const;
    
    // Evaluate DGLAP evovled xg(x, mu^2), will call an external DGLAP solver
    double xg(double x, double musqr, FitParameters parameters) const;
    
    
    
private:

    
    // General constraints, in which the model can be used
    double minx;
    double maxx;
    double minQ2;
    double maxQ2;
    // Flag to turn on/off saturation
    // If false, dipole amplitude ~ r^2
    bool saturation;
    
    double maxalphas;
    
};


#endif
