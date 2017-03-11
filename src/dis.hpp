#ifndef _DIS_H
#define _DIS_H

/*
 * Calculate gamma-p cross section and reduced cross section
 * Takes all the parameters in the MINUIT format, and calls IPsat
 * class with required parameters
 *
 * Heikki Mäntysaari <mantysaari@bnl.gov>, 2017
 */

#include <vector>
#include <Minuit2/FCNBase.h>
#include <Minuit2/MnUserParameterState.h>
#include "wave_function.hpp"
#include "virtual_photon.hpp"
#include "ipsat.hpp"
#include "data.hpp"

using namespace ROOT::Minuit2;
using namespace std;

class DISFitter : public FCNBase
{
public:
    // MINUIT functions
    double operator() (const vector<double>& par) const; // Calculate Chi^2
    double Up() const {return 1.;}
    
    // Initialize based on MINUIT parameters
    DISFitter(MnUserParameters parameters_);

    
    // Total gamma-p cross section with given polarization and parton
    double  ProtonPhotonCrossSection(const double Qsqr, const double xbj, const Polarization pol,const Parton p, const double mass, FitParameters fitparams) const;
    
    double  ReducedCrossSection(const double Qsqr, const double xbj, const double sqrts, const Parton p, const double mass, FitParameters fitparams) const;
    
    void AddDataset(Data& d);
private:
    IPsat dipole;

    MnUserParameters parameters;
    
    vector<Data*> datasets;

};



#endif
