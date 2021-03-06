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

// Parallerize chi^2 calculation
// DISFitter is designed to be thread safe. Note that also
// VirtualPhoton.PsiSqr_intz and IPsat dipole amplitude evaluations must be
// too.
// Note that MINUIT also runs in parallel (# of threads depends on the number of
// parameters 
#define PARALLEL_CHISQR



/*
 *
 * List of parameter strings:
 * light_mass       light quark mass in GeV
 * charm_mass       charm quark mass in GeV
 * bottom_mass      b quark mass in GEV
 *
 * B_G              proton size in GeV^(-2)
 *
 * lambda_g         x^(-lambda_g) term in the gluon distribution at the initial condition
 * A_g              Normalization of the gluon distribution in DGLAP solver
 * mu_0             Parametrizes scale, mu^2 = mu_0^2 + 4/r^2, in GeV
 * 
 */

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

    
    // Total gamma-p cross section with given wave function (describing polarization and parton)
    double  ProtonPhotonCrossSection(const double Qsqr, const double xbj, const Polarization pol,const VirtualPhoton* wf , FitParameters fitparams) const;
    
    double  ReducedCrossSection(const double Qsqr, const double xbj, const double sqrts, const VirtualPhoton* wf, FitParameters fitparams) const;
    
    double F2(double Q2, double x, FitParameters fitparams) const;
    double FL(double Q2, double x, FitParameters fitparams) const;
    
    void AddDataset(Data& d);
    
    void SetSaturation(bool s) { dipole.SetSaturation(s); }
    void SetSinglet(bool s) { dipole.SetSinglet(s); }
    
    IPsat &GetDipole() { return dipole; }
    
    DGLAP_Solver GetDGLAPSolver() { return dglapsolver; }
    void SetDGLAPSolver(DGLAP_Solver s) { dglapsolver = s; dipole.SetDGLAPSolver(s); }

    double MAXR;
private:
    IPsat dipole;

    MnUserParameters parameters;
    
    DGLAP_Solver dglapsolver;
    
    vector<Data*> datasets;
    

};

const double ALPHAS_MUR = 1;    // Scale at which we initialize alpha_s [GeV]



#endif
