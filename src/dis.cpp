

#include "dis.hpp"
#include "ipsat.hpp"
#include "data.hpp"
#include "virtual_photon.hpp"
#include <vector>
#include <iostream>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>
#include <Minuit2/MnUserParameterState.h>
#include <iomanip>
#include <sstream>

#include "dglap_cpp/EvolutionLO.h"
#include "dglap_cpp/EvolutionLO_nocoupling.h"
#include "dglap_cpp/AlphaStrong.h"

using namespace std;

const double MINR = 1e-6;
const double MAXR =50*5.068;
// Hera data is very accurate, so eventually one needs to use better accuracy
const double RINTACCURACY = 0.000001;

const int INTEGRATIONDEPTH = 200;


// xg = x*gluon
extern "C"
{
    //double lo_evol_(double *x, double* Q2, double *gluon, int* coupling, double* Ag, double* lambdag  );
    double alphas_(double *mu);
    void initalphas_(int *iord, double *fr2, double *mur, double* asmur, double *mc, double *mb, double* mt);
    void init_(); // Init Mellin momenta
};

string PrintVector(vector<double> v);
double alphas_helper(double asmur, void* p);

// Initialize alphas with given masses such that as(M_z) = experimental value
double InitAlphasMur(FitParameters *par);

/*
 * Calculate chi^2
 * Go trough all datasets
 */




double DISFitter::operator()(const std::vector<double>& par) const
{
    double chisqr = 0;
    
    double light_mass = par[ parameters.Index("light_mass")];
    double charm_mass = par[ parameters.Index("charm_mass")];
    double bottom_mass = par[ parameters.Index("bottom_mass")];
    double lambdag = par[ parameters.Index("lambda_g")];
    double Ag = par[ parameters.Index("A_g")];
    double mu0 =par[parameters.Index("mu_0")];
    double Cscale = par[parameters.Index("C")];
    
    // Force som limits in case MINUIT does not handle these properly
    if (light_mass < 0 or light_mass > 1 or charm_mass < 1.1 or charm_mass > 10
            or lambdag < -10 or lambdag > 10 or Ag < 0 or mu0 < 0
            or Cscale <= 0 or Cscale > 1e6)
        return 9999;

    FitParameters fitparams;
    fitparams.values = &par;
    fitparams.parameter = &parameters;
    
    
    
    // Init alphas
    // Now we take mu_0 as the initial scale of xg
    // InitAlphasMur initializes Alphas() at the initial scale mu_0
    // such that we keep alphas(M_z) = 0.1184
    // Not thread safe!
    EvolutionLO_gluon *cppdglap = NULL;
#pragma omp critical
{
    if (dglapsolver == PIA)
    {
        double asmur = InitAlphasMur(&fitparams);
        fitparams.alphas_mur =asmur;
    }
    else if (dglapsolver == SARTRE)
    {
        dipole.InitializeDGLAP(fitparams);
    }
    else if (dglapsolver == CPPPIA)
    {
        // Initialize alpha_s(M_Z=91.1876)=0.1183
        AlphaStrong *alphas = new AlphaStrong(0, 1.0, 91.1876, 0.1300, charm_mass, bottom_mass, 175);
        // DGLAP_Solver will take care of deleting alphas when it is deleted
        cppdglap = new EvolutionLO_gluon(alphas);
        
        fitparams.cppdglap = cppdglap;
        fitparams.alpha_strong = alphas;
    }
}
    
    
    // Initialize wave functions with current quark masses
    VirtualPhoton wf_lightquark;
    VirtualPhoton wf_charm;
    VirtualPhoton wf_bottom;
    wf_lightquark.SetQuark(LIGHT, light_mass);
    wf_charm.SetQuark(C, charm_mass);
    wf_bottom.SetQuark(B, bottom_mass);
#ifdef USE_INTERPOLATOR
    wf_lightquark.InitializeZintInterpolators();
    wf_charm.InitializeZintInterpolators();
    wf_bottom.InitializeZintInterpolators();
#endif
    
    
    int points=0;
    
    
    int totalpoints = 0;
    for (unsigned int dataset=0; dataset<datasets.size(); dataset++)
        totalpoints += datasets[dataset]->NumOfPoints();
    
    // These loops are trivially parallerizable
    // We only parallerize the inner loop where we have about
    // 250 points (total sigmar) and 50 points (charm)
    for (unsigned int dataset=0; dataset<datasets.size(); dataset++)
    {
#ifdef PARALLEL_CHISQR
    #pragma omp parallel for schedule(dynamic) reduction(+:chisqr) reduction(+:points)
#endif
        for (int i=0; i<datasets[dataset]->NumOfPoints(); i++)
        {
            double x = datasets[dataset]->xbj(i);
            double y = datasets[dataset]->y(i);
            double Q2 = datasets[dataset]->Qsqr(i);
            double sigmar = datasets[dataset]->ReducedCrossSection(i);
            double sigmar_err = datasets[dataset]->ReducedCrossSectionError(i);
            bool onlycharm = datasets[dataset]->OnlyCharm(i);
            
            double sqrts = sqrt( Q2 / (x * y) );
            
            double charmx = x * (1.0 + 4.0*charm_mass*charm_mass / Q2);
            
            double charmx_limitmass = x * (1.0 + 4.0*1.42*1.42 / Q2);
            
            // Include only datapoints where charm contribution can be calculated
            // Use here charmx corresponding to m_c = 1.42 GeV, to keep # of datapoints the same!
            if (charmx_limitmass > 0.01)
                continue;

            double theory_charm = ReducedCrossSection(Q2, charmx, sqrts, &wf_charm, fitparams);
            double theory_light=0, theory_bottom=0;
            if (!onlycharm) // All quarks
            {
                theory_light = ReducedCrossSection(Q2, x, sqrts, &wf_lightquark, fitparams);
                
                // Include bottom conribution if kinematically allowed
                double bottomx = x * (1.0 + 4.0*bottom_mass*bottom_mass/Q2);
                if (bottomx < 0.1) //0.01: chi2 1.399  0.1: 1.409, no b: 1.27
                    theory_bottom = ReducedCrossSection(Q2, bottomx, sqrts, &wf_bottom, fitparams);
                
                
            }
            double theory = theory_light + theory_charm + theory_bottom;
            
            
            if (isnan(theory) or isinf(theory))
            {
                cerr << "Warning: theory result " << theory << " with parameters " << PrintVector(par) << endl;
                theory = 99999999;
            }

            chisqr += datasets[dataset]->Weight(i)*SQR( (theory - sigmar) / sigmar_err );
            points = points + datasets[dataset]->Weight(i);

            // Output for plotting
            
            //cout << setw(10) << x << " " << setw(10)  << Q2 << " " << setw(10) << y << " " << setw(10) << sigmar << " " <<  " " << setw(10)  << sigmar_err << " " << setw(10) << theory_light << " " << setw(10) << theory_charm << " " << setw(10) << theory_bottom << endl;

            
            
        }
    }
    
    cout << "# Calculated chi^2/N = " << chisqr/points << " (N=" << points << "), parameters (" << PrintVector(par) << ")" << endl;
    
    if (dglapsolver == CPPPIA)
    {
        delete cppdglap;
    }
    return chisqr/points;
}




/*
 * Compute total gamma-p cross section
 * by default p=LIGHT which means that we sum over light quarks (u,d,s)
 * If p is something else, use the given quark
 *
 * Default value for the mass is -1 which means that the default values
 * for the quark masses (from VirtualPhoton class) are used
 */

/*
 * Virtual photon-proton cross sections
 * Separately for transversially and longitudinally polarized photons
 */
struct Inthelper_totxs
{
    const IPsat* N;
    Polarization pol;    // L (longitudinal) or T (transverse)
    double Qsqr,xbj;
    const VirtualPhoton* wf;
    FitParameters fitparameters;
    int config;     // To support fluctuations
};

double Inthelperf_totxs(double lnr, void* p)
{
    double r = exp(lnr);
    Inthelper_totxs* par = (Inthelper_totxs*)p;
    
    double result = 0;
    //if (par->fitparameters.values->at( par->fitparameters.parameter->Index("A")) == 1)
        result = r*par->N->DipoleAmplitude_bint(r,par->xbj, par->fitparameters, par->config);
    //else
    //    result = r*par->N->DipoleAmplitude_bint_lumpyA(r,par->xbj, par->fitparameters, par->config);
    
    result *= r;    // Jacobian, as we integrate ln r
    
    if (par->pol==LONGITUDINAL)    // Longitudinal
        result *= par->wf->PsiSqr_L_intz(par->Qsqr, r);
    else if (par->pol==TRANSVERSE)   // Transverse
        result *= par->wf->PsiSqr_T_intz(par->Qsqr, r);
    else
        cerr << "Invalid polarization " << par->pol << " at " << LINEINFO << endl;
    
    return result;
}



 double DISFitter::ProtonPhotonCrossSection(const double Qsqr, const double xbj, const Polarization pol, const VirtualPhoton *wf, FitParameters fitparams) const
{
    Inthelper_totxs par; par.N=&dipole;
    par.wf=wf;
    par.pol=pol; par.Qsqr=Qsqr; par.xbj=xbj;
    par.fitparameters = fitparams;
    par.config = -1;    // no fluctuations
    
    
    gsl_function fun; fun.function=Inthelperf_totxs;
    fun.params=&par;
    
    
    
    double result,abserr;
    gsl_integration_workspace* ws = gsl_integration_workspace_alloc(INTEGRATIONDEPTH);
    int status = gsl_integration_qag(&fun, log(MINR), log(MAXR), 0, RINTACCURACY,
                                     INTEGRATIONDEPTH, GSL_INTEG_GAUSS51, ws, &result, &abserr);
    gsl_integration_workspace_free(ws);
    
    if(status){ std::cerr<< "r integral in ProtonPhotonCrossSection failed with code "
        << status << " (Qsqr=" << Qsqr << ", xbj=" << xbj << " result " << result
        << " relerr=" << abserr/result << "), parameters " << PrintVector(*fitparams.values) <<" at " << LINEINFO << std::endl;
    }
    
    return 2.0*2.0*M_PI*result; //2\pi from \theta integral (angular integral over the dipole orientation)
    // Additional 2 in front is part of gamma^*p cross section
}


/*
 * Reduced cross section from gamma-p cross section
 */
 double DISFitter::ReducedCrossSection(const double Qsqr, const double xbj, const double sqrts, const VirtualPhoton *wf, FitParameters fitparams) const
{
    double kin_y = Qsqr/(sqrts*sqrts*xbj);   // inelasticity, not rapidity
    double xs_l = ProtonPhotonCrossSection(Qsqr, xbj, LONGITUDINAL, wf, fitparams);
    double xs_t = ProtonPhotonCrossSection(Qsqr, xbj, TRANSVERSE, wf, fitparams);
    double f2 = Qsqr/(4.0*SQR(M_PI)*ALPHA_e)*(xs_l+xs_t);
    double fl = Qsqr/(4.0*SQR(M_PI)*ALPHA_e)*xs_l;
    
    return f2 - SQR(kin_y) / ( 1.0 + SQR(1.0-kin_y) ) * fl;  // reduced cross section
}



/*
 * F2
 */
double DISFitter::F2(double Q2, double xbj, FitParameters fitparams ) const
{
    InitAlphasMur(&fitparams);
    VirtualPhoton photon;
    double ml = fitparams.values->at( fitparams.parameter->Index("light_mass"));
    double mc = fitparams.values->at( fitparams.parameter->Index("charm_mass"));
    double mb = fitparams.values->at( fitparams.parameter->Index("bottom_mass"));
    photon.SetQuark(LIGHT, ml );
    double xs_light_l = ProtonPhotonCrossSection(Q2, xbj, LONGITUDINAL, &photon, fitparams);
    double xs_light_t = ProtonPhotonCrossSection(Q2, xbj, TRANSVERSE, &photon, fitparams);
    
    double xc = xbj * (1.0 + 4.0*mc*mc/Q2);
    
    photon.SetQuark(C, mc);
    double xs_charm_l =ProtonPhotonCrossSection(Q2, xc, LONGITUDINAL, &photon, fitparams);
    double xs_charm_t =ProtonPhotonCrossSection(Q2, xc, TRANSVERSE, &photon, fitparams);
    
    double xb = xbj*(1.0 + 4.0*mb*mb/Q2);
    double xs_bottom_l=0;
    double xs_bottom_t = 0;
    if (xb < 0.01  )
    {
        photon.SetQuark(B, mb);
        double xs_bottom_l = ProtonPhotonCrossSection(Q2, xb, LONGITUDINAL, &photon, fitparams);
        double xs_bottom_t = ProtonPhotonCrossSection(Q2, xb, TRANSVERSE, &photon, fitparams);
    }

    return Q2 / (4.0*M_PI*M_PI*ALPHA_e) * (xs_light_l + xs_light_t + xs_charm_l + xs_charm_t + xs_bottom_t + xs_bottom_l);
}

double DISFitter::FL(double Q2, double xbj, FitParameters fitparams ) const
{
    InitAlphasMur(&fitparams);
    VirtualPhoton photon;
    double ml = fitparams.values->at( fitparams.parameter->Index("light_mass"));
    double mc = fitparams.values->at( fitparams.parameter->Index("charm_mass"));
    double mb = fitparams.values->at( fitparams.parameter->Index("bottom_mass"));
    photon.SetQuark(LIGHT, ml );
    double xs_light_l = ProtonPhotonCrossSection(Q2, xbj, LONGITUDINAL, &photon, fitparams);
    
    double xc = xbj * (1.0 + 4.0*mc*mc/Q2);
    
    photon.SetQuark(C, mc);
    double xs_charm_l =ProtonPhotonCrossSection(Q2, xc, LONGITUDINAL, &photon, fitparams);
    
    double xb = xbj*(1.0 + 4.0*mb*mb/Q2);
    double xs_bottom_l=0;
    if (xb < 0.01 )
    {
        photon.SetQuark(B, mb);
        double xs_bottom_l = ProtonPhotonCrossSection(Q2, xb, LONGITUDINAL, &photon, fitparams);
    }
    
    return Q2 / (4.0*M_PI*M_PI*ALPHA_e) * (xs_light_l + xs_charm_l + xs_bottom_l);
}




void DISFitter::AddDataset(Data& d)
{
    datasets.push_back(&d);
}


DISFitter::DISFitter(MnUserParameters parameters_)
{
    parameters = parameters_;
    dipole.SetSaturation(true);
    
    int A = parameters.Value("A");
    if (A>1)
        dipole.InitNucleus(A);
    
    
    // Init dglap
    dglapsolver = CPPPIA;
    dipole.SetDGLAPSolver(CPPPIA);

    
    
}

string PrintVector(vector<double> v)
{
    stringstream ss;
    for (int i=0; i <v.size(); i++)
    {
        ss << v[i];
        if (i < v.size()-1)
            ss <<", ";
    }
    return ss.str();
}


//////////// ALPHA_S INITIALIZATION
// Used when Fortran DGLAP solver is selected


double InitAlphasMur(FitParameters *par)
{
    // Init Alpha_s
    // Find alphas_mur
    gsl_function f; f.function = &alphas_helper;
    f.params = par;
    
    const gsl_root_fsolver_type *T = gsl_root_fsolver_brent;
    gsl_root_fsolver *s = gsl_root_fsolver_alloc (T);
    double min = 0.1; double max = 0.5;
    gsl_root_fsolver_set (s, &f, min, max);
    
    
    int iter=0; int maxiter = 100;
    double asmur = 0.4; int status;
    do
    {
        iter++;
        gsl_root_fsolver_iterate (s);
        asmur = gsl_root_fsolver_root (s);
        min = gsl_root_fsolver_x_lower (s);
        max = gsl_root_fsolver_x_upper (s);
        status = gsl_root_test_interval (min, max,   0, 0.001);
        
    }
    while (status == GSL_CONTINUE && iter < maxiter);
    
    gsl_root_fsolver_free(s);
    
    if (iter >= maxiter)
        cerr << "Initializing alphas with parameters " << PrintVector(*par->values) << " did not succeed";
    //cout << "--- initialized alphas at " << par->values->at( par->parameter->Index("mu_0")) << " GeV to " << asmur << endl;
    
    
    
    return asmur;

}
// Helper function to solve asmur s.t. alphas(M_z)=0.11884
double alphas_helper(double asmur, void* p)
{
    //cout << "Trying " << asmur << endl;
    // Init alphas
    FitParameters *par = (FitParameters*)p;
    int iord=0;
    double fr2 = 1.0;
    double mur = ALPHAS_MUR;
    double mc = par->values->at( par->parameter->Index("charm_mass"));
    double mb = par->values->at( par->parameter->Index("bottom_mass"));
    double mt = 175;
    if (mur > mc)
    {
        cerr << "mu_r > m_c at alphas_helperf(), parameters: " << PrintVector(*par->values) << endl;
        return 100;
    }
    
    initalphas_(&iord, &fr2, &mur, &asmur, &mc, &mb, &mt);
    double mz=91.1876;
    return alphas_(&mz) - 0.1183; // From HERA II 1506.06042
}
