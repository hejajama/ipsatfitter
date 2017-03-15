

#include "dis.hpp"
#include "ipsat.hpp"
#include "data.hpp"
#include "virtual_photon.hpp"
#include <vector>
#include <iostream>
#include <gsl/gsl_integration.h>
#include <Minuit2/MnUserParameterState.h>
#include <sstream>

using namespace std;

const double MINR = 1e-6;
const double MAXR = 50;
const double RINTACCURACY = 0.005;

const int INTEGRATIONDEPTH = 50;


// LO DGLAP solver
//SUBROUTINE LO_evol(X, Q2, gluon, coupling, Ag, lambdag)
// xg = x*gluon
extern "C"
{
    double lo_evol_(double *x, double* Q2, double *gluon, int* coupling, double* Ag, double* lambdag  );
    double alphas_(double mu);
    void init_(); // Init Mellin momenta
};


/*
 * Calculate chi^2
 * Go trough all datasets
 */

string PrintVector(vector<double> v) {
    stringstream ss; for (int i=0; i <v.size(); i++) ss << v[i] << ", ";
    return ss.str();
}

double DISFitter::operator()(const std::vector<double>& par) const
{
    double chisqr = 0;
    
    double light_mass = par[ parameters.Index("light_mass")];
    double heavy_mass = par[ parameters.Index("heavy_mass")];
    double lambdag = par[ parameters.Index("lambda_g")];
    double Ag = par[ parameters.Index("A_g")];
    
    // Init dglap
    init_();
    // evaluate xg once to call necessary initialization
    double x=0.01; double q2=10; double gluon=0; int coupling=0;
    lo_evol_(&x,&q2,&gluon,&coupling,&Ag,&lambdag);
    //cout << "as * xg = " << gluon << endl;
    
    // If quark masses do not change, then this is not necessary.
    // However, now this possible duplication of work keeps the code
    // thread-safe

    VirtualPhoton wf_lightquark;
    VirtualPhoton wf_heavyquark;
    wf_lightquark.SetQuark(LIGHT, light_mass);
    wf_heavyquark.SetQuark(C, heavy_mass);
#ifdef USE_INTERPOLATOR
    wf_lightquark.InitializeZintInterpolators();
    wf_heavyquark.InitializeZintInterpolators();
#endif
    
    FitParameters fitparams;
    fitparams.values = &par;
    fitparams.parameter = &parameters;
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
    #pragma omp parallel for reduction(+:chisqr) reduction(+:points)
#endif
        for (int i=0; i<datasets[dataset]->NumOfPoints(); i++)
        {
            double x = datasets[dataset]->xbj(i);
            double y = datasets[dataset]->y(i);
            double Q2 = datasets[dataset]->Qsqr(i);
            double sigmar = datasets[dataset]->ReducedCrossSection(i);
            double sigmar_err = datasets[dataset]->ReducedCrossSectionError(i);
            bool onlycharm = datasets[dataset]->OnlyCharm();
            
            double sqrts = sqrt( Q2 / (x * y) );
            
            double charmx = x * (1.0 + 4.0*heavy_mass*heavy_mass / Q2);
            
            if (charmx > 0.01)
                continue;

            double theory_charm = ReducedCrossSection(Q2, charmx, sqrts, &wf_heavyquark, fitparams);
            double theory_light;
            if (onlycharm)
                theory_light = 0;
            else
                theory_light = ReducedCrossSection(Q2, x, sqrts, &wf_lightquark, fitparams);
            
            double theory = theory_light + theory_charm;
            
            
            

            chisqr += datasets[dataset]->Weight()*SQR( (theory - sigmar) / sigmar_err );
            points = points + datasets[dataset]->Weight();

/* #pragma omp critical
            {
                cout << x << " " << Q2 << " theory " << theory << " exp " << sigmar << " +/- " << sigmar_err << endl;
            if (points % 50 == 0)
                cout << "Done approximately " << points << " / " << totalpoints << endl;
            
            }
    */
            
            
        }
    }
    
    cout << "Calculated chi^2/N = " << chisqr/points << " (N=" << points << "), parameters " << PrintVector(par) << endl;
    
    exit(1);
    return chisqr;
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
    
    double result = r*par->N->DipoleAmplitude_bint(r,par->xbj, par->fitparameters, par->config);
    
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
    //int status = gsl_integration_qng(&fun, MinR(), MaxR(),
    //0, 0.001,  &result, &abserr, &eval);
    
    if(status){ std::cerr<< "r integral in ProtonPhotonCrossSection failed with code "
        << status << " (Qsqr=" << Qsqr << ", xbj=" << xbj << " result " << result
        << " relerr=" << abserr/result << ") at " << LINEINFO << std::endl;
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





void DISFitter::AddDataset(Data& d)
{
    datasets.push_back(&d);
}


DISFitter::DISFitter(MnUserParameters parameters_)
{
    parameters = parameters_;
    
    
}

