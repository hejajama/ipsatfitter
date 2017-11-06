
/*
 * Solve DGLAP equation with given parametrization
 * and evaluate dipole amplitude
 *
 * Uses exactly the same DGLAP solver (LO_evolution_routine.f and
 * alphaS.f) as is used in the fit code
 *
 * H. Mäntysaari and P. Zurita, 2017
 */

#include "dipoleamplitude.hpp"
#include <string>
#include <cmath>
#include <iostream>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>
#include <sstream>
#include <iomanip>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_bessel.h>

struct inthelper_xg
{
    DipoleAmplitude* dipole;
    double x;
    double Q2; // Upper limit, external
    double q;
    double r;
    double b;
};
const double INTACCURCY = 0.05;
const int INTDEPTH = 20;
double inthelperf_b(double b, void* p);
double inthelperf_q(double q2, void* p);
double inthelperf_r(double r, void* p);

double xg_from_ipsat(double x, double Q2, DipoleAmplitude& dipole)
{
    gsl_function fun; fun.function=inthelperf_r;
    inthelper_xg par;
    par.x=x; par.Q2=Q2;
    par.dipole = &dipole;
    fun.params=&par;
    
    double result,abserr;
    gsl_integration_workspace* ws = gsl_integration_workspace_alloc(INTDEPTH);
    int status = gsl_integration_qag(&fun, 0, 99, 0, INTACCURCY,
                                     INTDEPTH, GSL_INTEG_GAUSS21, ws, &result, &abserr);
    if (status)
        cerr << "bintegral failed at Q2=" << Q2 <<", result " << result << " relerror " << abserr/result << endl;
    gsl_integration_workspace_free(ws);
    
    return result;
}

double inthelperf_b(double b, void* p)
{
    inthelper_xg *par = (inthelper_xg*) p;
    par->b=b;
    gsl_function fun; fun.function=inthelperf_r;
    fun.params=par;
    double result,abserr;
    gsl_integration_workspace* ws = gsl_integration_workspace_alloc(INTDEPTH);
    int status = gsl_integration_qag(&fun, 0, 99, 0, INTACCURCY,
                                     INTDEPTH, GSL_INTEG_GAUSS21, ws, &result, &abserr);
    if (status)
        cerr << "qintegral failed at b=" << b <<", result " << result << " relerror " << abserr/result << endl;
    gsl_integration_workspace_free(ws);
    
    return result;
}

double inthelperf_r(double r, void* p)
{
    inthelper_xg *par = (inthelper_xg*) p;
    par->r=r;
    gsl_function fun; fun.function=inthelperf_q;
    fun.params=par;
    double result,abserr;
    gsl_integration_workspace* ws = gsl_integration_workspace_alloc(INTDEPTH);
    int status = gsl_integration_qag(&fun, 0, par->Q2, 0, INTACCURCY,
                                     INTDEPTH, GSL_INTEG_GAUSS21, ws, &result, &abserr);
    if (status)
        cerr << "qintegral failed at r=" << r <<", result " << result << " relerror " << abserr/result << endl;
    gsl_integration_workspace_free(ws);
    
    return result;
}

double inthelperf_q(double q2_internal, void* p)
{
    inthelper_xg *par = (inthelper_xg*) p;
    
    // dq^2 q^2 dr r J_0(k*r) S(r)
    return q2_internal *  par->r * gsl_sf_bessel_J0(sqrt(q2_internal) * par->r) * (1.0 - par->dipole->N_bint(par->r, par->x));
}


