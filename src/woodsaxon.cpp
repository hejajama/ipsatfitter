/*
 * From AmplitudeLib v2 
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011-2017
 */

#include <string>
#include <iostream>
#include <sstream>
#include <cmath>
#include <vector>
#include <gsl/gsl_integration.h>
#include "woodsaxon.hpp"
// Some helpers
#include "wave_function.hpp"
using namespace std;

/*
 * Nuclear density profile
 * Normalization: \int d^2 b T_A(b)=1
 * All lenght units: 1/GeV!
 */

const double FMGEV = 5.068;

double WoodsSaxon::W_S(double r)
{
	// R = 1.12 fm A^(1/3) - 0.86 fm * A^(-1/3)
	double ra = 1.12 * std::pow(A, 1.0/3.0) - 0.86 * std::pow(A, -1.0/3.0);
	double delta = 0.54 * FMGEV;
	ra *= FMGEV;	// fm => 1/GeV

	return w_s_normalization / (std::exp((r - ra)/delta)+1);
}

struct inthelper_ta
{
    WoodsSaxon *WS;
	double b;
	int A;
};

// Return W_S(\sqrt{ b^2 + z^2 } )
double inthelperf_ta(double z, void* p)
{
	inthelper_ta* par = (inthelper_ta*)p;
	return par->WS->W_S(std::sqrt( SQR(par->b) + SQR(z) ));
}

double WoodsSaxon::T_A(double b)
{
    const int INTERVALS = 20;
	gsl_integration_workspace *w = gsl_integration_workspace_alloc(INTERVALS);
	double res, abserr;
    inthelper_ta par; par.b=b; par.WS=this;
	gsl_function f; f.params=&par; f.function=inthelperf_ta;
	int status = gsl_integration_qag(&f, 0, 200, 0, 0.0001, INTERVALS, GSL_INTEG_GAUSS61,
		w, &res, &abserr);
    gsl_integration_workspace_free(w);
	if(status)
		cerr << "T_A integration failed at " << LINEINFO <<", result " << res
			<< ", relerr " << std::abs(res-abserr)/res <<", A=" << A << endl;
	return 2.0*res;	// 2.0 as we integrate z in [0,\infty]
}

// WS initialization, calculates normalization s.t. 
// \int d^3 b WS(b)=1
double inthelperf_ws(double r, void* p)
{
	return r*r*((WoodsSaxon*)p)->W_S(r);
}
	
WoodsSaxon::WoodsSaxon(int A_)
{
    A=A_;
    w_s_normalization = 1.0;
	gsl_integration_workspace *w = gsl_integration_workspace_alloc(10);
	double res, abserr;
	gsl_function f; f.params=this; f.function=inthelperf_ws;
	int status = gsl_integration_qag(&f, 0, 99, 0, 0.001, 10, GSL_INTEG_GAUSS61,
		w, &res, &abserr);
	if (status)
		cerr << "WS normalization integral failed at " << LINEINFO <<", result "
			<< res << " relerr " << std::abs(res-abserr)/res << endl;
	
	res *= 4.0*M_PI;
	
	w_s_normalization = 1.0/res;
    gsl_integration_workspace_free(w);
}

