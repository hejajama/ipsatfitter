/*
 * Virtual photon wave function
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2010-2016
 */
 
#include "virtual_photon.hpp"
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

using namespace std;

const int Nc=3;

const double ZINTACCURACY=0.001;
const int MAXITER_ZINT=500;
const double MINZ=1e-8;  // Integration limits
const double MAXZ=1.0-MINZ;

VirtualPhoton::VirtualPhoton()
{
	// Initialize with light quarks
	// Parameters from IPsat2012 fit
    // Quark charges, 0=u, 1=d, 2=s, 3=c
    e_f.push_back(2.0/3.0); e_f.push_back(-1.0/3.0); e_f.push_back(-1.0/3.0); e_f.push_back(2.0/3.0);
    m_f.push_back(0.01); m_f.push_back(0.01); m_f.push_back(0.01); m_f.push_back(1.4);
    
#ifdef USE_INTERPOLATOR
    interpolator_ready = false;
#endif
 
}

VirtualPhoton::~VirtualPhoton()
{
}
/*
 * Transversially polarized component of the overlap
 */



const double VirtualPhoton::PsiSqr_T(double Qsqr, double r, double z) const
{
    double result=0;
    
    // Optimize unnecessary K_0 evaluations if masses are the same
    double saved_mass = -1;
    double saved_K0 = -1;
    double saved_K1 = -1;
    for (unsigned int f=0; f<e_f.size(); f++)     // Sum quark flavors
    {
        
        double epstmp=Epsilon(Qsqr,z,f);
        
        gsl_sf_result result_k1;
        gsl_sf_result result_k0;
        double bessel_k1; double bessel_k0;
        
        if ( abs(m_f[f] - saved_mass) < 0.001 ) // Same mass as on the previous round
        {
            bessel_k0 = saved_K0;
            bessel_k1 = saved_K1;
        }
        else
        {
            if (gsl_sf_bessel_K1_e(epstmp*r, &result_k1) == GSL_EUNDRFLW)
                bessel_k1 = 0;
            else
                bessel_k1 = result_k1.val;
            
            if (gsl_sf_bessel_K0_e(epstmp*r, &result_k0) == GSL_EUNDRFLW)
                bessel_k0 = 0;
            else
                bessel_k0 = result_k0.val;
        }
        
        result += SQR(e_f[f])*(
            (SQR(z)+SQR(1.0-z))*SQR(epstmp)
                    *SQR(bessel_k1)
                + SQR(m_f[f]*bessel_k0)
            );
        
        saved_K0 = bessel_k0;
        saved_K1 = bessel_k1;
        saved_mass = m_f[f];
    }
    result *= Nc/(2.0*SQR(M_PI))*ALPHA_e;

    return result;
}

/*
 * Longitudinally polarized component of the overlap
 */

const double VirtualPhoton::PsiSqr_L(double Qsqr, double r, double z) const
{
    // Avoid unnecessary special function evaluation
    double saved_mass = -1;
    double saved_K0 = -1;
    
    double result=0;
    for (unsigned int f=0; f<e_f.size(); f++)     // Sum quark flavors
    {
        double bessel_k0;
        if ( abs(m_f[f] - saved_mass) < 0.0001)
            bessel_k0 = saved_K0;
        else
        {
            double epstmp=Epsilon(Qsqr,z,f);
            gsl_sf_result result_k0;
            
            
            if (gsl_sf_bessel_K0_e(epstmp*r, &result_k0) == GSL_EUNDRFLW)
                bessel_k0 = 0;
            else
                bessel_k0 = result_k0.val;
        }

        
        result += SQR(e_f[f])* SQR( bessel_k0 );
        
        saved_mass = m_f[f];
        saved_K0 = bessel_k0;
    }
    result *= 2.0*Nc/SQR(M_PI)*ALPHA_e*Qsqr*SQR(z)*SQR(1.0-z);

    return result;
}

/* 
 * Wave function overlap integrated over z=[0,1]
 * PsiSqr_T/L is quite a smooth function so there is 
 * nothing tricky to do here, just use gsl_integration_qng
 */
 
/* As we have to integrate a member function of this class by GSL,
 * we need some helper structures
 */
struct zinthelper{
    const VirtualPhoton *vm_p;
    double  Qsqr;
    double  r;
};

double  zhelperfuncT(double z, void * p){
  return ((zinthelper*)p)->vm_p->PsiSqr_T(
							((zinthelper*)p)->Qsqr,
							((zinthelper*)p)->r,
							z);
}

double  zhelperfuncL(double z, void * p){
  return ((zinthelper*)p)->vm_p->PsiSqr_L(
							((zinthelper*)p)->Qsqr,
							((zinthelper*)p)->r,
							z);
}

const double VirtualPhoton::PsiSqr_T_intz(double Qsqr, double r) const
{
#ifdef USE_INTERPOLATOR
    if (interpolator_ready and Qsqr < interpolator_maxQ2 and Qsqr > interpolator_minQ2 and r < interpolator_maxr and r  > interpolator_minr)
    {
        double interpolated = transverse_zint_interpolator.Evaluate(log(Qsqr), log(r));
        return interpolated;
        
    }
#endif
    double result,abserr;
    struct zinthelper zintpar;
    zintpar.vm_p=this;
    zintpar.Qsqr=Qsqr;
    zintpar.r=r;
    gsl_function int_helper;
    int_helper.function=&zhelperfuncT;
    int_helper.params=&zintpar;
    
    gsl_integration_workspace* ws = gsl_integration_workspace_alloc(MAXITER_ZINT);
    int status = gsl_integration_qag(&int_helper, MINZ, MAXZ, 0, ZINTACCURACY,
        MAXITER_ZINT, GSL_INTEG_GAUSS51, ws, &result, &abserr);
    gsl_integration_workspace_free(ws);

    if(status){ std::cerr<< "z integral in Photon failed with code " 
        << status << " (transverse, Qsqr=" << Qsqr << ", r=" << r 
        << "relerr=" << abserr/result << ") at " << LINEINFO << std::endl;}
  

    return result;
}

const double VirtualPhoton::PsiSqr_L_intz(double Qsqr, double r) const
{
    double result,abserr;
    struct zinthelper zintpar;
    zintpar.vm_p=this;
    zintpar.Qsqr=Qsqr;
    zintpar.r=r;
    gsl_function int_helper;
    int_helper.function=&zhelperfuncL;
    int_helper.params=&zintpar;
    
    gsl_integration_workspace* ws = gsl_integration_workspace_alloc(MAXITER_ZINT);
    int status = gsl_integration_qag(&int_helper, MINZ, MAXZ, 0, ZINTACCURACY,
        MAXITER_ZINT, GSL_INTEG_GAUSS51, ws, &result, &abserr);
    gsl_integration_workspace_free(ws);
    
    if(status){ std::cerr<< "z integral in VirtualPhoton failed: code " 
        << status << " (longitudinal, Qsqr=" << Qsqr << ", r=" << r 
        << "relerr=" << abserr/result << ") at " << LINEINFO << std::endl;}

    return result;
}



double VirtualPhoton::Epsilon(double Qsqr, double z, int f) const
{
    return std::sqrt( z*(1.0-z)*Qsqr+SQR(m_f[f]) );
}


/*
 * Change wave function to the specific quark, supported quarks are u,d,s,c,b
 * Note: Only this quark is used once this function is called
 * 
 * A good way to check that everything works is to compute e.g. F2 using the
 * default quark content (u,d,s) and these quarks separately
 * 
 * If mass<0 (default value), then use standard literature value for the masses
 */

void VirtualPhoton::SetQuark(Parton p, double mass)
{
	// Clear 
	e_f.clear();
	m_f.clear();
	double m;
	switch(p)
	{
        case LIGHT:
            m=0.14;
            if (mass>=0)
                m=mass;
            m_f.push_back(m);
            e_f.push_back(2.0/3.0);
            m_f.push_back(m);
			e_f.push_back(-1.0/3.0);
            m_f.push_back(m);
			e_f.push_back(-1.0/3.0);
            break;
		case U:
			m_f.push_back(0.14);
            if (mass>=0) m_f[0]=mass;
			e_f.push_back(2.0/3.0);
			break;
		case D:
			m_f.push_back(0.14);
            if (mass>=0) m_f[0]=mass;
			e_f.push_back(-1.0/3.0);
			break;
		case S:
			m_f.push_back(0.14);
            if (mass>=0) m_f[0]=mass;
			e_f.push_back(-1.0/3.0);
			break;
		case C:
			m_f.push_back(1.27);
            if (mass>=0) m_f[0]=mass;
			e_f.push_back(2.0/3.0);
			break;
		case B:
			m_f.push_back(4.2);
            if (mass>=0) m_f[0]=mass;
			e_f.push_back(-1.0/3.0);
			break;
		default:
			cerr << "Unknown parton " << p << " at " << LINEINFO << endl;
	}
    
    // When quark content is changed, interpolator can not be ready!
#ifdef USE_INTERPOLATOR
    interpolator_ready = false;
#endif

}

std::string VirtualPhoton::GetParamString()
{
    std::stringstream str;
    for (unsigned int f=0; f<e_f.size(); f++)
		str << "e_f[" << f << "]=" << e_f[f]  << ", m_f[" << f <<"]=" << m_f[f] << " ";
    return str.str();
}

std::ostream& operator<<(std::ostream& os, VirtualPhoton& ic)
{
    return os << " Virtual photon wave function. Params: "
        << ic.GetParamString() << " .";
        
}

#ifdef USE_INTERPOLATOR

/*
 * Initialize zint grid
 */
void VirtualPhoton::InitializeZintInterpolators()
{
    double minlnQ2 = log(0.1);
    double maxlnQ2 = log(200);
    double lnQ2step = 0.1;
    double minlnr = log(1e-5);
    double maxlnr = log(1e2);
    double lnrstep = 0.1;
    
    vector<double> lnQ2vals; vector<double> lnrvals;
    vector<double> zint_t; vector<double> zint_l;
    for (double lnQ2=minlnQ2; lnQ2<=maxlnQ2; lnQ2+=lnQ2step) lnQ2vals.push_back(lnQ2);
    for (double lnr=minlnr; lnr<=maxlnr; lnr+=lnrstep) lnrvals.push_back(lnr);
    
    // Fill datagrid, Q2=x, r=y
    for (double lnr=minlnr; lnr<=maxlnr; lnr+=lnrstep)
    {
        for (double lnQ2=minlnQ2; lnQ2<=maxlnQ2; lnQ2+=lnQ2step)
        {
        
            double r = exp(lnr);
            double Q2 = exp(lnQ2);
            zint_t.push_back(PsiSqr_T_intz(Q2,r));
            zint_l.push_back(PsiSqr_L_intz(Q2,r));
            
        }
    }
    transverse_zint_interpolator.Initialize(lnQ2vals, lnrvals, zint_t);
    longitudinal_zint_interpolator.Initialize(lnQ2vals, lnrvals, zint_l);
    
    interpolator_maxr = exp(lnrvals[lnrvals.size()-1]);
    interpolator_minr = exp(lnrvals[0]);
    interpolator_maxQ2 = exp(lnQ2vals[lnQ2vals.size()-1]);
    interpolator_minQ2 = exp(lnQ2vals[0]);
    
    interpolator_ready = true;
}

#endif

