#ifndef WAVE_FUNCTION_H
#define WAVE_FUNCTION_H

/*
 * General class for wave functions
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2010
 */

// NOTE: Not used at this point in the IPsatfitter, except to define Partons/Polarizations/constants
#include <string>
#include <cmath>

#define LINEINFO __FILE__ << ":" << __LINE__


const int VM_MODE_TOT=1; const int VM_MODE_L=2; const int VM_MODE_T=3;

enum Polarization
{
    TRANSVERSE,
    LONGITUDINAL
};

enum Parton
{
    U,D,S,C,B,LIGHT    // LIGHT = U+D+S
};

typedef double REAL;

class WaveFunction{
    public:
        WaveFunction();
        virtual  REAL PsiSqr_T(REAL Qsqr, REAL r, REAL z) const = 0;
        virtual  REAL PsiSqr_L(REAL Qsqr, REAL r, REAL z) const = 0;
        virtual  REAL PsiSqr_T_intz(REAL Qsqr, REAL r) const = 0;
        virtual  REAL PsiSqr_L_intz(REAL Qsqr, REAL r) const = 0;
        virtual std::string GetParamString()=0;
        REAL PsiSqr_tot(REAL Qsqr, REAL r, REAL z);
        REAL PsiSqr_tot_intz(REAL Qsqr, REAL r);
        REAL PsiSqr_intz(REAL Qsqr, REAL r);
        void SetMode(int m);

        virtual REAL MesonMass();   // Get vector meson mass
    private:
        int mode;   // What to return when PsiSqr_intz
};



// eps(z,Q,r) = sqrt(z(1-z)*Q^2 + m^2)
double epsfunsqr(double z, double Qsqr, double msqr);

double epsfun(double z, double Qsqr, double msqr);

const double ALPHA_e = 1.0/137.035999679;
const double e = sqrt(4.0*M_PI*ALPHA_e);

const int NC = 3;

inline double ABS(double x) { return std::abs(x); }
inline double SQR(double x) { return x*x; }


#endif  // WAVE_FUNCTION_H
