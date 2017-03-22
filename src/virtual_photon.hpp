#ifndef VirtualPhoton_H
#define VirtualPhoton_H

/* AmplitudeLib
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2010-2016
 *
 * Difference to AmlitudeLib: By default charm is included, and light quark masses are set to ~0
 */

#include "wave_function.hpp"

#include <iostream>
#include <string>
#include <vector>

// Interpolate zintegrals
// This is a flag, as in some older GSL installations the necessary 2D interpolation
// functions are not available, in that case one can just comment out this
// Also, removing this allows one to easily check if interpolation works
// (results should not change)
//#define USE_INTERPOLATOR


#ifdef USE_INTERPOLATOR
    #include "interpolation2d.hpp"
#endif



/**
 * \class VirtualPhoton virtual_photon.hpp amplitudelib/virtual_photon.hpp
 * Virtual photon-q bar q overlap
 * 
 * Ref: Kowalski, Motyka and Watt, see arXiv: hep-ph/0606272v2
 * 
 * By default uses 3 light quarks (u,d,s)
 */
class VirtualPhoton  {
    public:
        VirtualPhoton();
        ~VirtualPhoton();
    
        /**
         * Overlap between q bar
         * q and transverse photon
         *
         * @param Qsqr photon virtuality [GeV^2]
         * @param r dipole size
         * @param z longitudinal momentum fraction of the quark
         */
        const double PsiSqr_T(double Qsqr, double r, double z) const;

        /**
         * Overlap between q bar q and longitudinal photon
         *
         * @param Qsqr photon virtuality [GeV^2]
         * @param r dipole size
         * @param z longitudinal momentum fraction of the quark
         */
        const double PsiSqr_L(double Qsqr, double r, double z) const;
        
        // Overlap wave functions integrated over z=[0,1]
        /**
         * Overlap between q bar q and transverse photon integrated over z
         *
         * @param Qsqr photon virtuality [GeV^2]
         * @param r dipole size
         */
        const double PsiSqr_T_intz(double Qsqr, double r) const;

        /**
         * Overlap between q bar q and longitudinal photon integrated over z
         *
         * @param Qsqr photon virtuality [GeV^2]
         * @param r dipole size
         */
        const double PsiSqr_L_intz(double Qsqr, double r) const;
        
        
        
        std::string GetParamString();
        
        /**
         * Select quar flavor
         *
         * If not set, sums over light quarks (u,d,s) with default mass 0.14 GeV.
         *
         * @param p quark flavor
         * @param mass quark mass, optional, default values used if not given
         */
        void SetQuark(Parton p, double mass=-1);
    
    
#ifdef USE_INTERPOLATOR
        /**
         * Initialize zint interpolations, note that these depent on used quark flavors and masses!
         *
         */
        void InitializeZintInterpolators();
#endif
    
    private:
        // Parameters
        std::vector<double> e_f;   // Quark charges
        std::vector<double> m_f;   // Quark masses (GeV)

        double Epsilon(double Qsqr, double z, int f) const;
    
#ifdef USE_INTERPOLATOR
        Interpolator2D transverse_zint_interpolator;
        Interpolator2D longitudinal_zint_interpolator;
        bool interpolator_ready;
        double interpolator_maxQ2;
        double interpolator_minQ2;
        double interpolator_maxr;
        double interpolator_minr;
#endif
    
        
};

std::ostream& operator<<(std::ostream& os, VirtualPhoton& ic);



#endif
