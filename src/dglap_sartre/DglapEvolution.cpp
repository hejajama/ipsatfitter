//==============================================================================
//  DglapEvolution.cpp
//
//  Copyright (C) 2010-2013 Tobias Toll and Thomas Ullrich
//
//  This file is part of Sartre version: 1.1
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation.
//  This program is distributed in the hope that it will be useful,
//  but without any warranty; without even the implied warranty of
//  merchantability or fitness for a particular purpose. See the
//  GNU General Public License for more details.
//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//  Author: Tobias Toll
//  Last update:
//  $Date: 2013-05-30 01:56:19 +0530 (Thu, 30 May 2013) $
//  $Author: thomas.ullrich@bnl.gov $
//==============================================================================
#include "DglapEvolution.h"
//#include "TableGeneratorSettings.h"
#include <iostream>
#include <cmath>
#include "Parameters_bSat.h"
//#include "Parameters_bNonSat.h"
#include "AlphaStrong.h"
//extern "C" {
//#include "laguerre.h"
//#include "dglap.h"
//}

using namespace std;

DglapEvolution* DglapEvolution::mInstance = 0;
double DglapEvolution::mAg = 0;
double DglapEvolution::mLambdaG = 0;
double DglapEvolution::mMu02 = 0;
double DglapEvolution::mDerivative = 0;
double DglapEvolution::mCharmThreshold = 0;
double DglapEvolution::mBeautyThreshold = 0;
grid_dist* DglapEvolution::mDist3 = 0;
grid_dist* DglapEvolution::mDist4 = 0;
grid_dist* DglapEvolution::mDist5 = 0;

DglapEvolution::DglapEvolution()
{
    mAg = bSat_Ag;
    mLambdaG = bSat_lambdaG;
    mMu02 = bSat_Mu02;
    mCharmThreshold=quarkMass[3]*quarkMass[3];
    mBeautyThreshold=quarkMass[4]*quarkMass[4];
}

DglapEvolution& DglapEvolution::instance()
{
    if(!mInstance) {
        mInstance = new DglapEvolution;
    }
    return *mInstance;
}

DglapEvolution::~DglapEvolution()
{
    if( mInstance ) {
        delete mInstance;
        mInstance = 0;
    }
}

void DglapEvolution::setS(double val)
{
    mS=val;
}

double DglapEvolution::qns0(double)
{
    return 0.0; // no light quarks
}

double DglapEvolution::qs0(double)
{
    return 0.0; // no strange quarks
}

double DglapEvolution::g03(double x)
{
    return mAg*pow(x,-mLambdaG)*pow((1-x),5.6);  //
}

double DglapEvolution::g04(double x)
{
    if(x<=mDist3->xmin) return 0;
    value_dist val = PDF_interpolated(x, mCharmThreshold, mDist3);
    return val.g;
}

double DglapEvolution::g05(double x)
{
    if(x<=mDist4->xmin) return 0;
    value_dist val = PDF_interpolated(x, mBeautyThreshold, mDist4);
    return val.g;
}

double DglapEvolution::G(double x, double Q2)
{
    //
    //  To speed up thing during all these integrations
    //  we keep the last value. Helps a bit.
    //
    static bool init = true;
    
    int N=100;         // maximum Laguerre order
    
    double xmin=1.0e-8;
    //Create one table for 3, 4, 5 number of flavours respectively
    int Nf3=3;
    double q2i3=mMu02;  // initial value of Q^2
    double q2f3=mS;//mCharmThreshold;
    tab_Pij *tt_pij3;
    tab_dist *dist03;
    evol_pair pair3;
    tab_evol *tt_evol3;
    //////////////////////////////////////////
    int Nf4=4;
    double q2i4=mCharmThreshold;
    if(mMu02>mCharmThreshold) q2i4=mMu02;
    double q2f4=mS;//mBeautyThreshold;
    tab_Pij *tt_pij4;
    tab_dist *dist04;
    evol_pair pair4;
    tab_evol *tt_evol4;
    ///////////////////////////////////////////
    int Nf5=5;
    double q2i5=mBeautyThreshold;
    if(mMu02>mBeautyThreshold) q2i5=mMu02;
    double q2f5=mS;     // max/final value of Q^2
    tab_Pij *tt_pij5;
    tab_dist *dist05;
    evol_pair pair5;
    tab_evol *tt_evol5;
    ////////////////////////////////////////////
    
    if (init) {
        // The procedure:
        // Create separate PDFs for each Nf=3, 4, 5
        // For each mass threshold, the previous PDF serves as starting distribution for the next.
        //
        // In order to get the coupling at thresholds, create a temporary instance of alpha_s,
        // then match the coupling at thresholds using this.
        //
        AlphaStrong as;
        as.init(bSat_coupling, bSat_asOrder);
        
        double as_at_charmThreshold=as.alphaS(mCharmThreshold);
        double as_at_beautyThreshold=as.alphaS(mBeautyThreshold);
        
        //////////////////////////////////////////////////////////////
        if(mMu02<mCharmThreshold){
            cout<<endl<< "DglapEvolution::G(): initializing DGLAP evolution engine for Nf=3" << endl;
            tt_pij3=create_Lag_Pij_table(Nf3,N);
            dist03=Lag_dist(N,qns0,qs0,g03);
            dist03->q2=q2i3;
            
            set_LO();
            set_NO_QUARKS();
            compute_GLUON_DERIVATIVE(); // computes  val.dg = d(xG)/dx
            
            pair3=DGLAP_evol(dist03,q2f3,4000,tt_pij3, mCharmThreshold, as_at_charmThreshold);
            tt_evol3=pair3.pdf;
            mDist3=PDF_grid_calc(2000,xmin,tt_evol3);
        }
        ///////////////////////////////////////////////////////////////
        if(mMu02<mBeautyThreshold){
            cout<<endl<< "DglapEvolution::G(): initializing DGLAP evolution engine for Nf=4" << endl;
            tt_pij4=create_Lag_Pij_table(Nf4,N);
            if(mMu02<mCharmThreshold)
            dist04=Lag_dist(N,qns0,qs0,g04);
            else
            dist04=Lag_dist(N,qns0,qs0,g03);
            dist04->q2=q2i4;
            
            set_LO();
            set_NO_QUARKS();
            compute_GLUON_DERIVATIVE(); // computes  val.dg = d(xG)/dx
            
            pair4=DGLAP_evol(dist04,q2f4,4000,tt_pij4, mBeautyThreshold, as_at_beautyThreshold);
            tt_evol4=pair4.pdf;
            mDist4=PDF_grid_calc(2000,xmin,tt_evol4);
        }
        ///////////////////////////////////////////////////////////////
        cout<<endl<< "DglapEvolution::G(): initializing DGLAP evolution engine for Nf=5" << endl;
        tt_pij5=create_Lag_Pij_table(Nf5,N);
        if(mMu02<mBeautyThreshold)
        dist05=Lag_dist(N,qns0,qs0,g05);
        else
        dist05=Lag_dist(N,qns0,qs0,g03);
        dist05->q2=q2i5;
        
        set_LO();
        set_NO_QUARKS();
        compute_GLUON_DERIVATIVE(); // computes  val.dg = d(xG)/dx
        
        double MZ2=91.188*91.188;
        pair5=DGLAP_evol(dist05,q2f5,4000,tt_pij5,MZ2, bSat_coupling);
        tt_evol5=pair5.pdf;
        mDist5=PDF_grid_calc(2000,xmin,tt_evol5);
        
        init = false;
    }
    
    double result=0;
    if(Q2<mCharmThreshold){
        value_dist val = PDF_interpolated(x, Q2, mDist3);
        result = val.g;  // x*G
        mDerivative = val.dg;   // keep
    }
    else if(Q2<mBeautyThreshold and Q2>=mCharmThreshold){
        value_dist val = PDF_interpolated(x, Q2, mDist4);
        result = val.g;  // x*G  
        mDerivative = val.dg;   // keep
    }
    else if(Q2>=mBeautyThreshold){
        value_dist val = PDF_interpolated(x, Q2, mDist5);  
        result = val.g;  // x*G  
        mDerivative = val.dg;   // keep
    }
    
    result /= x;            // G  
    
    
    return result;   
}  

double DglapEvolution::dxGdxOfLastCall() const 
{
    return mDerivative;
}


