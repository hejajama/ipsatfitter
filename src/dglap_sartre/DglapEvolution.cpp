//==============================================================================
//  DglapEvolution.cpp
//
//  Copyright (C) 2010-2016 Tobias Toll and Thomas Ullrich 
//
//  This file is part of Sartre. 
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
//  $Date: 2016-11-01 13:24:45 -0400 (Tue, 01 Nov 2016) $
//  $Author: ullrich $
//==============================================================================
#include "Dglap.h"
#include "DglapEvolution.h"
#include "TableGeneratorSettings.h"  
#include "DipoleModelParameters.h"
#include <iostream>
#include <cmath>

using namespace std;  

DglapEvolution* DglapEvolution::mInstance = 0;  
double DglapEvolution::mAg = 0;  
double DglapEvolution::mLambdaG = 0;  
double DglapEvolution::mMu02 = 0; 
double DglapEvolution::mDerivative = 0;
  
DglapEvolution::DglapEvolution()   
{  
    DipoleModelParameters parameters(TableGeneratorSettings::instance());

    //
    //  For speed purposes the key parameters are held as data member.
    //  Here we assign the proper ones depending on the model and the
    //  parameters set choosen.
    //
    mAg = parameters.Ag();
    mLambdaG = parameters.lambdaG();
    mMu02 = parameters.mu02();
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
  
double DglapEvolution::g0(double x)  
{  
    return mAg*pow(x,-mLambdaG)*pow((1-x),5.6);  //  
}  
  
double DglapEvolution::G(double x, double Q2)  
{  
    //  
    //  To speed up thing during all these integrations  
    //  we keep the last value. Helps a bit.  
    //  
    static bool init = true;  
      
    int N=100;         // maximum Laguerre order  
    double q2i=mMu02;  // initial value of Q^2  
    double q2f=mS;     // max/final value of Q^2   
      
    int Nc=3;  
    tab_Pij *tt_pij;  
    tab_dist *dist0;  
    evol_pair pair;  
    tab_evol *tt_evol;  
    static grid_dist *dist;      
      
    if (init) {  
        cout << "DglapEvolution::G(): initializing DGLAP evolution engine" << endl;  
        
        tt_pij=create_Lag_Pij_table(Nc,N);

        dist0=Lag_dist(N,qns0,qs0,g0);
        dist0->q2=q2i;   
        
        set_LO();
        set_NO_QUARKS();  
        compute_GLUON_DERIVATIVE(); // computes  val.dg = d(xG)/dx  
          
        pair=DGLAP_evol(dist0,q2f,4000,tt_pij);
        tt_evol=pair.pdf;  
        dist=PDF_grid_calc(2000,1.0e-8,tt_evol);  
 
        init = false;
    }  
      
    value_dist val = PDF_interpolated(x, Q2, dist);  
    double result = val.g;  // x*G  
    result /= x;            // G  
      
    mDerivative = val.dg;   // keep
    
    return result;   
}  

double DglapEvolution::dxGdxOfLastCall() const 
{
    return mDerivative;
}
