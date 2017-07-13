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
//#include "TableGeneratorSettings.h"
//#include "DipoleModelParameters.h"
#include <iostream>
#include <cstdlib>
#include <cmath>

using namespace std;  

DglapEvolution* DglapEvolution::mInstance = 0;  

DglapIC dglap_initial_condition;


double DglapEvolution::mDerivative = 0;

// An ugly hack: this is actually modified in Init() const method, but ensuring thread safety
// manually.
//grid_dist *dist;
  
DglapEvolution::DglapEvolution()
{  
    //DipoleModelParameters parameters(TableGeneratorSettings::instance());

    //
    //  For speed purposes the key parameters are held as data member.
    //  Here we assign the proper ones depending on the model and the
    //  parameters set choosen.
    //
    /*
    mAg = parameters.Ag();
    mLambdaG = parameters.lambdaG();
    mMu02 = parameters.mu02();
     */
    //dist = NULL;
    dglap_initial_condition.mAg = 0;
    dglap_initial_condition.mLambdaG = 0;
    dglap_initial_condition.mMu02 = 0;
    
    mS = 1e6;
    
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
  
double qns0(double)
{  
    return 0.0; // no light quarks  
}  
  
double qs0(double)  
{  
    return 0.0; // no strange quarks  
}  
  
double g0(double x)
{  
    return dglap_initial_condition.mAg*pow(x,-dglap_initial_condition.mLambdaG)*pow((1-x),5.6);  //
}


  
grid_dist *dist;
bool distgrid_initialized = false;
double DglapEvolution::G(double x, double Q2)  const
{
    if (Q2 > mS or Q2 < dglap_initial_condition.mMu02)
    {
	    return 0;
	}
    
    value_dist val = PDF_interpolated(x, Q2, dist);  
    double result = val.g;  // x*G  

	if (val.g==0)
			cerr << "xg=0 at x=" << x << "Q2=" << Q2 << endl;


    result /= x;            // G  
      
    // mDerivative = val.dg;   // keep  - can't keep, as this should be thread safe!
    
    return result;   
    
}


// H.M. 201707: move initialization to a separated method, so
// DglapEvolution::G becomes thread safe (hopefully), and can be made const method
int DglapEvolution::Init(double Ag, double lambdag, double mu02) const
{
	if (distgrid_initialized == true)
	{
		delete dist;
	}
    int N=40;         // maximum Laguerre order
    double q2i=mu02;  // initial value of Q^2
    double q2f=mS;     // max/final value of Q^2
    
    dglap_initial_condition.mAg = Ag;
    dglap_initial_condition.mLambdaG = lambdag;
    dglap_initial_condition.mMu02 = mu02;
    
    //static grid_dist *dist;
    
    int Nc=3;
    tab_Pij *tt_pij;
    tab_dist *dist0;
    evol_pair pair;
    tab_evol *tt_evol;
    

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
    
    delete tt_pij;
    delete dist0;
    delete tt_evol;
   
    return 0;

}

double DglapEvolution::dxGdxOfLastCall() const 
{
    cerr << "DglapEvolution::dxGdxOfLastCall() is not working, thanks to thread-safe-changes by Heikki!" << endl;
    return 0;
    //return mDerivative;
}
