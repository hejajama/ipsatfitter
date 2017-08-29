//==============================================================================
//  AlphaStrong.cpp
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
//  Author: Thomas Ullrich
//  Last update: 
//  $Date: 2013-05-30 01:56:19 +0530 (Thu, 30 May 2013) $
//  $Author: thomas.ullrich@bnl.gov $
//==============================================================================
//       
//  Original code from Pythia8:     
//  Copyright (C) 2010 Torbjorn Sjostrand.     
//  PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.     
//     
//  The code was modified for Sartre and made multithread safe. Several      
//  methods that are not needed in Sartre where removed.     
//==============================================================================    
#include "AlphaStrong.h"     
//#include "Constants.h"
#include "Parameters_bSat.h"
#include <iostream>     
#include <cmath>
#include <cstdlib>
#define ALPHASTRONG_REMEMBERS_LAST_CALL     

#define PR(x) cout << #x << " = " << (x) << endl;

using namespace std;     
   
//const double AlphaStrong::mMassC   = 1.5;          
const double AlphaStrong::mMassC   = 1.27;
const double AlphaStrong::mMassB   = 4.75;
const double AlphaStrong::mMassZ   = 91.188;           
const double AlphaStrong::mSafetyMargin1 = 1.07;          
const double AlphaStrong::mSafetyMargin2 = 1.33;          

//AlphaStrong::AlphaStrong(double as, int order)
AlphaStrong::AlphaStrong()      
{     
  mIsInit = false;     
  //    init(as, order) ;     
}     
   
void AlphaStrong::init(double value, int order){
  cout<<"Initialising the strong coupling with alpha_s(Mz2)="<<value<<", at order "<<order<<"."<<endl;
  mValueRef  = value;     
  mOrder     = order;
  if(order<0 || order>3){
    cout<<"AlphaStrong::init: chosen order does not exits: "<<order<<", stopping..."<<endl;
    exit(1);
  }
  else if(order==0) { mIsInit=true; }
  else if(order >= 1 && order<= 3){
    mIsInit=true;
    //
    // Need to find the necessary lambda
    // such that as(Mz2)=mValueRef
    // Use bisection method.
    //
    double scale2=mMassZ*mMassZ;
    double delta=90;
    double LambdaA2=1e-5;
    double LambdaB2=1.;
    double LambdaC2=0.5*(LambdaA2+LambdaB2);
    int count=0;
    while(delta>1e-7){
      mLambda2=LambdaA2;
      double asA=alphaS(scale2);
      mLambda2=LambdaC2;
      double asC=alphaS(scale2);
      if((asA-mValueRef)*(asC-mValueRef)<0.)LambdaB2=LambdaC2;
      else LambdaA2=LambdaC2;
      delta=fabs(asC-mValueRef);
      LambdaC2=0.5*(LambdaA2+LambdaB2);
      count++;
      if(count>1e8){
	cout<<"AlphaStrong::init: failed to find LambdaQCD, stopping..."<<endl;
	exit(1);
      }
    }
    //    mLambda2=LambdaC2;
    mLambda2 = mMassZ*exp(-6.*M_PI/(23.*mValueRef));     

    double lam5=sqrt(mLambda2);
    double lam4=sqrt(mMassZ*exp(-6.*M_PI/(25.*mValueRef)));
    double lam3=sqrt(mMassZ*exp(-6.*M_PI/(27.*mValueRef)));
    PR(lam3);
    PR(lam4);
    PR(lam5);
    cout<<"AlphaStrong initialised, AlphaStrong("<<mMassZ*mMassZ<<")="<<alphaS(mMassZ*mMassZ)<<", LambdaQCD="<<sqrt(mLambda2)<<endl;
  }
}

double AlphaStrong::alphaS(double scale2){
  double LambdaQCD2=mLambda2;

  int nf=3;
  if(scale2>mMassC*mMassC) nf=4;
  if(scale2>mMassB*mMassB) nf=5;
  
  //Textbook:
  return mValueRef/(1.+mValueRef/(12.*M_PI) * (33.-2.*nf) * log(scale2/(mMassZ*mMassZ))); 
}


