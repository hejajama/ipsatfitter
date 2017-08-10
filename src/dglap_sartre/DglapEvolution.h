//==============================================================================
//  DglapEvolution.h
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
#ifndef DglapEvolution_h
#define DglapEvolution_h
extern "C" {  
#include "laguerre.h"  
#include "dglap.h"  
}  
class Parameters{
 public:
  double       mAg;
  double       mLambdaG;
  double       mMu02;
  double       mBG;
  double       mRmax;

  void setParameters(double, double, double, double, double); //Ag, lambda_g, mu02, BG, rmax

};

class DglapEvolution {
public:
    static DglapEvolution& instance();
    ~DglapEvolution();
    
    double  G(double x, double Q2);
    static  double qns0(double);
    static  double qs0(double);
    static  double g03(double);
    static  double g04(double);
    static  double g05(double);
    void    setS(double);

    double  dxGdxOfLastCall() const;  

    Parameters params;
    
//private:
    DglapEvolution();
    
private:
    double mS;
    static double mAg;
    static double mLambdaG;
    static double mMu02;
    static double mDerivative;
    static DglapEvolution* mInstance;
    static double mCharmThreshold;
    static double mBeautyThreshold;

    static grid_dist* mDist3;      
    static grid_dist* mDist4;
    static grid_dist* mDist5;      

};
 
#endif
