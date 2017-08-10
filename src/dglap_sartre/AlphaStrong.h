//==============================================================================
//  AlphaStrong.h
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
//  The code was modified for Sartre and made multithread safe and several          
//  methods that are not needed in Sartre were removed.         
//  The formula this code is based on is from:         
//  Particle Data Group, W.-M. Yao et al., J. Phys. G33 (2006) 1         
//         
//  Notes on alpha_s values at M_Z:         
//  In the MRST alpha_s code the following values are used:         
//  LO: 0.13939         
//  NLO: 0.12018         
//  NNLO: 0.11707         
//         
//  PDG (2011) quotes: 0.1184         
//         
//  Pythia uses values varying from 0.12-0.13939 depending on where it is used.         
//           
//  This version with order=1 is identical to MRST LO,         
//  order=2 is close to MRST NLO for Q>2 GeV. order=0 is constant.         
//==============================================================================       
#ifndef AlphaStrong_h         
#define AlphaStrong_h         
         
class AlphaStrong {         
public:         
  //  AlphaStrong(double alphasAtMZ = 0.13, int order = 1);
  AlphaStrong();//double, int); 
             
    void init(double alphasAtMZ, int order);         
    double alphaS(double Q2);  // Q2 is the Q^2 scale in GeV^2         
            
private:         
    double pow2(double) const;         
             
private:             
    static const int    mMaxIter = 10;         
    static const double mMassC;         
    static const double mMassB;         
    static const double mMassZ;          
    static const double mSafetyMargin1;         
    static const double mSafetyMargin2;         
             
private:             
    bool   mIsInit;         
    int    mOrder;         
    double mValueRef;         
    double mValueNow, mScale2Now, mScale2Min;         
    double mMassC2, mMassB2;         
    double mLambda3Save, mLambda4Save, mLambda5Save;         
    double mLambda3Save2, mLambda4Save2, mLambda5Save2;         

    double mLambda3, mLambda4;
    double mLambda2;
};         
         
inline double AlphaStrong::pow2(double x) const {return x*x;}         
         
#endif         
