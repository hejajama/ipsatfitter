//==============================================================================
//  DglapEvolution.h
//
//  Copyright (C) 2010-2013 Tobias Toll and Thomas Ullrich 
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
//  $Date: 2016-11-01 13:24:37 -0400 (Tue, 01 Nov 2016) $
//  $Author: ullrich $
//==============================================================================
#ifndef DglapEvolution_h
#define DglapEvolution_h

class DglapEvolution {
public:
    static DglapEvolution& instance();
    ~DglapEvolution();
    
    double  G(double x, double Q2);
    static  double qns0(double);
    static  double qs0(double);
    static  double g0(double);
    void    setS(double);

    double  dxGdxOfLastCall() const;  

private:
    DglapEvolution();
    
private:
    double mS;
    static double mAg;
    static double mLambdaG;
    static double mMu02;
    static double mDerivative;
    static DglapEvolution* mInstance;
};
 
#endif
