//==============================================================================
//  Parameters_bSat.h
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
#ifndef parameters_bSat_h  
#define parameters_bSat_h  
  
//Table 3 in KMW  
/*
const double bSat_BG      = 4.;  
const double bSat_Mu02    = 1.17;  
const double bSat_lambdaG = 0.02;  
const double bSat_Ag      = 2.55;  
*/

/*
//KT:
const double bSat_BG      = 4.25;  
const double bSat_Mu02    = 0.77;  
const double bSat_lambdaG = -0.12;  
const double bSat_Ag      = 3.47;  
const double quarkMass[6] = {0.05, 0.05, 0.05, 1.4, 4.67, 172.}; //d, u, s, c from KT; b, t from DPG
const double as_qqg_fixed =0.14;
const double bSat_C       =4.0; //default=4

const double bSat_coupling=0.1184; 
const int bSat_asOrder=1;

const double WS_R0=0;
const double WS_d=0;
//const double bSat_coupling=0.1222; 
*/
/*
//Using HERA PDF NNLO:
const double quarkMass[6] = {0.0, 0.0, 0.0, 1.43, 4.5, 173.}; 
const double bSat_BG      = 4.;

const double bSat_C=16.; //default=4
const double     C2=0.; //Multiplies Q2
const double     C3=4;   //Constant term



const double bSat_Mu02    = 8.5;  

const double bSat_coupling=0.1184; //World average
const int bSat_asOrder=3;


//////////////////////////////
const double bSat_lambdaG = 0.058;  
const double bSat_Ag      = 2.298;  
*/


//RSKV 1212.297:
//const double bSat_BG      = 3.75;  
//const double bSat_BG      = 2.9;  
//const double bSat_BG      = 25.;  
//const double bSat_BG      = 19.55;
enum  ProtonShape {Gauss, HardSphere, WoodsSaxon, Laplace};

const ProtonShape bSat_Shape=Gauss;

const double quarkMass[6] = {0.0, 0.0, 0.0, 1.44, 4.18, 172.}; // b and t not in fit

const double bSat_rmax    = 1.37;   //fm confinement scale
const double bSat_Mu02    = 2.07; //Default 1.51
const double bSat_Ag      = 1.86;  //default: 2.308
const double bSat_lambdaG = 0.119; //Default 0.058  
const double bSat_BG      = 13.; //GeV^-2


const double bSat_C       =4.; //default=4
const double bSat_coupling=0.1184; //World average
const int bSat_asOrder=1;

const double bSat_Ctilde  =48.0; 
//const double as_qqg_fixed =0.101;
const double as_qqg_fixed =0.071;


const double bSat_Rhs =10.1; //GeV-1

const double bSat_R0WS=1.0; //GeV^-1
const double bSat_dWS=1.0; //GeV^-1

const bool correctForRealAmplitude=false;

#endif  
