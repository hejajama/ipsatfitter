//==============================================================================
//  Laguerre.h
//
//  Copyright (C) 2016 Francois Gelis and Thomas Ullrich
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
//  Author:  F. Gelis/T. Ullrich
//  $Date$
//  $Author$
//==============================================================================
//
//  From the original C version:
//  This file is a translation into C of the original FORTRAN code by
//  Laurent Schoeffel for the LO and NLO splitting functions. The only
//  difference is that the numerical integration is handled by routines
//  of the GNU Scientific Library.
//
//  This version was rewritten in C++(11) by T. Ullrich. The names staid
//  the same to maintain some backwards compatibility. All allocation is
//  now done in the table classes via new and delete.
//==============================================================================
#ifndef Laguerre_h
#define Laguerre_h

class tab_Pij {
public:
    tab_Pij();
    tab_Pij(int sz, int n, int nm);
    ~tab_Pij();
    
public:
    double nf;   // number of flavors
    double nmax; // maximum Laguerre order
    // the following are all tables of type double[nmax+1]
    // with indices going from 0 to nmax. Only the pointers
    // to the first cell of the table are stored in the
    // structure. All routines that return this type must
    // take care to allocate memory for the tables.
    double *pns0;      // LO non singlet
    double *pns1;      // NLO non singlet
    double *psqq0;     // LO singlet:  Pqq
    double *psqq1;     // NLO singlet: Pqq
    double *psqg0;     // LO singlet:  Pqg
    double *psqg1;     // NLO singlet: Pqg
    double *psgq0;     // LO singlet:  Pgq
    double *psgq1;     // NLO singlet: Pgq
    double *psgg0;     // LO singlet:  Pgg
    double *psgg1;     // NLO singlet: Pgg
};

// The same for the Wilson coefficients

class tab_b2 {
public:
    tab_b2();
    tab_b2(int sz, int n, int nm);
    ~tab_b2();
    
public:
    double nf;
    double nmax;
    double *b2q;
    double *b2g;
};

//
// An internal structure used to pass extra parameters to the
// integration routines
//
class data_pair {
public:
    double nf;         // number of flavors
    int nlag;          // Laguerre order
};

double Lag(int n, double x);

tab_Pij *create_Lag_Pij_table(double Nf, int N);
tab_Pij *create_Lag_Pij_table_DLA(double Nf, int N);
tab_b2 *create_Lag_b2_table(double Nf, int N);

#endif
