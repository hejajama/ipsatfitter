//==============================================================================
//  Dglap.h
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
//  This code is based on the paper "An elegant and fast method to
//  solve QCD evolution equations" by Laurent Schoeffel. The resolution
//  of the ODE relies on routines of the GNU Scientific Library.
//
//  This version was rewritten in C++(11) by T. Ullrich. The names staid
//  the same to maintain some backwards compatibility. All allocation is
//  now done in the table classes via new and delete.
//==============================================================================
#ifndef Dglap_h
#define Dglap_h
#include "Laguerre.h"

//
// Internal structure used to pass arguments to the routine that
// calculates Laguerre coefficients
//
class lag_pair {
public:
    int nlag;                      // Laguerre order
    double (*function)(double);    // a function of a real variable
};

//
// A structure that holds the Laguerre transform of some PDFs
//
class tab_dist {
public:
    tab_dist(int sz, int nm, double qsq, bool om = true);
    ~tab_dist();
    
public:
    int nmax;            // maximum Laguerre order
    double q2;           // the value of Q^2
    double *qns;         // quark non singlet xQns(x)
    double *qs;          // quark singlet     xQs(x)
    double *g;           // gluon             xG(x)
    
    bool mOwnMemory;
};

//
// A structure that holds the values of the PDFs at some x as a
// function of Q^2
//
class tab_dist_x {
public:
    tab_dist_x(int sz, int nqsq, double xx, double nff, double *qsq);
    ~tab_dist_x();

public:
    int nq2;             // number of values of Q^2
    double x;            // value of x
    double nf;           // number of flavors in the evolution
    double *q2;          // table of values of Q^2  q2[0..nq2]
    double *qns;         // corresponding values of xQns(x)
    double *qs;          // ....................... xQs(x)
    double *g;           // ....................... xG(x)
    double *dg;          // ....................... d[xG(x)]/dx
};

//
// A structure that holds the values of the PDFs on a grid in x and
// Q^2.
//
class grid_dist {
public:
    grid_dist(int nsq, int nxx, double nff, double xm);
    ~grid_dist();

public:
    int nq2;             // number of values of Q^2
    int nx;              // number of values of x
    double nf;           // number of flavors in the evolution
    double xmin;         // minimum value of x (the max value is x=1.0)
    double *q2;          // table of values of Q^2 q2[0..nq2]
    double *x;           // table of values of x    x[0..nx]
    // Then, we have the values of the PDFs at the corresponding points.
    // They are 2d tables table[0..nx][0..nq2] such that table[i][j]
    // is the value of the PDF at x=x[i] and Q^2=q2[j]. BEWARE that since
    // they are declared as pointers to a double and not as a 2d table, the
    // compiler does not know the length of a row in order to make any sense
    // of the notation table[i][j]. Instead, one should use table[i*(1+nq2)+j]
    double *qns;         // values of xQns(x,Q^2)     qns[0..nx][0..nq2]
    double *qs;          // values of xQs(x,Q^2)       qs[0..nx][0..nq2]
    double *g;           // values of xG(x,Q^2)         g[0..nx][0..nq2]
    double *dg;          // values of d[xG(x,Q^2)]/dx  dg[0..nx][0..nq2]
};

//
// A structure that holds the values of the PDFs at a single point (x,Q^2)
//
class value_dist {
public:
    double x;            // value of x
    double q2;           // value of Q^2
    double qns;          // value of xQns(x,Q^2)
    double qs;           // value of xQs(x,Q^2)
    double g;            // value of xG(x,Q^2)
    double dg;           // value of d[xG(x,Q^2)]/dx
};

//
// A structure that holds the solution of the differential equations
// for the Laguerre coefficients of the PDFs
//
class tab_evol {
public:
    tab_evol(int nm, int nqsq, int nff);
    ~tab_evol();
    
public:
    int nmax;           // maximum Laguerre order
    int nq2;            // number of values of Q^2
    double nf;          // value of Nf used in the evolution
    double *q2;         // table of values of Q^2  q2[0..nq2]
    // The corresponding values of the Laguerre coefficients of the PDFs.
    // table[i][j] is the j-th Laguerre coefficient of the PDF at Q^2=q2[i].
    // Because they are declared as pointers to a double, these tables
    // cannot be accessed as table[i][j], but instead as table[i*(1+nmax)+j]
    double *qns;        // table qns[0..nq2][0..nmax]
    double *qs;         // table  qs[0..nq2][0..nmax]
    double *g;          // table   g[0..nq2][0..nmax]
};

//
// A structure that holds the running coupling constant as a function
// of Q^2
//
class tab_alpha {
public:
    tab_alpha(int sz, int nqsq);
    ~tab_alpha();

public:
    int nq2;            // number of values of Q^2
    double *q2;         // table of values of Q^2:       q2[0..nq2]
    double *alpha;      // table of values of alpha:  alpha[0..nq2]
};

//
// Internal data structure for the ODEs
//
class ode_pair {
public:
    tab_Pij *pij;
    tab_alpha *alpha;
};

//
// A structure that holds a pointers to a variable of type tab_evol
// and one to a variable of type tab_alpha
//
class evol_pair {
public:
    tab_evol *pdf;
    tab_alpha *alpha;
};

double Lag_coeff(int n,double (*f)(double));
double Lag_inv(double x,int N,double *coeff);
double diff_Lag_inv(double x,int N,double *coeff);

tab_dist *Lag_dist(int N,double (*qns)(double),double (*qs)(double),double (*g)(double));

value_dist Lag_dist_inv(double x,tab_dist *tt);

double alpha_s_LO(double Q2,double Nf);
tab_alpha *alpha_s_evol(double,double,double,double,int,double);
double alpha_s(double,tab_alpha);

void set_LO();
void set_NO_QUARKS();
void set_NO_RUNNING();
void compute_GLUON_DERIVATIVE();


evol_pair DGLAP_evol(tab_dist *dist0,double Q2max,int NQ2,tab_Pij *pij);

tab_dist *get_step(int j,tab_evol *tt_evol);
tab_dist *PDF_at_q2(double q2,tab_evol *tt_evol);

tab_dist_x *PDF_at_x(double x,tab_evol *tt_evol);
grid_dist *PDF_grid_calc(int Nx, double xmin,tab_evol *tt_evol);

void PDF_grid_write(grid_dist *tt,char *name);
grid_dist *PDF_grid_read(char *name);

value_dist PDF_interpolated(double x,double q2,grid_dist *tt);

#endif
