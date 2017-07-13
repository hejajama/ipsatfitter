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

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <string.h>
#include <string>
#include <signal.h>
#include <unistd.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_sf_dilog.h>
#include <gsl/gsl_integration.h>
#include "Dglap.h"

using namespace std;

#define nullptr 0

//
//  Constants:
//  Nc=3 is fixed once for all, Nf remains free
//
#define Nc 3.0
#define Cf 1.333333333333333333333
#define Cg 3.0
#define Tr 0.5

//
// this counter is used by the function sigalarm_handler() to display
// something about the progress of the calculation
//
static int counter;

//
// this global variable is a flag to tell other routines whether we do
// LO or NLO calculations. It should be set by the function set_LO()
//
static double NLO=1.0;  // by default, we do NLO calculations 

//
// this global variable is used to force a "pure glue" evolution of
// xG(x,Q^2) in order to test the hypothetical situation where the
// quarks would be infinitely massive
//
static double QUARKS=1.0; // by default, we keep the quarks

//
// this global variable is used to force a constant alpha_s
//
static double RUNNING=1.0; // by default, alpha_s is running

//
// global variable to compute the x-derivative of the gluon distribution
//
static int GLUON_DERIVATIVE=0; // by default, we do not compute the derivative

//
//  Constructors and destructors
//
tab_dist::tab_dist(int sz, int nm, double qsq, bool om)
{
    if (sz) {
        qns = new double[sz];
        qs = new double[sz];
        g = new double[sz];
    }
    else {
        qns = nullptr;
        qs = nullptr;
        g = nullptr;
    }
    nmax = nm;
    q2 = qsq;
    
    mOwnMemory = om;   // does memory belong to me or not
}

tab_dist::~tab_dist()
{
    if (mOwnMemory) {
        delete [] qns;
        delete [] qs;
        delete [] g;
    }
}

tab_alpha::tab_alpha(int sz, int nqsq)
{
    q2 = new double[sz];
    alpha =  new double[sz];
    nq2 = nqsq;
}

tab_alpha::~tab_alpha()
{
    delete [] q2;
    delete [] alpha;
}

tab_evol::tab_evol(int nm, int nqsq, int nff)
{
    nmax=nm;
    nq2=nqsq;
    nf=nff;
    q2 = new double[nq2+1];
    qns = new double[(nq2+1)*(nmax+1)];
    qs = new double[(nq2+1)*(nmax+1)];
    g = new double[(nq2+1)*(nmax+1)];
}

tab_evol::~tab_evol()
{
    delete [] q2;
    delete [] qns;
    delete [] qs;
    delete [] g;
}

tab_dist_x::tab_dist_x(int sz, int nqsq, double xx, double nff, double *qsq)
{
    nq2=nqsq;
    x=xx;
    q2=qsq;
    nf=nff;
    qns = new double[sz];
    qs = new double[sz];
    g  = new double[sz];
    dg = new double[sz];
}

tab_dist_x::~tab_dist_x()
{
    delete [] qns;
    delete [] qs;
    delete [] g;
    delete [] dg;
}

grid_dist::grid_dist(int nsq, int nxx, double nff, double xm)
{
    nq2 = nsq;
    nx = nxx;
    nf = nff;
    xmin = xm;
    
    q2 = new double[1+nq2];
    x = new double[1+nx];
    qns = new double[(1+nq2)*(1+nx)];
    qs = new double[(1+nq2)*(1+nx)];
    g  = new double[(1+nq2)*(1+nx)];
    dg = new double[(1+nq2)*(1+nx)];
}

grid_dist::~grid_dist()
{
    delete [] q2;
    delete [] x;
    delete [] qns;
    delete [] qs;
    delete [] g;
    delete [] dg;
}

//
// A function which is called every 10 seconds and prints something
// about the state of the calculation, for the routines that are
// particularly slow
//
static void sigalarm_handler(int sig)
{
    (void) sig;  // silence clang warnings
    cout << counter << "% done" << endl;
    alarm(10); /* let's request another alarm 5 seconds from now */
#ifdef not_linux // On SOLARIS at least, one needs to reinstall the
    // signal handler after every use. Just give the extra option
    // "-Dnot_linux" to the compiler when compiling this file
    signal( SIGALRM, sigalarm_handler);
#endif
}

//
// Set up for LO calculations
//
void set_LO()
{
    NLO=0.0;
}

//
// Set things up so that quarks are neglected in the evolution of gluons
//
void set_NO_QUARKS()
{
    QUARKS=0.0;
}

//
// Set things up so that alpha_s is taken to be constant
//
void set_NO_RUNNING()
{
    RUNNING=0.0;
}

//
// Set things up to compute the x-derivative of xG(x,Q^2)
//
void compute_GLUON_DERIVATIVE()
{
    GLUON_DERIVATIVE=1;
}

//
// An intermediate function used in the computation of the Laguerre
// coefficients of a given function. It is not meant to be called
// directly.
//
double Lag_integrand(double z,void *params)
{
    double x=exp(-z);
    lag_pair *data = reinterpret_cast<lag_pair*>(params);
    int n=data->nlag;
    double (*f)(double)=data->function;
    
    return x*f(x)*Lag(n,z);
}

//
// Calculate the n-th Laguerre coefficient for the function f(x)
//
double Lag_coeff(int n,double (*f)(double))
{
    double epsilon=1.0e-12;
    gsl_integration_workspace *wksp=gsl_integration_workspace_alloc(200);
    double result,abserr;
    lag_pair params;
    gsl_function F;
    
    params.nlag=n;
    params.function=f;
    
    F.function=Lag_integrand;
    F.params=&params;
    
    gsl_integration_qag(&F,0.0,100.0,epsilon,epsilon,200,6,wksp,&result,&abserr);
    gsl_integration_workspace_free(wksp);
    
    return result;
}

//
// Reconstruct a function at point x from its Laguerre coefficients
//
double Lag_inv(double x, int N, double *coeff)
{
    // x: point at which we want the value of the function
    // N: maximum Laguerre order
    // coeff: the table[0..N] of the Laguerre coefficients
    double t=-log(x);
    double *lag0 = new double[N+2];
    double sum=0.0;
    
    lag0[0]=1.0;
    lag0[1]=1.0-t;
    for (int i=2;i<=N;i++){
        double a1=(static_cast<double>(2*i-1)-t)/i;
        double a2=(static_cast<double>(i-1))/i;
        lag0[i]=a1*lag0[i-1]-a2*lag0[i-2];
    }
    
    for (int i=N;i>=0;i--) sum+=coeff[i]*lag0[i];
    
    delete [] lag0;
    
    return sum;
}

//
// Reconstruct the derivative of a function at point x from its Laguerre coefficients
//
double diff_Lag_inv(double x,int N,double *coeff)
{
    // x: point at which we want the derivative of the function
    // N: maximum Laguerre order
    // coeff: the table[0..N] of the Laguerre coefficients
    double t=-log(x);
    double *lag0 = new double[N+2];
    double *dlag0 = new double[N+2];
    double sum=0.0;
    
    lag0[0]  = 1.0;
    dlag0[0] = 0.0;
    lag0[1]  = 1.0-t;
    dlag0[1] = -1.0;
    for (int i=2;i<=N;i++){
        double a1  = (static_cast<double>(2*i-1)-t)/i;
        double da1 = -1.0/((double)(i));
        double a2  = static_cast<double>(i-1)/i;
        lag0[i]  = a1*lag0[i-1] - a2*lag0[i-2];
        dlag0[i] = da1*lag0[i-1] + a1*dlag0[i-1] - a2*dlag0[i-2];
    }
    
    for (int i=N;i>=0;i--) sum+=coeff[i]*dlag0[i];
    
    delete [] lag0;
    delete [] dlag0;
    
    return -sum/x;
}

//
// Calculate the Laguerre coefficients of some PDFs passed as arguments
//
tab_dist *Lag_dist(int N,double (*qns)(double),double (*qs)(double),double (*g)(double))
{
    // N: maximum Laguerre order
    // qns: pointer to the non singlet quark PDF
    // qs:  pointer to the singlet quark PDF
    // g:   pointer to the gluon PDF
    tab_dist *tt = new tab_dist(1+N, N, 1.0);
    
    cout << "Calculating Laguerre coefficients of the initial PDFs..." << flush;
    
    for (int i=0;i<=N;i++){
        (tt->qns)[i]=Lag_coeff(i,qns);
        (tt->qs)[i]=Lag_coeff(i,qs);
        (tt->g)[i]=Lag_coeff(i,g);
    }
    
    cout << "done" << endl;
    
    return tt;
}

//
// Calculate the PDFs at some x from a table of Laguerre coefficients
//
value_dist Lag_dist_inv(double x, tab_dist *tt)
{
    value_dist out;
    
    out.x=x;
    out.q2=tt->q2;
    out.qns=Lag_inv(x,tt->nmax,tt->qns);
    out.qs=Lag_inv(x,tt->nmax,tt->qs);
    out.g=Lag_inv(x,tt->nmax,tt->g);
    if (GLUON_DERIVATIVE) {
        out.dg=diff_Lag_inv(x,tt->nmax,tt->g);
    }
    else {
        out.dg=0.0;
    }
    
    return out;
}

//
// The jacobian function, empty here since we don't need it (we use
// only explicit ODE solvers).
//
int jac(double t,const double y[],double *dfdy,double dfdt[],void *params)
{
    // silence clang warnings
    (void) t;
    (void) y;
    (void) dfdy;
    (void) dfdt;
    (void) params;
    
    return GSL_SUCCESS;
}

//
// The r.h.s. for the ODE that controls the evolution of
// alpha_s. NOTE: I do not take into account the heavy quark masses in
// the evolution, nor the fact that the number of active flavors
// change with Q^2.
//
int rhs_alpha(double q2,const double y[],double f[],void *params)
{
    double Nf=*(reinterpret_cast<double*>(params));
    double beta0=11.0-2.0*Nf/3.0;
    double beta1=102.0-38.0*Nf/3.0;
    double four_PI=4.0*M_PI;
    double alpha2=y[0]*y[0]/four_PI;
    double alpha3=alpha2*y[0]/four_PI;
    
    // If the global variable NLO is 0.0, we include only the LO beta
    // function
    
    f[0]=-RUNNING*(beta0*alpha2+NLO*beta1*alpha3)/q2;
    
    return GSL_SUCCESS;
}


//
// The routine that actually solves the ODE in order to get the
// running of alpha_s
//
tab_alpha *alpha_s_evol(double q20,double alpha0,double q2min,double q2max,int NQ2,double Nf)
{
    // q20: the value of Q^2 at which alpha_s is given
    // alpha0: the corresponding value of alpha_s
    // q2min: the min value of Q^2 in the produced table
    // q2max: the maximum value of Q^2..................
    // NQ2: the number of value of Q^2
    // Nf: the number of flavors
    double params[1]={Nf};
    const gsl_odeiv_step_type *T= gsl_odeiv_step_rk8pd;
    gsl_odeiv_step *s=gsl_odeiv_step_alloc(T,1);
    double ode_eps=1.0e-6;
    gsl_odeiv_control *c=gsl_odeiv_control_y_new(ode_eps,ode_eps);
    gsl_odeiv_evolve *e=gsl_odeiv_evolve_alloc(1);
    gsl_odeiv_system sys={rhs_alpha,jac,1,params};
    double h;
    double y[1];
    double q2;
    double dlq2=exp(log(q2max/q2min)/NQ2);
    
    tab_alpha *tt = new tab_alpha(1+NQ2, NQ2);
    
    cout << "Calculating the evolution of alpha_s..." << flush;
    
    // initialize the table of the Q^2's
    
    q2=q2min;
    int j=-1;
    for (int i=0;i<=NQ2;i++) {
        tt->q2[i]=q2;
        if (q2<q20) j=i; // this determines the index in the table just
        // before q2[i]
        q2*=dlq2;
    }
    
    // Forward evolution of alpha_s
    
    q2=q20;
    y[0]=alpha0;
    h=1.0e-10;
    
    for (int i=j+1;i<=NQ2;i++) {//if Q0^2>Qmax^2, then j=NQ2 and this loop is skept
        while (q2<tt->q2[i]) gsl_odeiv_evolve_apply(e,c,s,&sys,&q2,tt->q2[i],&h,y);
        tt->alpha[i]=y[0];
    }
    
    // Backward evolution of alpha_s
    
    q2=q20;
    y[0]=alpha0;
    h=-1.0e-10;
    
    for (int i=j;i>=0;i--) {//if Q0^2<Qmin^2, then j=-1 and this loop is skept
        while (q2>tt->q2[i]) gsl_odeiv_evolve_apply(e,c,s,&sys,&q2,tt->q2[i],&h,y);
        tt->alpha[i]=y[0];
    }
    
    // Free the allocated workspace
    
    gsl_odeiv_evolve_free(e);
    gsl_odeiv_step_free(s);
    gsl_odeiv_control_free(c);
    
    cout << "done" << endl;
    
    return tt;
}

//
// Running alpha_s obtained by interpolation from a table of
// precomputed values
//
double alpha_s(double Q2,tab_alpha *tt)
{
    int NQ2=tt->nq2;
    int i=(int)floor((static_cast<double>(NQ2))*log(Q2/tt->q2[0])/log(tt->q2[NQ2]/tt->q2[0]));
    double f1,f2,x1,x2;
    
    x1=tt->q2[i];
    x2=tt->q2[i+1];
    f1=tt->alpha[i];
    f2=tt->alpha[i+1];
    
    return f1+(f2-f1)*(Q2-x1)/(x2-x1);
}

//
// LO running alpha_s -- direct calculation
//
double alpha_s_LO(double Q2,double Nf)
{
    double Mz=91.18; // Z0 mass in GeV
    double inv_alpha_s_MZ=1.0/0.113;
    double beta0=(11.0-2.0*Nf/3.0)/(4.0*M_PI);
    double dlog=log(Mz*Mz/Q2);
    
    return 1.0/(inv_alpha_s_MZ-beta0*dlog);
}

//
// The r.h.s. of the ODE in the non singlet case
// y[0]..y[N]=qns
//
int rhs_NS(double q2,const double y[],double f[],void *params)
{
    tab_Pij *tt = reinterpret_cast<ode_pair*>(params)->pij;
    tab_alpha *tta = reinterpret_cast<ode_pair*>(params)->alpha;
    double alpha_s_2pi0=alpha_s(q2,tta)/(2.0*M_PI);
    double alpha_s_2pi=NLO*alpha_s_2pi0;
    double sum;
    int N=tt->nmax;
    
    f[0]=alpha_s_2pi0*y[0]*((tt->pns0)[0]+alpha_s_2pi*(tt->pns1)[0])/q2;
    
    for (int i=1;i<=N;i++){
        sum=0.0;
        for (int j=0;j<i;j++) {
            sum+=y[j]*((tt->pns0)[i-j]+alpha_s_2pi*(tt->pns1)[i-j]\
                       -(tt->pns0)[i-j-1]-alpha_s_2pi*(tt->pns1)[i-j-1]);
        }
        sum+=y[i]*((tt->pns0)[0]+alpha_s_2pi*(tt->pns1)[0]);
        f[i]=sum*alpha_s_2pi0/q2;
    }
    
    return GSL_SUCCESS;
}

//
// The r.h.s. of the ODE in the singlet case
// y[0]..y[N]=qs
// y[N+1]..y[2N+1]=g
//
int rhs_S(double q2,const double y[],double f[],void *params)
{
    tab_Pij *tt = reinterpret_cast<ode_pair*>(params)->pij;
    tab_alpha *tta = reinterpret_cast<ode_pair*>(params)->alpha;
    double alpha_s_2pi0=alpha_s(q2,tta)/(2.0*M_PI);
    double alpha_s_2pi=NLO*alpha_s_2pi0;
    double sum_q;
    double sum_g;
    //double Nf=tt->nf;
    int N=tt->nmax;
    
    f[0]=alpha_s_2pi0*(y[0]*((tt->psqq0)[0]+alpha_s_2pi*(tt->psqq1)[0])\
                       +y[N+1]*((tt->psqg0)[0]+alpha_s_2pi*(tt->psqg1)[0]))/q2;
    f[N+1]=alpha_s_2pi0*(QUARKS*y[0]*((tt->psgq0)[0]+alpha_s_2pi*(tt->psgq1)[0])\
                         +y[N+1]*((tt->psgg0)[0]+alpha_s_2pi*(tt->psgg1)[0]))/q2;
    
    for (int i=1;i<=N;i++){
        sum_q=0.0;
        sum_g=0.0;
        for (int j=0;j<i;j++) {
            sum_q+=y[j]*((tt->psqq0)[i-j]+alpha_s_2pi*(tt->psqq1)[i-j]\
                         -(tt->psqq0)[i-j-1]-alpha_s_2pi*(tt->psqq1)[i-j-1])\
            +y[j+N+1]*((tt->psqg0)[i-j]+alpha_s_2pi*(tt->psqg1)[i-j]\
                       -(tt->psqg0)[i-j-1]-alpha_s_2pi*(tt->psqg1)[i-j-1]);
            sum_g+=QUARKS*y[j]*((tt->psgq0)[i-j]+alpha_s_2pi*(tt->psgq1)[i-j]\
                                -(tt->psgq0)[i-j-1]-alpha_s_2pi*(tt->psgq1)[i-j-1])\
            +y[j+N+1]*((tt->psgg0)[i-j]+alpha_s_2pi*(tt->psgg1)[i-j]\
                       -(tt->psgg0)[i-j-1]-alpha_s_2pi*(tt->psgg1)[i-j-1]);
        }
        sum_q+=y[i]*((tt->psqq0)[0]+alpha_s_2pi*(tt->psqq1)[0])\
        +y[i+N+1]*((tt->psqg0)[0]+alpha_s_2pi*(tt->psqg1)[0]);
        sum_g+=QUARKS*y[i]*((tt->psgq0)[0]+alpha_s_2pi*(tt->psgq1)[0])\
        +y[i+N+1]*((tt->psgg0)[0]+alpha_s_2pi*(tt->psgg1)[0]);
        f[i]=sum_q*alpha_s_2pi0/q2;
        f[i+N+1]=sum_g*alpha_s_2pi0/q2;
    }
    
    return GSL_SUCCESS;
}

//
// DGLAP evolution of the PDFs. This routine actually solves the ODEs
// for the Laguerre coefficients of the PDFs. It returns a structure
// of type evol_pair which contains the evolution with Q^2 of the
// Laguerre coefficients of the PDFs, and the evolution with Q^2 of
// alpha_s.
// NOTE: this routine relies heavily on the ODE solvers in the GNU
// scientific library. In principle, there is no problem to rewrite it
// with some "homemade" ODE solver, but that would make the code
// unecessarily messy...
//
evol_pair DGLAP_evol(tab_dist *dist0,double Q2max,int NQ2,tab_Pij *pij)
{
    // dist0: Laguerre coefficients of the initial PDFs. They must be
    // computed with the function Lag_dist. Note: the value of Q0^2 must
    // be defined in the field dist0->q2.
    // Q2max: value of Q^2 where we stop the evolution
    // NQ2: number of steps in the evolution. The intermediate values of
    // Q^2 will be equally spaced on a log scale.
    // pij: the table of Laguerre coefficients of the splitting
    // functions. This is calculated with the function
    // create_Lag_Pij_table in laguerre.c
    int N=dist0->nmax;
    ode_pair params;
    const gsl_odeiv_step_type *T= gsl_odeiv_step_rk8pd;
    gsl_odeiv_step *s_NS=gsl_odeiv_step_alloc(T,N+1);
    gsl_odeiv_step *s_S=gsl_odeiv_step_alloc(T,2*N+2);
    double ode_eps=1.0e-6;
    gsl_odeiv_control *c_NS=gsl_odeiv_control_y_new(ode_eps,ode_eps);
    gsl_odeiv_control *c_S=gsl_odeiv_control_y_new(ode_eps,ode_eps);
    gsl_odeiv_evolve *e_NS=gsl_odeiv_evolve_alloc(N+1);
    gsl_odeiv_evolve *e_S=gsl_odeiv_evolve_alloc(2*N+2);
    gsl_odeiv_system sys_NS={rhs_NS,jac,static_cast<size_t>(N+1),&params};
    gsl_odeiv_system sys_S={rhs_S,jac,static_cast<size_t>(2*N+2),&params};
    
    double *y_NS = new double[N+2];
    double *y_S = new double[2*N+4];
    double h_NS,h_S;
    double q2_NS,q2_S;
    double q2;
    double q2i=dist0->q2;
    double fact_q2=exp(log(Q2max/q2i)/NQ2);
    evol_pair out;
    double Mz=91.18; // The mass od the Z0
    
    params.pij=pij;

    tab_evol *tt= new tab_evol(N, NQ2, pij->nf);

    //
    // first, compute the running alpha_s as a function of Q^2. We set
    // the value of alpha_s in the MSBAR scheme at the Z0 mass.
    //
    params.alpha=alpha_s_evol(Mz*Mz,0.113,q2i,Q2max,NQ2,pij->nf);
    
    cout << "Solving the ODE..." << flush;
    
    //
    // Set the initial conditions
    //
    for (int i=0;i<=N;i++) {
        y_NS[i]=(dist0->qns)[i];
        y_S[i]=(dist0->qs)[i];
        y_S[N+1+i]=(dist0->g)[i];
    }
    
    q2_NS=q2i;
    q2_S=q2i;
    q2=q2i;
    h_NS=1.0e-10;
    h_S=1.0e-10;
    
    //
    // Store the initial condition
    //
    (tt->q2)[0]=q2;
    for (int i=0;i<=N;i++){
        (tt->qns)[i]=y_NS[i];
        (tt->qs)[i]=y_S[i];
        (tt->g)[i]=y_S[N+1+i];
    }
    
    //
    // The main loop: NQ2 steps
    //
    for (int j=1;j<=NQ2;j++){
        q2=q2*fact_q2;
        while (q2_NS<q2) gsl_odeiv_evolve_apply(e_NS,c_NS,s_NS,&sys_NS,&q2_NS,q2,&h_NS,y_NS);
        while (q2_S<q2) gsl_odeiv_evolve_apply(e_S,c_S,s_S,&sys_S,&q2_S,q2,&h_S,y_S);
        
        //
        // Store the result of this step. See the function "get_step"
        // below in order to understand why it is better to use
        // [i+(1+N)*j] instead of [j+(1+NQ2)*i].
        //
        (tt->q2)[j]=q2;
        for (int i=0;i<=N;i++){
            (tt->qns)[i+(1+N)*j]=y_NS[i];
            (tt->qs)[i+(1+N)*j]= y_S[i];
            (tt->g)[i+(1+N)*j]= y_S[N+1+i];
        }
    }
    
    out.alpha=params.alpha;
    out.pdf=tt;
    
    //
    // Free the allocated workspace
    //
    gsl_odeiv_evolve_free(e_NS);
    gsl_odeiv_evolve_free(e_S);
    gsl_odeiv_step_free(s_NS);
    gsl_odeiv_step_free(s_S);
    gsl_odeiv_control_free(c_NS);
    gsl_odeiv_control_free(c_S);
    
    cout << "done" << endl;
    
    delete [] y_NS;
    delete [] y_S;
    
    return out;
}

//
// extract step #j of the Q^2 evolution - useful if one needs one
// value of Q^2 that sits on the grid and many values of x
//
tab_dist *get_step(int j, tab_evol *tt_evol)
{
    //
    // j: the # of the step
    // tt_evol: the structure returned by DGLAP_evol
    //
    int N=tt_evol->nmax;
    
    //
    // We do not allocate any memory for the tables in tt. Instead, we
    // make them point to the correct location in the structure
    // tt_evol. This saves some memory and makes things faster, but
    // assumes that the variable tt_evol is not altered in any way
    // before one is done using the variable tt...
    //
    // TU: Outch, this is the place where things can go bitterly wrong
    //     hence for the C++ version we build in some protection.
    // tab_dist class last argument decides on the ownership (here false)
    // which prevents ~tab_dist to delete the memory.
    //

    tab_dist *tt = new tab_dist(0, N, (tt_evol->q2)[j], false);
    
    tt->qns=&((tt_evol->qns)[j*(1+N)]); // the data in the table tt_evol->qns
    tt->qs = &((tt_evol->qs)[j*(1+N)]); // is laid out in such a way that the
    tt->g  =  &((tt_evol->g)[j*(1+N)]); // Laguerre coefficients for the j-th
    // value of Q^2 form a sub-vector of
    // length 1+N that starts at the index
    // j*(1+N). Hence that trick to get
    // a pointer to this sub-vector.
    
    return tt;
}


//
// Get the Laguerre coefficients of the PDFs at some intermediate Q^2
// by interpolation -- useful if all one needs is a specific value of
// Q^2 and many values of x
//
tab_dist *PDF_at_q2(double q2,tab_evol *tt_evol)
{
    // q2: the value of Q^2
    // tt_evol: the structure returned by DGLAP_evol
    int N=tt_evol->nmax;
    int NQ2=tt_evol->nq2;
    double dq2;

    tab_dist *tt = new tab_dist(1+N, N, q2);
    
    if ((q2<=(tt_evol->q2)[0])||(q2>=(tt_evol->q2)[NQ2])) {
        cout << "ERROR: q2=" << q2 << " out of range\n" << endl;
        delete tt;
        return nullptr;
    }
    
    int j=floor(NQ2*log(q2/(tt_evol->q2)[0])/log((tt_evol->q2)[NQ2]/(tt_evol->q2)[0]));
    
    dq2=(q2-(tt_evol->q2)[j])/((tt_evol->q2)[j+1]-(tt_evol->q2)[j]);
    
    for (int i=0;i<=N;i++){
        (tt->qns)[i]=(1-dq2)*(tt_evol->qns)[j*(1+N)+i]\
        +dq2*(tt_evol->qns)[(j+1)*(1+N)+i];
        (tt->qs)[i]=(1-dq2)*(tt_evol->qs)[j*(1+N)+i]\
        +dq2*(tt_evol->qs)[(j+1)*(1+N)+i];
        (tt->g)[i]=(1-dq2)*(tt_evol->g)[j*(1+N)+i]\
        +dq2*(tt_evol->g)[(j+1)*(1+N)+i];
    }
    
    return tt;
}

//
// Calculate PDFs at some value of x as a function of Q^2
//
tab_dist_x *PDF_at_x(double x,tab_evol *tt_evol)
{
    int NQ2=tt_evol->nq2;
    int N=tt_evol->nmax;
   
    // x: the value of x
    // tt_evol: the structure returned by DGLAP_evol
    tab_dist_x *tt = new tab_dist_x(1+NQ2, NQ2, x, tt_evol->nf, tt_evol->q2);
    
    for (int j=0;j<=NQ2;j++){
        (tt->qns)[j]=Lag_inv(x,N,&((tt_evol->qns)[(1+N)*j]));
        (tt->qs)[j]= Lag_inv(x,N,&((tt_evol->qs)[(1+N)*j]));
        (tt->g)[j]=  Lag_inv(x,N,&((tt_evol->g)[(1+N)*j]));
        if (GLUON_DERIVATIVE) (tt->dg)[j]=  diff_Lag_inv(x,N,&((tt_evol->g)[(1+N)*j]));
    }
    
    return tt;
}

//
// Calculate the values of the PDFs on a rectangular grid in x and
// Q^2, from the table of Laguerre coefficients. In principle, one can
// compute values at any point x directly from the table of Laguerre
// coefficients. However, this is a rather slow process since one has
// to sum over the Laguerre polynomials. Therefore, it is better to
// precompute a table of values directly as a function of x. Then,
// this table can be stored in a file on disk, and re-read later.
// NOTE: the points x and Q^2 are linearly spaced on a log scale.
//
grid_dist *PDF_grid_calc(int Nx, double xmin, tab_evol *tt_evol)
{
    //
    // Nx: the number of points in the x direction
    // xmin: the minimum value of x (xmax is always 1.0)
    // tt_evol: the structure returned by tt_evol
    //
    int NQ2=tt_evol->nq2;
    int N=tt_evol->nmax;
    double dlx=exp(log(xmin)/((double)Nx));
    double x;
 
    //
    // Install a signal handler
    //
    signal( SIGALRM, sigalarm_handler );
    alarm(10);

    grid_dist *tt = new grid_dist(NQ2, Nx, tt_evol->nf, xmin);
    
    cout << "Calculating a grid of values..." << endl;
    
    //
    // Build the table
    //
    for (int i=0; i<1+NQ2; i++) tt->q2[i] = tt_evol->q2[i];
    
    x=1.0;
    for (int i=Nx;i>=0;i--) {
        (tt->x)[i]=x;
        for (int j=0;j<=NQ2;j++){
            counter=static_cast<int>(100.0*((Nx-i)*(1+NQ2)+j)/(static_cast<double>((1+NQ2)*(1+Nx))));
            (tt->qns)[i*(1+NQ2)+j]= Lag_inv(x,N,&((tt_evol->qns)[(1+N)*j]));
            (tt->qs)[i*(1+NQ2)+j]= Lag_inv(x,N,&((tt_evol->qs)[(1+N)*j]));
            (tt->g)[i*(1+NQ2)+j]= Lag_inv(x,N,&((tt_evol->g)[(1+N)*j]));
            if (GLUON_DERIVATIVE) (tt->dg)[i*(1+NQ2)+j]= diff_Lag_inv(x,N,&((tt_evol->g)[(1+N)*j]));
        }
        x*=dlx;
    }
    
    //
    // remove the signal handler
    //
    signal( SIGALRM, SIG_IGN);
    
    return tt;
}


//
// Write the grid on disk. NOTE: the file generated by this routine is
// a binary file in order to simplify the parsing when it is read by
// PDF_grid_read. This means that it can be read only on computers
// where the binary representation of floating point numbers uses the
// same byte ordering... This file can then be read with the function
// PDF_grid_read, which returns a pointer to a structure of type
// "grid_dist".
//
void PDF_grid_write(grid_dist *tt, char *name)
{
    //
    // tt: the structure returned by PDF_grid_calc
    // name: some identifier for the file. The filename will be
    // "table_PDF_name.dat"
    //
    int Nx=tt->nx;
    int NQ2=tt->nq2;
    double xmin=tt->xmin;
    double Nf=tt->nf;
    
    FILE *fd;
    
    string base_name("table_PDF_");
    string suffix(".dat");
    string filename = base_name + string(name) + suffix;
    
    fd=fopen(filename.c_str(),"w");
    fwrite(&GLUON_DERIVATIVE,sizeof(int),1,fd);
    fwrite(&Nx,sizeof(int),1,fd);
    fwrite(&NQ2,sizeof(int),1,fd);
    fwrite(&xmin,sizeof(double),1,fd);
    fwrite(&Nf,sizeof(double),1,fd);
    fwrite(tt->q2,sizeof(double),1+NQ2,fd);
    fwrite(tt->x,sizeof(double),1+Nx,fd);
    fwrite(tt->qns,sizeof(double),(1+NQ2)*(1+Nx),fd);
    fwrite(tt->qs,sizeof(double),(1+NQ2)*(1+Nx),fd);
    fwrite(tt->g,sizeof(double),(1+NQ2)*(1+Nx),fd);
    if (GLUON_DERIVATIVE) fwrite(tt->dg,sizeof(double),(1+NQ2)*(1+Nx),fd);
    
    fclose(fd);
}

//
// Read the grid from disk
//
grid_dist *PDF_grid_read(char *name)
{
    //
    // name: the name that has been used to identify the file in
    // PDF_dist_write.
    //
    int Nx;
    int NQ2;
    int der;
    double xmin,Nf;
    FILE *fd;
    
    string base_name("table_PDF_");
    string suffix(".dat");
    string filename = base_name + string(name) + suffix;
    
    fd=fopen(filename.c_str(),"r");
    fread(&der,sizeof(int),1,fd);
    fread(&Nx,sizeof(int),1,fd);
    fread(&NQ2,sizeof(int),1,fd);
    fread(&xmin,sizeof(double),1,fd);
    fread(&Nf,sizeof(double),1,fd);
    
    if (der==1) GLUON_DERIVATIVE=1;
    
    grid_dist *tt = new grid_dist(NQ2, Nx, Nf, xmin);
    
    //
    // Read the rest of the data
    //
    fread(tt->q2,sizeof(double),1+NQ2,fd);
    fread(tt->x,sizeof(double),1+Nx,fd);
    fread(tt->qns,sizeof(double),(1+NQ2)*(1+Nx),fd);
    fread(tt->qs,sizeof(double),(1+NQ2)*(1+Nx),fd);
    fread(tt->g,sizeof(double),(1+NQ2)*(1+Nx),fd);
    if (GLUON_DERIVATIVE) fread(tt->dg,sizeof(double),(1+NQ2)*(1+Nx),fd);
    
    fclose(fd);
    
    return tt;
}

//
// Get the distributions at some intermediate point (x,Q^2) by
// interpolation.
// FIXME: written as it is, this routine needs to have the complete
// table, which can be huge, loaded in memory first. One should write a
// version of this function that determines where in the file is the
// piece of data we need, and then read only the corresponding (small)
// part of the file...
//
value_dist PDF_interpolated(double x,double q2,grid_dist *tt)
{
    //
    // x: the value of x
    // q2: the value of Q^2
    // tt: the structure returned by PDF_grid_calc or by PDF_grid_read
    //
    value_dist out;
    int Nx=tt->nx;
    int NQ2=tt->nq2;
    double xmin=tt->xmin;
    double q2min=(tt->q2)[0];
    double q2max=(tt->q2)[NQ2];
    int i,j;
    double f1,f2,f3,f4,t,u;
    
    out.x=0.0;
    out.q2=0.0;
    out.qns=0.0;
    out.qs=0.0;
    out.g=0.0;
    out.dg=0.0;
    
    if ((q2<q2min)||(q2>=q2max)) {
        cout << "ERROR: q2=" << q2 << " out of range" << endl;
        return out;
    }
    
    if ((x>=208)||(x<xmin)) {
        cout << "ERROR: x=" << x << " out of range" << endl;
        return out;
    }
    
    out.x=x;
    out.q2=q2;
    
    //
    // Find the right spot in the table
    //
    i=(int)floor(Nx*log(x/xmin)/log(1.0/xmin));
    j=(int)floor(NQ2*log(q2/q2min)/log(q2max/q2min));
    
    //
    // Linear interpolation
    //
    t=(x-tt->x[i])/(tt->x[i+1]-tt->x[i]);
    u=(q2-tt->q2[j])/(tt->q2[j+1]-tt->q2[j]);
    
    f1=tt->qns[(1+NQ2)*i+j];
    f2=tt->qns[(1+NQ2)*(i+1)+j];
    f3=tt->qns[(1+NQ2)*(i+1)+j+1];
    f4=tt->qns[(1+NQ2)*i+j+1];
    
    out.qns=(1.0-t)*(1.0-u)*f1\
    +t*(1.0-u)*f2\
    +t*u*f3\
    +(1.0-t)*u*f4;
    
    f1=tt->qs[(1+NQ2)*i+j];
    f2=tt->qs[(1+NQ2)*(i+1)+j];
    f3=tt->qs[(1+NQ2)*(i+1)+j+1];
    f4=tt->qs[(1+NQ2)*i+j+1];
    
    out.qs=(1.0-t)*(1.0-u)*f1\
    +t*(1.0-u)*f2\
    +t*u*f3\
    +(1.0-t)*u*f4;
    
    f1=tt->g[(1+NQ2)*i+j];
    f2=tt->g[(1+NQ2)*(i+1)+j];
    f3=tt->g[(1+NQ2)*(i+1)+j+1];
    f4=tt->g[(1+NQ2)*i+j+1];
    
    out.g=(1.0-t)*(1.0-u)*f1\
    +t*(1.0-u)*f2\
    +t*u*f3\
    +(1.0-t)*u*f4;
    
    if (GLUON_DERIVATIVE) {
        f1=tt->dg[(1+NQ2)*i+j];
        f2=tt->dg[(1+NQ2)*(i+1)+j];
        f3=tt->dg[(1+NQ2)*(i+1)+j+1];
        f4=tt->dg[(1+NQ2)*i+j+1];
        
        out.dg=(1.0-t)*(1.0-u)*f1\
        +t*(1.0-u)*f2\
        +t*u*f3\
        +(1.0-t)*u*f4;
    }
    
    return out;
}

