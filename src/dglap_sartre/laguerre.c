//Last modified on November 15th, 2002 -- F. Gelis

// This file is a translation into C of the original FORTRAN code by
// Laurent Schoeffel for the LO and NLO splitting functions. The only
// difference is that the numerical integration is handled by routines
// of the GNU Scientific Library.

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_sf_dilog.h>
#include <gsl/gsl_integration.h>


// Nc=3 is fixed once for all, Nf remains free

#define Nc 3.0
#define Cf 1.333333333333333333333
#define Cg 3.0
#define Tr 0.5

double epsilon=1.0e-10; // Precision for the integrations


// An internal structure used to pass extra parameters to the
// integration routines

typedef struct {
    double nf;         // number of flavors
    int nlag;          // Laguerre order
} data_pair;


// A structure that holds the Laguerre coefficients of the splitting
// functions

typedef struct {
    double nf;         // number of flavors
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
} tab_Pij;


// The same for the Wilson coefficients

typedef struct {
    double nf;
    double nmax;
    double *b2q;
    double *b2g;
} tab_b2;


// A function used in some intermediate calculations. We use the
// implementation of the GNU scientific library

double spence(double x){
    return -gsl_sf_dilog(1.0-x);
}


// A routine that performs numerical integrations between two bounds

double gaussk(double (*f)(double,void *),void *params,double a,double b,double eps){
    // f: function to integrate. It must take two arguments: the first
    // one is the integration variable, and the second one is a pointer
    // to some extra parameters (it can be a table of parameters, or a
    // structure).
    // params: a pointer to the extra parameters to pass to the function f
    // a: the lower integration limit
    // b: the upper integration limit
    // eps: a parameter that controls the precision of the integration
    gsl_integration_workspace *wksp=gsl_integration_workspace_alloc(200);
    double result,abserr;
    double epsabs=eps;
    double epsrel=eps;
    gsl_function F;
    
    F.function=f;
    F.params=params;
    
    gsl_integration_qag(&F,a,b,epsabs,epsrel,200,6,wksp,&result,&abserr);
    gsl_integration_workspace_free(wksp);
    return result;
}


// Laguerre polynomial of order n at point x

double Lag(int n, double x){
    double *lag0 = (double*)malloc((n+2)*sizeof(double)); // TU+TT
    int nn=n+1;
    double dnn;
    int i,ii;
    double a1,a2;
    
    lag0[1]=1.0;
    lag0[2]=1.0-x;
    dnn=((double)nn);
    
    double result;
    
    if(dnn<2.5) {
        result = lag0[nn];
    }
    else {
        for (ii=3;ii<=nn;ii++){
            i=ii-1;
            a1=((double)(2*i-1)-x)/((double)(i));
            a2=((double)(i-1))/((double)(i));
            lag0[ii]=a1*lag0[ii-1]-a2*lag0[ii-2];
        }
        result = lag0[nn];
    }
    free(lag0);
    
    return result;
}


// The code for the splitting functions and Wilson coefficients. Taken
// verbatim from Laurent Schoeffel's original code

double Pns0i(double zp,void *params){
    data_pair *data=(data_pair *)params;
    int Nlag=data->nlag;
    double zp0=0.0;
    double zk=exp(-zp);
    double omz=1.0-zk;
    double opz2=1.0+zk*zk;
    double pns01=Cf*(opz2*zk*Lag(Nlag,zp)-2.0*Lag(Nlag,zp0))/omz;
    
    pns01=pns01*exp(-zp);
    return pns01;
}

double Pns1i(double zp,void *params){
    data_pair *data=(data_pair *)params;
    int Nlag=data->nlag;
    double Nf=data->nf;
    double sign=-1.0;
    double zp0=0.0;
    double zk=exp(-zp);
    double zk1=1.0/(1.0+zk);
    double zk2=zk/(1.0+zk);
    double sf1=spence(zk1);
    double sf2=spence(zk2);
    double opz=1.0+zk;
    double omz=1.0-zk;
    double opz2=1.0+zk*zk;
    double dlnz=log(zk);
    double dlno=log(1.0-zk);
    double dln1=log(zk1);
    double dln2=log(zk2);
    double pf=(-2.0*opz2/omz*dlnz*dlno-(3.0/omz+2.0*zk)*dlnz\
               -opz/2.0*dlnz*dlnz-5.0*omz)*zk*Lag(Nlag,zp);
    double pa=(2.0*opz2/opz*(sf2-sf1-1.0/2.0*(dln1*dln1-dln2*dln2))\
               +2.0*opz*dlnz+4.0*omz)*zk*Lag(Nlag,zp);
    double pgfz=opz2*(dlnz*dlnz+11.0/3.0*dlnz+67.0/9.0-M_PI*M_PI/3.0);
    double pgf1=2.0*(67.0/9.0-M_PI*M_PI/3.0);
    double pg=(pgfz/omz+2.0*opz*dlnz+40.0/3.0*omz)*zk*Lag(Nlag,zp)\
    -pgf1/omz*Lag(Nlag,zp0);
    double pnffz=opz2*(-dlnz-5.0/3.0);
    double pnff1=-10.0/3.0;
    double pnf=2.0/3.0*((pnffz/omz-2.0*omz)*zk*Lag(Nlag,zp)
                        -pnff1/omz*Lag(Nlag,zp0));
    double pns11=Cf*Cf*(pf+sign*pa)+Cf*Cg/2.0*(pg-sign*pa)+Cf*Tr*Nf*pnf;
    
    pns11=pns11*exp(-zp);
    
    return pns11;
}

double Psqq0i(double zp,void *params){
    data_pair *data=(data_pair *)params;
    int Nlag=data->nlag;
    double zp0=0.0;
    double zk=exp(-zp);
    double omz=1.0-zk;
    double opz2=1.0+zk*zk;
    double psqq01=Cf*(opz2*zk*Lag(Nlag,zp)-2.0*Lag(Nlag,zp0))/omz;
    
    psqq01=psqq01*exp(-zp);
    
    return psqq01;
}

double Psqq1i(double zp,void *params){
    data_pair *data=(data_pair *)params;
    int Nlag=data->nlag;
    double Nf=data->nf;
    double sign=1.0;
    double zp0=0.0;
    double zk=exp(-zp);
    double zk1=1.0/(1.0+zk);
    double zk2=zk/(1.0+zk);
    double sf1=spence(zk1);
    double sf2=spence(zk2);
    double opz=1.0+zk;
    double omz=1.0-zk;
    double opz2=1.0+zk*zk;
    double dlnz=log(zk);
    double dlno=log(1.0-zk);
    double dln1=log(zk1);
    double dln2=log(zk2);
    double pf=(-2.0*opz2/omz*dlnz*dlno-(3.0/omz+2.0*zk)*dlnz\
               -opz/2.0*dlnz*dlnz-5.0*omz)*Lag(Nlag,zp)*zk;
    double pa=(2.0*opz2/opz*(sf2-sf1-1.0/2.0*(dln1*dln1-dln2*dln2))\
               +2.0*opz*dlnz+4.0*omz)*Lag(Nlag,zp)*zk;
    double pgfz=opz2*(dlnz*dlnz+11.0/3.0*dlnz+67.0/9.0-M_PI*M_PI/3.0);
    double pgf1=2.0*(67.0/9.0-M_PI*M_PI/3.0);
    double pg=(pgfz/omz+2.0*opz*dlnz+40.0/3.0*omz)*Lag(Nlag,zp)*zk\
    -pgf1/omz*(Lag(Nlag,zp0));
    double pnffz=opz2*(-dlnz-5.0/3.0);
    double pnff1=-10.0/3.0;
    double pnf=2.0/3.0*((pnffz/omz-2.0*omz)*Lag(Nlag,zp)*zk\
                        -pnff1/omz*Lag(Nlag,zp0));
    double pns1s=Cf*Cf*(pf+sign*pa)+Cf*Cg/2.0*(pg-sign*pa)+Cf*Tr*Nf*pnf;
    //double sfopz=spence(1.0+zk);
    // double sfj=-M_PI*M_PI/12.0-dlnz*log(opz)+sfopz;
    double pns1p=pns1s;
    double fqq=(20.0/9.0/zk-2.0+6.0*zk-56.0/9.0*zk*zk\
                +(1.0+5.0*zk+24.0/9.0*zk*zk)*dlnz\
                -opz*dlnz*dlnz)*Lag(Nlag,zp)*zk;
    double psqq11=pns1p+2.0*Cf*Tr*Nf*fqq;
    
    psqq11=psqq11*exp(-zp);
    
    return psqq11;
}

double Psqg0i(double zp,void *params){
    data_pair *data=(data_pair *)params;
    int Nlag=data->nlag;
    double Nf=data->nf;
    //double zp0=0.0;
    double zk=exp(-zp);
    double omz=1.0-zk;
    //double opz2=1.0+zk*zk;
    
    // The factor 2Nf from Laurent Schoeffel's paper is included here
    
    double psqg01=2.0*Nf*Tr*(zk*zk+omz*omz)*zk*Lag(Nlag,zp);
    
    psqg01=psqg01*exp(-zp);
    
    return psqg01;
}

double Psqg1i(double zp,void *params){
    data_pair *data=(data_pair *)params;
    int Nlag=data->nlag;
    double Nf=data->nf;
    //double zp0=0.0;
    double zk=exp(-zp);
    double omz=1.0-zk;
    double opz=1.0+zk;
    double dlnz=log(zk);
    double dlno=log(1.0-zk);
    double sfopz=spence(1.0+zk);
    double sfj=-M_PI*M_PI/12.0-dlnz*log(opz)+sfopz;
    double fqg1=(-40.0/9.0/zk+4.0-50.0*zk+436.0/9.0*zk*zk\
                 -(2.0+16.0*zk+88.0/3.0*zk*zk)*dlnz\
                 +(2.0+4.0*zk)*dlnz*dlnz\
                 +8.0*zk*omz*dlno\
                 +(2.0-4.0*zk+4.0*zk*zk)*(dlno*dlno-M_PI*M_PI/6.0)\
                 -4.0*(1.0+2.0*zk+2.0*zk*zk)*sfj)*\
    zk*Lag(Nlag,zp);
    double fqg2=(-14.0+29.0*zk-20.0*zk*zk\
                 +(-3.0+4.0*zk-8.0*zk*zk)*dlnz\
                 -(1.0-2.0*zk+4.0*zk*zk)*dlnz*dlnz\
                 -8.0*zk*omz*dlno\
                 +(2.0-4.0*zk+4.0*zk*zk)\
                 *(-dlno*dlno+2.0*dlnz*dlno+M_PI*M_PI/3.0))*\
    zk*Lag(Nlag,zp);
    double psqg11=-Nf*(Cg*Tr*fqg1+Cf*Tr*fqg2);
    
    psqg11=psqg11*exp(-zp);
    
    return psqg11;
}

double Psgq0i(double zp,void *params){
    data_pair *data=(data_pair *)params;
    int Nlag=data->nlag;
    //double zp0=0.0;
    double zk=exp(-zp);
    double omz=1.0-zk;
    double psgq01=Cf*(1.0+omz*omz)/zk*(Lag(Nlag,zp)*zk);
    
    psgq01=psgq01*exp(-zp);
    
    return psgq01;
}

double Psgq1i(double zp,void *params){
    data_pair *data=(data_pair *)params;
    int Nlag=data->nlag;
    double Nf=data->nf;
    //double zp0=0.0;
    double zk=exp(-zp);
    //double opz=1.0+zk;
    //double omz=1.0-zk;
    double dlnz=log(zk);
    double dlno=log(1.0-zk);
    double sfopz=spence(1.0+zk);
    //double sfi=spence(zk);
    double sfj=-M_PI*M_PI/12.0-dlnz*log(1.0+zk)+sfopz;
    double fgq1=(5.0/2.0+7.0/2.0*zk-(2.0+7.0/2.0*zk)*dlnz\
                 +(1.0-zk/2.0)*dlnz*dlnz+(6.0/zk-6.0+5.0*zk)*dlno\
                 +(2.0/zk-2.0+zk)*dlno*dlno)*Lag(Nlag,zp)*zk;
    double fgq2=(-1.0/zk-19.0/9.0-37.0/9.0*zk-44.0/9.0*zk*zk\
                 +(12.0+5.0*zk+8.0/3.0*zk*zk)*dlnz\
                 -(2.0+zk)*dlnz*dlnz\
                 -1.0/3.0*(22.0/zk-22.0+17.0*zk)*dlno\
                 -(2.0/zk-2.0+zk)\
                 *(dlno*dlno-2.0*dlnz*dlno-M_PI*M_PI/6.0)\
                 +(4.0/zk+4.0+2.0*zk)*sfj)*Lag(Nlag,zp)*zk;
    double fgq3=(40.0/9.0/zk-40.0/9.0+32.0/9.0*zk\
                 +4.0/3.0*(2.0/zk-2.0+zk)*dlno)*Lag(Nlag,zp)*zk;
    double psgq11=-(Cf*Cf*fgq1+Cf*Cg*fgq2+Cf*Tr*Nf*fgq3);
    
    psgq11=psgq11*exp(-zp);
    
    return psgq11;
}

double Psgg0i(double zp,void *params){
    data_pair *data=(data_pair *)params;
    int Nlag=data->nlag;
    double zp0=0.0;
    double zk=exp(-zp);
    double omz=1.0-zk;
    double psgg01;
    
    // fprintf(stderr,"Psgg0i: Nf= %.2f nlag= %d\n",Nf,Nlag);
    
    psgg01=2.0*Cg*((zk*Lag(Nlag,zp)*zk-Lag(Nlag,zp0))/omz \
                   +(omz/zk+zk*omz)*Lag(Nlag,zp)*zk);
    psgg01=psgg01*exp(-zp);
    
    return psgg01;
}

double PsggDLAi(double zp,void *params){
    data_pair *data=(data_pair *)params;
    int Nlag=data->nlag;
    double zk=exp(-zp);
    double psggDLA;
    
    psggDLA=2.0*Cg*Lag(Nlag,zp)*zk;
    
    return psggDLA;
}


double Psgg1i(double zp,void *params){
    data_pair *data=(data_pair *)params;
    int Nlag=data->nlag;
    double Nf=data->nf;
    double zp0=0.0;
    double zk=exp(-zp);
    double opz=1.0+zk;
    double omz=1.0-zk;
    double dlnz=log(zk);
    double dlno=log(1.0-zk);
    double sfopz=spence(1.0+zk);
    double sfi=spence(zk);
    double sfj=-M_PI*M_PI/12.0-dlnz*log(1.0+zk)+sfopz;
    double fgg1=(-20.0/9.0/omz\
                 -2.0/9.0*(23.0/zk-29.0+19.0*zk-23.0*zk*zk)\
                 -4.0/3.0*opz*dlnz)*Lag(Nlag,zp)*zk\
    +20.0/9.0/omz*Lag(Nlag,zp0);
    double fgg2=(4.0/3.0/zk-16.0+8.0*zk+20.0/3.0*zk*zk\
                 -(6.0+10.0*zk)*dlnz-2.0*opz*dlnz*dlnz)*\
    Lag(Nlag,zp)*zk;
    double fgg3=(67.0/9.0/omz-8.0/9.0-zk/18.0\
                 -(47.0/6.0-25.0/6.0*zk+44.0/3.0*zk*zk)*dlnz\
                 +(1.0/omz+1.0/opz+5.0/4.0\
                   +6.0*zk-2.0*zk*zk)*dlnz*dlnz\
                 +(-11.0/2.0/zk+7.0/2.0+13.0/4.0*zk)*dlno\
                 +1.0/4.0*(-5.0/zk+5.0+2.0*zk)*dlno*dlno\
                 -(4.0/omz+4.0/zk-5.0+7.0*zk-4.0*zk*zk)*dlnz*dlno\
                 +opz*sfi+4.0*(1.0/opz-1.0/zk-2.0-zk-zk*zk)*sfj\
                 -M_PI*M_PI/3.0*(1.0/omz+1.0/zk-1.0+2.0*zk-zk*zk))*\
    Lag(Nlag,zp)*zk\
    +(-67.0/9.0+M_PI*M_PI/3.0)/omz*Lag(Nlag,zp0);
    double psgg11=Cg*Tr*Nf*fgg1+Cf*Tr*Nf*fgg2+Cg*Cg*fgg3;
    
    psgg11=psgg11*exp(-zp);
    
    return  psgg11;
}

double b2qi(double zp,void *params){
    data_pair *data=(data_pair *)params;
    int Nlag=data->nlag;
    double zp0=0.00;
    double zk=exp(-zp);
    double omz=1.0-zk;
    double opz2=1.0+zk*zk;
    double dlnz=log(zk);
    double dlno=log(1.0-zk);
    double fq=(opz2/omz*(-3.0/2.0-2.0*dlnz+2.0*dlno)\
               +(9.0+5.0*zk)/2.0)*zk*Lag(Nlag,zp)\
    +(3.0-4.0*dlno)/omz*Lag(Nlag,zp0);
    double b2q1=Cf*fq*exp(-zp);
    
    return  b2q1;
}

double b2gi(double zp,void *params){
    data_pair *data=(data_pair *)params;
    int Nlag=data->nlag;
    double zk=exp(-zp);
    double omz=1.0-zk;
    double dlnz=log(zk);
    double dlno=log(1.0-zk);
    double fg=((1.0-2.0*zk+2.0*zk*zk)*(dlno-dlnz)\
               -1.0+8.0*zk*omz)*zk*Lag(Nlag,zp);
    double b2g1=4.0*Tr*fg*exp(-zp);
    
    return b2g1;
}


double Pns0(int n,double Nf){
    int Nlag=n;
    double zpmax=100.0;
    double zp0=0.0;
    //double zta=1.2020569031595943;
    //double tils=-5.0*zta;
    double dlnomx=0.0;
    double pns0=0.0;
    double pns0da;
    data_pair params;
    
    params.nf=Nf;
    params.nlag=Nlag;
    
    pns0=gaussk(Pns0i,&params,zp0,zpmax,epsilon);
    pns0da=1.5*Cf*(1.0+1.33333333333333333333*dlnomx);
    pns0=pns0+pns0da*Lag(Nlag,zp0);
    
    return pns0;
}


double Pns1(int n,double Nf){
    int Nlag=n;
    double zpmax=100.0;
    double zp0=0.0;
    double zta=1.20205690315959430;
    double tils=-5.0*zta;
    double dlnomx=0.00;
    double pns1;
    double dfun1=Cf*Cf*(3.0/8.0-M_PI*M_PI/2.0+zta-tils)\
    +1.0/2.0*Cf*Cg*(17.0/12.0+11.0*M_PI*M_PI/9.0-zta+tils)\
    -Cf*Tr*Nf*(1.0/6.0+2.0*M_PI*M_PI/9.0);
    double pgap=Cf*Cg*(67.0/9.0-M_PI*M_PI/3.0)*dlnomx;
    double pnfap=Cf*Tr*Nf*2.0/3.0*(-10.0/3.0)*dlnomx;
    data_pair params;
    
    params.nf=Nf;
    params.nlag=Nlag;
    
    pns1=gaussk(Pns1i,&params,zp0,zpmax,epsilon);
    pns1=pns1+(dfun1+pgap+pnfap)*Lag(Nlag,zp0);
    
    return pns1;
}

double Psqq0(int n,double Nf){
    int Nlag=n;
    double zpmax=100.0;
    double zp0=0.0;
    double psqq0=0.00;
    double pqq0da=3.0/2.0*Cf;
    data_pair params;
    
    params.nf=Nf;
    params.nlag=Nlag;
    
    psqq0=gaussk(Psqq0i,&params,zp0,zpmax,epsilon);
    psqq0=psqq0+pqq0da*Lag(Nlag,zp0);
    
    return  psqq0;
}
double Psqq1(int n,double Nf){
    int Nlag=n;
    double zpmax=100.00;
    double zp0=0.00;
    double zta=1.20205690315959430;
    double tils=-5.0*zta;
    double dlnomx=0.00;
    double psqq1=0.00;
    double dfun1=Cf*Cf*(3.0/8.0-M_PI*M_PI/2.0+zta-tils)\
    +1.0/2.0*Cf*Cg*(17.0/12.0+11.0*M_PI*M_PI/9.0-zta+tils)\
    -Cf*Tr*Nf*(1.0/6.0+2.0*M_PI*M_PI/9.0);
    double pgap=Cf*Cg*(67.0/9.0-M_PI*M_PI/3.0)*dlnomx;
    double pnfap=Cf*Tr*Nf*2.0/3.0*(-10.0/3.0)*dlnomx;
    data_pair params;
    
    params.nf=Nf;
    params.nlag=Nlag;
    
    psqq1=gaussk(Psqq1i,&params,zp0,zpmax,epsilon);
    psqq1=psqq1+(dfun1+pgap+pnfap)*Lag(Nlag,zp0);
    
    return psqq1;
}

double Psqg0(int n,double Nf){
    int Nlag=n;
    double zpmax=100.00;
    double zp0=0.00;
    //double zta=1.20205690315959430;
    //double tils=-5.0*zta;
    //double dlnomx=0.00;
    double psqg0=0.00;
    data_pair params;
    
    params.nf=Nf;
    params.nlag=Nlag;
    
    psqg0=gaussk(Psqg0i,&params,zp0,zpmax,epsilon);
    
    return psqg0;
}

double Psqg1(int n,double Nf){
    int Nlag=n;
    double zpmax=100.00;
    double zp0=0.00;
    //double zta=1.20205690315959430;
    //double tils=-5.0*zta;
    //double dlnomx=0.00;
    double psqg1=0.00;
    data_pair params;
    
    params.nf=Nf;
    params.nlag=Nlag;
    
    psqg1=gaussk(Psqg1i,&params,zp0,zpmax,epsilon);
    
    return psqg1;
}

double Psgq0(int n,double Nf){
    int Nlag=n;
    double zpmax=100.00;
    double zp0=0.00;
    //double zta=1.20205690315959430;
    //double tils=-5.0*zta;
    //double dlnomx=0.00;
    double psgq0=0.00;
    data_pair params;
    
    params.nf=Nf;
    params.nlag=Nlag;
    
    psgq0=gaussk(Psgq0i,&params,zp0,zpmax,epsilon);
    
    return psgq0;
}

double Psgq1(int n,double Nf){
    int Nlag=n;
    double zpmax=100.00;
    double zp0=0.00;
    //double zta=1.20205690315959430;
    //double tils=-5.0*zta;
    //double dlnomx=0.00;
    double psgq1=0.00;
    data_pair params;
    
    params.nf=Nf;
    params.nlag=Nlag;
    
    psgq1=gaussk(Psgq1i,&params,zp0,zpmax,epsilon);
    
    return psgq1;
}

double Psgg0(int n,double Nf){
    int Nlag=n;
    double zpmax=100.0;
    double zp0=0.0;
    double psgg0=0.0;
    double pgg0da;
    data_pair params;
    
    params.nf=Nf;
    params.nlag=Nlag;
    
    psgg0=gaussk(Psgg0i,&params,zp0,zpmax,epsilon);
    pgg0da=2.0*Cg*(11.0/12.0-1.0/3.0*Nf*Tr/Cg);
    psgg0=psgg0+pgg0da*Lag(Nlag,zp0);
    
    return psgg0;
}

double PsggDLA(int n,double Nf){
    int Nlag=n;
    double zpmax=100.0;
    double zp0=0.0;
    double psggDLA;
    data_pair params;
    
    params.nf=Nf;
    params.nlag=Nlag;
    
    psggDLA=gaussk(PsggDLAi,&params,zp0,zpmax,epsilon);
    
    return psggDLA;
}

double Psgg1(int n,double Nf){
    int Nlag=n;
    double zpmax=100.00;
    double zp0=0.00;
    double zta=1.20205690315959430;
    double tils=-5.0*zta;
    double dlnomx=0.00;
    double psgg1=0.00;
    double fgg1da=-Cg*Tr*Nf*(20.0/9.0*dlnomx+4.0/3.0);
    double fgg2df=-Cf*Tr*Nf;
    double fgg3da=Cg*Cg*((67.0/9.0-M_PI*M_PI/3.0)*dlnomx\
                         -tils/2.0+zta/2.0+8.0/3.0);
    data_pair params;
    
    params.nf=Nf;
    params.nlag=Nlag;
    
    psgg1=gaussk(Psgg1i,&params,zp0,zpmax,epsilon);
    psgg1=psgg1+(fgg1da+fgg2df+fgg3da)*Lag(Nlag,zp0);
    
    return psgg1;
}

double b2q(int n,double Nf){
    int Nlag=n;
    double zpmax=100.00;
    double zp0=0.00;
    double dlnomx=0.00;
    double b2q;
    double b2qdap=Cf*(-9.0-2.0/3.0*M_PI*M_PI\
                      -3.0*dlnomx+2.0*dlnomx*dlnomx);
    data_pair params;
    
    params.nf=Nf;
    params.nlag=Nlag;
    
    b2q=gaussk(b2qi,&params,zp0,zpmax,epsilon);
    b2q=b2q+b2qdap*Lag(Nlag,zp0);
    
    return b2q;
}

double b2g(int n,double Nf){
    int Nlag=n;
    double zpmax=100.00;
    double zp0=0.00;
    double b2g=0.0;
    data_pair params;
    
    params.nf=Nf;
    params.nlag=Nlag;
    
    b2g=gaussk(b2gi,&params,zp0,zpmax,epsilon);
    
    return b2g;
}


// A routine that calculates the Laguerre coefficients of the
// splitting functions. It returns a pointer to a structure of type
// tab_Pij, with properly allocated tables.

tab_Pij *create_Lag_Pij_table(double Nf, int N){
    tab_Pij *tt=(tab_Pij *)malloc(sizeof(tab_Pij));
    int i;
    
    fprintf(stderr,"Calculating Laguerre transform of the Pij...");
    
    // Allocate some memory
    
    tt->pns0=(double *)malloc((1+N)*sizeof(double));
    tt->pns1=(double *)malloc((1+N)*sizeof(double));
    tt->psqq0=(double *)malloc((1+N)*sizeof(double));
    tt->psqq1=(double *)malloc((1+N)*sizeof(double));
    tt->psqg0=(double *)malloc((1+N)*sizeof(double));
    tt->psqg1=(double *)malloc((1+N)*sizeof(double));
    tt->psgq0=(double *)malloc((1+N)*sizeof(double));
    tt->psgq1=(double *)malloc((1+N)*sizeof(double));
    tt->psgg0=(double *)malloc((1+N)*sizeof(double));
    tt->psgg1=(double *)malloc((1+N)*sizeof(double));
    
    tt->nf=Nf;
    tt->nmax=N;
    
    // Calculate the Laguerre coefficients
    
    for(i=0;i<=N;i++){
        (tt->pns0)[i]=Pns0(i,Nf);
        (tt->pns1)[i]=Pns1(i,Nf);
        (tt->psqq0)[i]=Psqq0(i,Nf);
        (tt->psqq1)[i]=Psqq1(i,Nf);
        (tt->psqg0)[i]=Psqg0(i,Nf);
        (tt->psqg1)[i]=Psqg1(i,Nf);
        (tt->psgq0)[i]=Psgq0(i,Nf);
        (tt->psgq1)[i]=Psgq1(i,Nf);
        (tt->psgg0)[i]=Psgg0(i,Nf);
        (tt->psgg1)[i]=Psgg1(i,Nf);
    }
    
    fprintf(stderr,"done\n");
    
    return tt;
}


// Same routine for the DLA approximation of DGLAP

tab_Pij *create_Lag_Pij_table_DLA(double Nf, int N){
    tab_Pij *tt=(tab_Pij *)malloc(sizeof(tab_Pij));
    int i;
    
    fprintf(stderr,"Calculating Laguerre transform of the Pij...");
    
    // Allocate some memory
    
    tt->pns0=(double *)malloc((1+N)*sizeof(double));
    tt->pns1=(double *)malloc((1+N)*sizeof(double));
    tt->psqq0=(double *)malloc((1+N)*sizeof(double));
    tt->psqq1=(double *)malloc((1+N)*sizeof(double));
    tt->psqg0=(double *)malloc((1+N)*sizeof(double));
    tt->psqg1=(double *)malloc((1+N)*sizeof(double));
    tt->psgq0=(double *)malloc((1+N)*sizeof(double));
    tt->psgq1=(double *)malloc((1+N)*sizeof(double));
    tt->psgg0=(double *)malloc((1+N)*sizeof(double));
    tt->psgg1=(double *)malloc((1+N)*sizeof(double));
    
    tt->nf=Nf;
    tt->nmax=N;
    
    // Calculate the Laguerre coefficients
    
    for(i=0;i<=N;i++){
        (tt->pns0)[i]=0.0;
        (tt->pns1)[i]=0.0;
        (tt->psqq0)[i]=0.0;
        (tt->psqq1)[i]=0.0;
        (tt->psqg0)[i]=0.0;
        (tt->psqg1)[i]=0.0;
        (tt->psgq0)[i]=0.0;
        (tt->psgq1)[i]=0.0;
        (tt->psgg0)[i]=PsggDLA(i,Nf);
        (tt->psgg1)[i]=0.0;
    }
    
    fprintf(stderr,"done\n");
    
    return tt;
}


// The same for the Wilson coefficients

tab_b2 *create_Lag_b2_table(double Nf, int N){
    tab_b2 *tt=(tab_b2 *)malloc(sizeof(tab_b2));
    int i;
    
    fprintf(stderr,"Calculating Laguerre transform of Wilson coefficients...");
    
    tt->b2q=(double *)malloc((1+N)*sizeof(double));
    tt->b2g=(double *)malloc((1+N)*sizeof(double));
    
    tt->nf=Nf;
    tt->nmax=N;
    
    for(i=0;i<=N;i++){
        (tt->b2q)[i]=b2q(i,Nf);
        (tt->b2g)[i]=b2g(i,Nf);
    }
    
    fprintf(stderr,"done\n");
    
    return tt;
}


// free the memory allocated for a structure of type tab_Pij

void free_Pij_table(tab_Pij *tt){
    free(tt->pns0);
    free(tt->pns1);
    free(tt->psqq0);
    free(tt->psqq1);
    free(tt->psqg0);
    free(tt->psqg1);
    free(tt->psgq0);
    free(tt->psgq1);
    free(tt->psgg0);
    free(tt->psgg1);
    
    free(tt);
}


// free the memory allocated for a structure of type tab_b2

void free_b2_table(tab_b2 *tt){
    free(tt->b2q);
    free(tt->b2g);
    
    free(tt);
}





