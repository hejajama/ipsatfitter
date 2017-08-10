// Last modified on February 5th, 2010 -- F. Gelis


// prototypes of the data structures and functions of dglap.c that are
// meant to be used from another program. Their definition is
// explained in the comments of the file dglap.c

typedef struct {
  int nmax;
  double q2;
  double *qns;
  double *qs;
  double *g;
} tab_dist;

typedef struct {
  int nq2;
  double x;
  double nf;
  double *q2;
  double *qns;
  double *qs;
  double *g;
  double *dg;
} tab_dist_x;

typedef struct {
  int nq2;
  int nx;
  double nf;
  double xmin;
  double *q2;
  double *x;
  double *qns;
  double *qs;
  double *g;
  double *dg;
} grid_dist;

typedef struct {
  double x;
  double q2;
  double qns;
  double qs;
  double g;
  double dg;
} value_dist;

typedef struct {
  int nmax;          
  int nq2;          
  double nf;
  double *q2;    
  double *qns;    
  double *qs;   
  double *g;          
} tab_evol;

typedef struct {
  int nq2;
  double *q2;
  double *alpha;
} tab_alpha;

typedef struct {
  tab_evol *pdf;
  tab_alpha *alpha;
} evol_pair;


double Lag_coeff(int n,double (*f)(double));
double Lag_inv(double x,int N,double *coeff);
double diff_Lag_inv(double x,int N,double *coeff);

tab_dist *Lag_dist(int N,double (*qns)(double),double (*qs)(double),double (*g)(double));

value_dist Lag_dist_inv(double x,tab_dist *tt);

double alpha_s_LO(double Q2,double Nf);
tab_alpha *alpha_s_evol(double,double,double,double,int,double);
double alpha_s(double,tab_alpha);

void set_LO(void);
void set_NO_QUARKS(void);
void set_NO_RUNNING(void);
void compute_GLUON_DERIVATIVE(void);


evol_pair DGLAP_evol(tab_dist *dist0,double Q2max,int NQ2,tab_Pij *pij, double refScale, double refValue);

tab_dist *get_step(int j,tab_evol *tt_evol);
tab_dist *PDF_at_q2(double q2,tab_evol *tt_evol);

tab_dist_x *PDF_at_x(double x,tab_evol *tt_evol);
grid_dist *PDF_grid_calc(int Nx, double xmin,tab_evol *tt_evol);

void PDF_grid_write(grid_dist *tt,char *name);
grid_dist *PDF_grid_read(char *name);

value_dist PDF_interpolated(double x,double q2,grid_dist *tt);

void free_grid_dist(grid_dist *tt);
void free_tab_evol(tab_evol *tt);
