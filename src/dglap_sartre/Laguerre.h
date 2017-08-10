//Last modified on November 15th, 2002 -- F. Gelis


// prototypes of the structures and functions of laguerre.c that are
// meant to be used from other programs


typedef struct {
    double nf;
    double nmax;
    double *pns0;
    double *pns1;
    double *psqq0;
    double *psqq1;
    double *psqg0;
    double *psqg1;
    double *psgq0;
    double *psgq1;
    double *psgg0;
    double *psgg1;
} tab_Pij;


typedef struct {
    double nf;
    double nmax;
    double *b2q;
    double *b2g;
} tab_b2;

tab_Pij *create_Lag_Pij_table(double Nf, int N);
tab_Pij *create_Lag_Pij_table_DLA(double Nf, int N);
tab_b2 *create_Lag_b2_table(double Nf, int N);

void free_Pij_table(tab_Pij *tt);
void free_b2_table(tab_b2 *tt);


