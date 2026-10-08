#ifndef HERADATA_H_
#define HERADATA_H_

/*
 * Loads HERA DIS data on reduced cross section
 * Keeps also track whether the data is for reduced cross section
 * or for only charm
 *
 * Datafile format is
 * Q^2 x y sigmared error
 *
 * Heikki Mäntysaari <mantysaari@bnl.gov>, 2017
 */

#include <vector>
#include <string>
#include "wave_function.hpp"
using namespace std;

enum DataType
{
    TOTAL,      // total reduced cross section
    CHARM,       // only charm contribution
    BOTTOM,      // only b
    UDS          // Light quarks u,d,s
};

enum Observable
{
    SIGMAR,     // reduced cross section
    STRUCTURE_F2    // structure function F2
};

class Data
{
public:
    Data();
    int LoadData(string filename, DataType type, double weight = 1.0);
    // Load F2 values from a csv file (header "Q2","x","F2"), uncertainties are
    // taken as the relative uncorrelated (stat+uncor) uncertainty of the HERA
    // reduced cross section measurement in herafile at the same (Q2,x)
    int LoadF2CSV(string filename, string herafile, DataType type, double weight = 1.0);
    int NumOfPoints() const;
    double Qsqr(unsigned int n) const;
    double xbj(unsigned int n) const;
    double y(unsigned int n) const;
    double ReducedCrossSection(unsigned int n) const;
    double ReducedCrossSectionError(unsigned int n) const;
    DataType DataPointType(unsigned int n) const;
    Observable PointObservable(unsigned int n) const { return observables[n]; }
    
    void SetMinQsqr(double q2) { minQ2 = q2; }
    void SetMaxQsqr(double q2) { maxQ2 = q2; }
    void SetMinX(double x) { minx = x; }
    void SetMaxX(double x) { maxx = x; }
    double Weight(unsigned int i) const { return weights[i]; }
    
    
    
private:
    
    double weight;  // Individual weight given to this dataset, default=1
    
    double minx, maxx, minQ2, maxQ2;
    
    vector<double> Qsqrvals;
    vector<double> xbjvals;
    vector<double> yvals;
    vector<double> sigmarvals;
    vector<double> errors;
    vector<DataType> point_type;
    vector<Observable> observables;
    vector<double> weights;
    
};

double StrToReal(std::string str);

#endif
