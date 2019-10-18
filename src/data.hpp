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
    UDS,          // Light quarks u,d,s
    INC_DIFFRACTIVE_TOTAL,  // Total inclusive diffraction
};

class Data
{
public:
    Data();
    int LoadData(string filename, DataType type, double weight = 1.0);
    int NumOfPoints() const;
    double Qsqr(unsigned int n) const;
    double xbj(unsigned int n) const;
    double y(unsigned int n) const;
    double beta(unsigned int n) const;
    double ReducedCrossSection(unsigned int n) const;
    double ReducedCrossSectionError(unsigned int n) const;
    DataType DataPointType(unsigned int n) const;
    
    void SetMinQsqr(double q2) { minQ2 = q2; }
    void SetMaxQsqr(double q2) { maxQ2 = q2; }
    void SetMinX(double x) { minx = x; }
    void SetMaxX(double x) { maxx = x; }
    void SetMinBeta(double b) { minbeta = b; }
    void SetMaxBeta(double b) { maxbeta = b; }
    double Weight(unsigned int i) const { return weights[i]; }
    
    
    
private:
    
    double weight;  // Individual weight given to this dataset, default=1
    
    double minx, maxx, minQ2, maxQ2, maxbeta, minbeta;
    
    vector<double> Qsqrvals;
    vector<double> xbjvals;
    vector<double> yvals;       // used in structure functions
    vector<double> betavals;    // used in inclusive diffraction
    vector<double> sigmarvals;
    vector<double> errors;
    vector<DataType> point_type;
    vector<double> weights;
    
};

double StrToReal(std::string str);

#endif
