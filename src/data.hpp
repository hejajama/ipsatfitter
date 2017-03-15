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

class Data
{
public:
    Data();
    int LoadData(string filename, double only_charm);
    int NumOfPoints() const;
    double Qsqr(unsigned int n) const;
    double xbj(unsigned int n) const;
    double y(unsigned int n) const;
    double ReducedCrossSection(unsigned int n) const;
    double ReducedCrossSectionError(unsigned int n) const;
    bool OnlyCharm(unsigned int n) const;
    
    void SetMinQsqr(double q2) { minQ2 = q2; }
    void SetMaxQsqr(double q2) { maxQ2 = q2; }
    void SetMinX(double x) { minx = x; }
    void SetMaxX(double x) { maxx = x; }
    void SetWeight(double w) { weight = w; }
    double Weight() const { return weight; }
    
    
    
private:
    
    double weight;  // Individual weight given to this dataset, default=1
    
    double minx, maxx, minQ2, maxQ2;
    
    vector<double> Qsqrvals;
    vector<double> xbjvals;
    vector<double> yvals;
    vector<double> sigmarvals;
    vector<double> errors;
    vector<bool> only_charm;
    
};

#endif
