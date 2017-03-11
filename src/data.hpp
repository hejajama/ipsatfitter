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
    const int NumOfPoints();
    const double Qsqr(unsigned int n);
    const double xbj(unsigned int n);
    const double y(unsigned int n);
    const double ReducedCrossSection(int n);
    const double ReducedCrossSectionError(int n);
    
    const bool OnlyCharm() { return onlycharm; }
    
private:
    bool onlycharm;
    
    double minx, maxx, minQ2, maxQ2;
    
    vector<double> Qsqrvals;
    vector<double> xbjvals;
    vector<double> yvals;
    vector<double> sigmarvals;
    vector<double> errors;
    
};

#endif
