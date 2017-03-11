
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
#include "data.hpp"
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>

using namespace std;

double StrToReal(std::string str)
{
    std::stringstream buff(str);
    double tmp;
    buff >> tmp;
    return tmp;
}

// Load data from given file, return 0 if no errors
int Data::LoadData(string filename, double only_charm)
{
    onlycharm = only_charm;
    
    ifstream file(filename.c_str());
    if (!file.is_open())
    {
        cerr << "ERROR! Coudn't read file " << file << endl;
        return -1;
    }
    
    while(!file.eof() )
    {
        string line;
        getline(file, line);
        if (line.substr(0, 1)=="#")
        continue;
        string x,qsqr,y,sigmar,err;
        stringstream l(line);
        l >> qsqr; l>>x; l>>y; l>>sigmar; l>>err;
        
        if (StrToReal(x)>maxx or StrToReal(x)<minx or StrToReal(qsqr)<minQ2 or StrToReal(qsqr)>maxQ2) continue;
        Qsqrvals.push_back(StrToReal(qsqr)); xbjvals.push_back(StrToReal(x));
        yvals.push_back(StrToReal(y));
        sigmarvals.push_back(StrToReal(sigmar)); errors.push_back(StrToReal(err));
    }
    
    cout << "# Loaded " << sigmarvals.size() << " datapoints from " << filename << " in Q2 range " << minQ2 << " - " << maxQ2 << " GeV^2" << endl;
    
    return 0;
}
const int Data::NumOfPoints()
{
    return sigmarvals.size();
}
                   
const double Data::Qsqr(unsigned int n)
{
    return Qsqrvals[n];
}
const double Data::xbj(unsigned int n)
{
    return xbjvals[n];
}
const double Data::y(unsigned int n)
{
    return yvals[n];
}
const double Data::ReducedCrossSection(int n)
{
    return sigmarvals[n];
}
const double Data::ReducedCrossSectionError(int n)
{
    return errors[n];
}

Data::Data()
{
    minx=1e-99;
    maxx=0.01;
    minQ2=4;
    maxQ2=20;
}
