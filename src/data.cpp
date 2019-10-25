
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
#include <cstdlib>

using namespace std;

double StrToReal(std::string str)
{
    std::stringstream buff(str);
    double tmp;
    buff >> tmp;
    return tmp;
}

// Load data from given file, return 0 if no errors
// By default weight = 1
int Data::LoadData(string filename, DataType type, double weight)
{
    
    ifstream file(filename.c_str());
    if (!file.is_open())
    {
        cerr << "ERROR! Coudn't read file " << filename << endl;
        return -1;
    }
    
    int points=0;
    
    while(!file.eof() )
    {
        string line;
        getline(file, line);
        if (line.substr(0, 1)=="#")
        continue;
        string x,qsqr,y,sigmar,err;
        stringstream l(line);
        
        // Inclusive case
        if (type != INC_DIFFRACTIVE_TOTAL)
        {
            l >> qsqr; l>>x; l>>y; l>>sigmar; l>>err;
            
            if (StrToReal(x)>maxx or StrToReal(x)<minx or StrToReal(qsqr)<minQ2 or StrToReal(qsqr)>maxQ2) continue;
            points++;
            yvals.push_back(StrToReal(y));
            errors.push_back(StrToReal(err));
            
        }
        else if (type == INC_DIFFRACTIVE_TOTAL)
        {
            string beta; string tmp;
            l >> qsqr; l>>beta; l>>x; l>>sigmar;
            l>>err;
            if (StrToReal(x)>maxx or StrToReal(x)<minx or StrToReal(qsqr)<minQ2 or StrToReal(qsqr)>maxQ2 or StrToReal(beta)<minbeta or StrToReal(beta)>maxbeta)
                continue;
            
            betavals.push_back(StrToReal(beta));
            errors.push_back(StrToReal(err)*StrToReal(sigmar)/100.0);
            points++;
            
        }
        Qsqrvals.push_back(StrToReal(qsqr)); xbjvals.push_back(StrToReal(x));
        point_type.push_back(type);
        weights.push_back(weight);
        sigmarvals.push_back(StrToReal(sigmar));
    }
    
    cout << "# Loaded " << points << " datapoints from " << filename << " with weight " << weight << " in Q2 range " << minQ2 << " - " << maxQ2 << " GeV^2, now we have in total " << sigmarvals.size() << " points " << endl;
    
    return 0;
}
int Data::NumOfPoints() const
{
    return sigmarvals.size();
}
                   
double Data::Qsqr(unsigned int n) const
{
    if (n >= NumOfPoints())
    {
        cerr << "Index " << n << " out of range!" << endl;
        exit(1);
    }
    return Qsqrvals[n];
}
 double Data::xbj(unsigned int n) const
{
    return xbjvals[n];
}
 double Data::y(unsigned int n) const
{
    if (point_type[n] == INC_DIFFRACTIVE_TOTAL)
        return 0;   // y not used in this case
    return yvals[n];
}
double Data::beta(unsigned int n) const
{
    if (point_type[n] != INC_DIFFRACTIVE_TOTAL)
        return 0;
    return betavals[n];
}
 double Data::ReducedCrossSection(unsigned int n) const
{
    return sigmarvals[n];
}
double Data::ReducedCrossSectionError(unsigned int n) const
{
    return errors[n];
}


DataType Data::DataPointType(unsigned int n) const
{
    return point_type[n];
}





Data::Data()
{
    minx=1e-99;
    maxx=0.01;
    minQ2=0.75;
    maxQ2=650;
    minbeta=0;
    maxbeta=1.0;
    weight = 1.0;
}
