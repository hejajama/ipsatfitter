
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
#include <cmath>
#include <vector>

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
        l >> qsqr; l>>x; l>>y; l>>sigmar; l>>err;
        
        if (StrToReal(x)>maxx or StrToReal(x)<minx or StrToReal(qsqr)<minQ2 or StrToReal(qsqr)>maxQ2) continue;
        points++;
        Qsqrvals.push_back(StrToReal(qsqr)); xbjvals.push_back(StrToReal(x));
        yvals.push_back(StrToReal(y));
        sigmarvals.push_back(StrToReal(sigmar)); errors.push_back(StrToReal(err));
        point_type.push_back(type);
        observables.push_back(SIGMAR);
        weights.push_back(weight);
    }
    
    cout << "# Loaded " << points << " datapoints from " << filename << " with weight " << weight << " in Q2 range " << minQ2 << " - " << maxQ2 << " GeV^2, now we have in total " << sigmarvals.size() << " points " << endl;
    
    return 0;
}
// Split line at commas and remove surrounding quotes and whitespace
static vector<string> SplitCSV(const string& line)
{
    vector<string> res;
    stringstream ss(line);
    string field;
    while (getline(ss, field, ','))
    {
        size_t b = field.find_first_not_of(" \t\r\n\"");
        size_t e = field.find_last_not_of(" \t\r\n\"");
        res.push_back(b == string::npos ? "" : field.substr(b, e-b+1));
    }
    return res;
}

int Data::LoadF2CSV(string filename, string herafile, DataType type, double weight)
{
    // Read HERA data: Q2 x y sigmar totalerror StatErr[%] Uncor[%] ...
    // If there are multiple measurements at the same (Q2,x) (different sqrt(s)),
    // use the one with the highest sqrt(s)
    struct HeraPoint { double q2, x, y, relerr, sqrts; };
    vector<HeraPoint> hera;
    ifstream hfile(herafile.c_str());
    if (!hfile.is_open())
    {
        cerr << "ERROR! Coudn't read file " << herafile << endl;
        return -1;
    }
    string line;
    while (getline(hfile, line))
    {
        if (line.empty() or line.substr(0, 1)=="#")
            continue;
        stringstream l(line);
        double q2, x, y, sigmar, toterr, staterr, uncor;
        if (!(l >> q2 >> x >> y >> sigmar >> toterr >> staterr >> uncor))
            continue;
        HeraPoint hp;
        hp.q2 = q2; hp.x = x; hp.y = y;
        // Neglect correlated systematic uncertainties
        hp.relerr = std::sqrt(staterr*staterr + uncor*uncor)/100.0;
        hp.sqrts = std::sqrt(q2/(x*y));
        hera.push_back(hp);
    }

    ifstream file(filename.c_str());
    if (!file.is_open())
    {
        cerr << "ERROR! Coudn't read file " << filename << endl;
        return -1;
    }

    // Skip everything until the header line containing column names Q2, x and F2
    int q2col=-1, xcol=-1, f2col=-1;
    while (getline(file, line))
    {
        vector<string> cols = SplitCSV(line);
        for (unsigned int i=0; i<cols.size(); i++)
        {
            if (cols[i]=="Q2") q2col=i;
            else if (cols[i]=="x") xcol=i;
            else if (cols[i]=="F2") f2col=i;
        }
        if (q2col>=0 and xcol>=0 and f2col>=0)
            break;
        q2col=xcol=f2col=-1;
    }
    if (q2col<0)
    {
        cerr << "ERROR! Could not find header \"Q2\",\"x\",\"F2\" in " << filename << endl;
        return -1;
    }

    int points=0;
    while (getline(file, line))
    {
        vector<string> cols = SplitCSV(line);
        if ((int)cols.size() <= std::max(q2col, std::max(xcol, f2col)) or cols[q2col].empty())
            continue;
        double q2 = StrToReal(cols[q2col]);
        double x = StrToReal(cols[xcol]);
        double f2 = StrToReal(cols[f2col]);

        if (x>maxx or x<minx or q2<minQ2 or q2>maxQ2) continue;

        const HeraPoint* match = NULL;
        for (unsigned int i=0; i<hera.size(); i++)
        {
            if (std::abs(hera[i].q2/q2 - 1.0) < 1e-4 and std::abs(hera[i].x/x - 1.0) < 1e-4)
            {
                if (match == NULL or hera[i].sqrts > match->sqrts)
                    match = &hera[i];
            }
        }
        if (match == NULL)
        {
            cerr << "# Warning: no HERA datapoint at Q2=" << q2 << ", x=" << x << " in " << herafile << ", skipping F2 point" << endl;
            continue;
        }

        points++;
        Qsqrvals.push_back(q2); xbjvals.push_back(x);
        yvals.push_back(match->y);
        sigmarvals.push_back(f2); errors.push_back(match->relerr * f2);
        point_type.push_back(type);
        observables.push_back(STRUCTURE_F2);
        weights.push_back(weight);
    }

    cout << "# Loaded " << points << " F2 datapoints from " << filename << " (relative uncertainties from " << herafile << ") with weight " << weight << " in Q2 range " << minQ2 << " - " << maxQ2 << " GeV^2, now we have in total " << sigmarvals.size() << " points " << endl;

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
    return yvals[n];
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
    weight = 1.0;
}
