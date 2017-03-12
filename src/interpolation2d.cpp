/*
 * C++ wrapper to GSL 2d interpolator
 * Heikki Mäntysaari <mantysaari@bnl.gov>, 2017
 */

#include <gsl/gsl_bspline.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>
#include <cmath>
#include <iostream>
#include <gsl/gsl_spline2d.h>
#include "interpolation2d.hpp"

using namespace std;
typedef  unsigned int uint;

#ifndef LINEINFO
    #define LINEINFO __FILE__ << ":" << __LINE__
#endif

double Interpolator2D::Evaluate(double x, double y) const
{
    return gsl_spline2d_eval(spline, x, y, xacc, yacc);
}

// values is xpoints*ypoints array, dtapoint = values[yinx*xpoints + xind]
void Interpolator2D::Initialize(vector<double> xpoints, vector<double> ypoints, vector<double> values)
{
    Clear();
    const gsl_interp2d_type *T = gsl_interp2d_bilinear;
    int nx = xpoints.size();
    int ny = ypoints.size();
    spline = gsl_spline2d_alloc(T, nx, ny);
    
    xacc = gsl_interp_accel_alloc();
    yacc = gsl_interp_accel_alloc();
    
    double *xa = new double[nx];
    double *ya = new double[ny];
    for (int i=0; i<nx; i++) xa[i] = xpoints[i];
    for (int i=0; i<ny; i++) ya[i] = ypoints[i];
    
    // set datapoints
    double *data = new double[nx*ny];
    
    for (int yind=0; yind<ny; yind++)
    {
        for (int xind=0; xind<nx; xind++)
        {
            gsl_spline2d_set(spline, data, xind, yind, values[yind*nx+xind]);
        }
    }
    gsl_spline2d_init(spline, xa, ya, data, nx, ny);
    initialized=true;
    
    delete[] data;
    delete[] xa;
    delete[] ya;
    
    
    
}

void Interpolator2D::Clear()
{
    if (initialized)
    {
        gsl_spline2d_free(spline);
        gsl_interp_accel_free(xacc);
        gsl_interp_accel_free(yacc);
        initialized = false;
    }
}


// Initialize, data format: [xind][yind]
// Grid is rectangular, that is, grid(x,y) = (grid[x], grid[y])
Interpolator2D::Interpolator2D()
{
    initialized=false;
}

Interpolator2D::~Interpolator2D()
{
    Clear();

}


// Copy data from given class and initialize this, as this is
// the copy constructor
Interpolator2D::Interpolator2D(Interpolator2D& inter)
{
    cerr << "Interpolator2D copy constructor may not work???? " << LINEINFO << endl;
}
