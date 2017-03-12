/*
 * C++ wrapper to GSL 2d interpolator
 * Heikki Mäntysaari <mantysaari@bnl.gov>, 2017
 */


#ifndef _INTERPOLATION2D_H
#define _INTERPOLATION2D_H



#include <gsl/gsl_bspline.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <vector>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

using namespace std;

class Interpolator2D
{
public:
    // data format [xind][yind]
    // Grid is rectangular, that is, grid(x,y) = (grid[x], grid[y])
    /**
     * Initialize 2D interpolator
     * Assumes that the interpolation grid is rectangular, and
     * coordinates are given in the grid vector.
     * Grid point corresponding to indexes (x,y) is
     * grid(x,y) = (grid[x], grid[y])
     * @param grid Vector of x and y coordinates
     * @param data 2D vector of datapoints, data[i][j] = f(grid[i], grid[j])
     * @see Interpolator
     */
    Interpolator2D();
    Interpolator2D(Interpolator2D& inter);
    // values is xpoints*ypoints array, dtapoint = values[yinx*xpoints + xind]
    void Initialize(vector<double> xpoints, vector<double> ypoints, vector<double> values);
    ~Interpolator2D();
    void Clear();
    
    double Evaluate(double x, double y) const;
    
    
private:
    gsl_spline2d *spline;
    gsl_interp_accel *xacc;
    gsl_interp_accel *yacc;
    
    bool initialized;
    
    


};




#endif
