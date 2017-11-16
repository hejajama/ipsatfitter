#ifndef WOODS_SAXON_H
#define WOODS_SAXON_H
/*
 * Simple Woods-Saxon distribution class
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2018
 */

class WoodsSaxon
{
public:
    WoodsSaxon(int A_);
    double W_S(double r);
    double T_A(double b);

private:
    int A;
    double w_s_normalization=1;
};



#endif

