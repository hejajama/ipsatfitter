#include "AlphaStrong.h"
#include <iostream>
#include <iomanip>

using namespace std;

int main()
{
    int order = 0;
    AlphaStrong alphaHandler(order, 1, 1, 0.5, 1.4, 4.75, 1e10);
    cout << "Order =" << order << endl;
    for (int i=1; i<200; i+=2) {
        cout << "alpha_s(" << i << " GeV) = " << setprecision(17) << alphaHandler.value(static_cast<double>(i)) << endl;
    }
    return 0;
}

