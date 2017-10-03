#include "AlphaStrong.h"
#include "EvolutionLO.h"
#include <iostream>
#include <iomanip>

using namespace std;

int main()
{
    int order = 0;
    AlphaStrong *myAlphas = new AlphaStrong(order, 1, 1, 0.5, 1.4, 4.75, 1e10);
    
    EvolutionLO myEvol(myAlphas);
    
    for (int i=1; i <=1000; i+=10) {
        for (int j=1; j <=1000; j+=10) {
            
            double x = i/2000.;
            double Q2 = j;
            double mu0 = 1;
            double Ag =  2.29831;
            double lambdag = 0.09708;
            double As = 0;
            double lambdas = 0;
            int coupling = 0;
    
            double value = myEvol.alphasxG(x, Q2, mu0, coupling, Ag, lambdag,  As,  lambdas);
            cout << setprecision(5) << "x=" << x << ", Q2=" << Q2 << setprecision(12)
                 << ", value=" << value << endl;
        }
    }
    return 0;
}
