
#include "utils.h"

using namespace std;
using namespace NTL;


void print_poly (ZZ_pX& P)
{
    long degree = deg(P);
    if (-1 == degree) {
        cout << "0";
        return;
    }
    for(long i=0; i<=degree; i++) {
        cout << coeff(P, i);
        if (i==1)
            cout << "X";
        else if(i>1)
            cout << "X^" << i;
        if (i<degree) {
            cout << " + ";
        }
    }
//    cout << endl << "random poly:" << endl << P << endl;
}
