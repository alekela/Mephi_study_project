#include <iostream>
#include <float.h>
#include <cmath>

using namespace std;

int main() {
    float eps;
    eps = 1;
    while (1 + eps > 1){
        eps = eps / 2;
    }
    cout << "Machine epsilon, float: " << eps * 2 << "\n";
    cout << "Theoretical machine epsilon, float: " << FLT_EPSILON << "\n\n";

    double eps2;
    eps2 = 1;
    while (1 + eps2 > 1) {
        eps2 = eps2 / 2;
    }
    cout << "Machine epsilon, double: " << eps2 * 2 << "\n";
    cout << "Theoretical machine epsilon, double: " << DBL_EPSILON << "\n\n";

    long double eps3;
    eps3 = 1;
    while (1 + eps3 > 1) {
        eps3 = eps3 / 2;
    }
    cout << "Machine epsilon, long double: " << eps3 * 2 << "\n";
    cout << "Theoretical machine epsilon, long double: " << LDBL_EPSILON << "\n-----------------------------------\n";


    float dwarf, real_dwarf;
    dwarf = 1;
    while (dwarf > 0){
        real_dwarf = dwarf;
        dwarf = dwarf / 2;
    }
    cout << "Machine dwarf, float: " << real_dwarf << "\n";
    cout << "Theoretical machine min, float: " << FLT_MIN << "\n\n";

    double dwarf2, real_dwarf2;
    dwarf2 = 1;
    while (dwarf2 > 0){
        real_dwarf2 = dwarf2;
        dwarf2 = dwarf2 / 2;
    }
    cout << "Machine dwarf, double: " << real_dwarf2 << "\n";
    cout << "Theoretical machine min, double: " << DBL_MIN << "\n\n";

    long double dwarf3, real_dwarf3;
    dwarf3 = 1;
    while (dwarf3 > 0){
        real_dwarf3 = dwarf3;
        dwarf3 = dwarf3 / 2;
    }
    cout << "Machine dwarf, long double: " << real_dwarf3 << "\n";
    cout << "Theoretical machine min, long double: " << LDBL_MIN << "\n-----------------------------------\n";

    float maximum, pre_maximum;
    maximum = 1;
    while (maximum != (float) INFINITY){
        pre_maximum = maximum;
        maximum = maximum * 2;
    }
    cout << "Maximum, float: " << pre_maximum << "\n";
    cout << "Theoretical maximum, float: " << FLT_MAX << "\n\n";

    double maximum2, pre_maximum2;
    maximum2 = 1;
    while (maximum2 != (double) INFINITY){
        pre_maximum2 = maximum2;
        maximum2 = maximum2 * 2;
    }
    cout << "Maximum, double: " << pre_maximum2 << "\n";
    cout << "Theoretical maximum, double: " << DBL_MAX << "\n\n";

    long double maximum3, pre_maximum3;
    maximum3 = 1;
    while (maximum3 != (long double) INFINITY){
        pre_maximum3 = maximum3;
        maximum3 = maximum3 * 2;
    }
    cout << "Maximum, long double: " << pre_maximum3 << "\n";
    cout << "Theoretical maximum, long double: " << LDBL_MAX << "\n\n";
    
}
