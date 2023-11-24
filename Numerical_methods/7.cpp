#include <iostream>

using namespace std;


void tridiag_matrix_algorithm(double* a, double* b, double* c, double* d, int n) {
    c[0] /= b[0];
    d[0] /= b[0];

    for (int i = 1; i < n; i++) {
        c[i] /= b[i] - a[i]*c[i-1];
        d[i] = (d[i] - a[i]*d[i-1]) / (b[i] - a[i]*c[i-1]);
    }

    for (int i = n - 1; i >= 0; i--) {
        d[i] -= c[i]*d[i+1];
    }
}


int main() {
	int  n = 4;
	double a[4] = { 0, -1, -1, -1 };
	double b[4] = { 4,  4,  4,  4 };
	double c[4] = {-1, -1, -1,  0 };
	double d[4] = { 5,  5, 10, 23 };
	// results    { 2,  3,  5, 7  }
	tridiag_matrix_algorithm(a,b,c,d,n);
	for (int i = 0; i < n; i++) {
		cout << d[i] << endl;
	}
	cout << endl << "n= " << n << endl;
	return 0;
}