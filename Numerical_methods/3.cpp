#include <iostream>
#include <cmath>

using namespace std;


double func(double x) {
    return sqrt(2 - x * x);
}


double rect(int n, double left, double right) {
    double h = (right - left) / n;
    double sum = 0;
    for (int i = 0; i < n; i++) {
        sum += func(left + i * h + h / 2.);
    }
    return sum * h;
}


double trapezoid(int n, double left, double right) {
    double h = (right - left) / n;
    double sum = 0;
    double a = left;
    for (int i = 0; i < n; i++) {
        sum += (func(left + h * (i + 1)) + func(left + i * h));
    }
    return sum * h / 2;
}


double simpson(int n, double left, double right) {
    double h = (right - left) / n;
    double sum = 0;
    for(int i = 0; i < n - 1; i++) {
        sum += (func(left + i * h) + 4.0*func(left + (i + 1) * h) + func(left + (i + 2) * h));
    }

    return sum * h / 6.;
}


int main() {
    double left, right;
    left = 0;
    right = 1;
    int n;
    cout << "Введите количество шагов:" << endl;
    cin >> n;
    cout << rect(n, left, right) << endl;
    cout << trapezoid(n, left, right) << endl;
    cout << simpson(n, left, right) << endl;

}