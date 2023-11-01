#include <iostream>
#include <vector>
#include <cmath>

using namespace std;


double func(double x) {
    return sin(x);
}


double centr_der(int i, double *data, int data_size, double h) {
    if (i == 0) {
        return (4 * data[1] - data[2] - 3 * data[0]) / 2 / h;
    }
    else if (i == data_size - 1) {
        return -(4 * data[data_size - 2] - data[data_size - 3] - 3 * data[data_size - 1]) / 2 / h;
    }
    return (data[i + 1] - data[i - 1]) / 2 / h;
}


double bisection(double left, double right, double eps) {
    double x;
    int n = log2((right - left) / eps) + 1;
    for (int i = 0; i < n; i++) {
        x = (right + left) / 2.;
        if (func(x) * func(right) < 0) {
            left = x;
        }
        else {
            right = x;
        }
    }
    return x;
}


int main() {
    double left0 = 1, right0 = 8;
    const int n0 = 8;
    double h = (right0 - left0) / (n0 - 1);

    double data[n0];
    for (int i = 0; i < n0; i++) {
        data[i] = func(left0 + i * h);
    }

    double der_data[n0];
    for (int i = 0; i < n0; i++) {
        der_data[i] = centr_der(i, data, n0, h);
    }
    vector<int> starts;
    vector<int> ends;
    starts.push_back(0);
    for (int i = 1; i < n0; i++) {
        if (der_data[i - 1] * der_data[i] < 0) {
            starts.push_back(i);
            ends.push_back(i -1);
        }
    }
    ends.push_back(n0 - 1);
    vector <double> clear_starts, clear_ends;
    for (int i = 0; i < starts.size(); i++) {
        if (data[starts[i]] * data[ends[i]] < 0) {
            clear_starts.push_back(left0 + starts[i] * h);
            clear_ends.push_back(left0 + ends[i] * h);
        }
    }

    double eps;
    cout << "Введите точность:" << endl;
    cin >> eps;

    for (int i = 0; i < clear_starts.size(); i++) {
        double left = clear_starts[i];
        double right = clear_ends[i];

        double x = bisection(left, right, eps);
        cout << "Найденный корень из промежутка от " << clear_starts[i] << " до " << clear_ends[i] << " равен " << x << endl;
    }

}