#include <iostream>
#include <cmath>
#include <sstream>
#include <fstream>

using namespace std;


double func(double x) {
    return cosh(exp(-pow(x, 2)));
}


double right_der(int i, double *data, int data_size, double h) {
    if (i != data_size - 1) {
        return (data[i + 1] - data[i]) / h;
    }
    return (data[i] - data[i - 1]) / h;
}


double left_der(int i, double *data, int data_size, double h) {
    if (i != 0) {
        return (data[i] - data[i - 1]) / h;
    }
    return (data[i + 1] - data[i]) / h;
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


double second_der2(int i, double* data, int data_size, double h) {
    if (i == 0) {
        return (data[0] - 2 * data[1] + data[2]) / h / h;
    }
    else if (i == data_size - 1) {
        return -(data[data_size - 1] - 2 * data[data_size - 2] + data[data_size - 3]) / h / h;
    }
    return (data[i + 1] - 2 * data[i] + data[i - 1]) / h / h;
}


double second_der4(int i, double* data, int data_size, double h) {
    if (i == 0 || i == 1 || i == data_size - 1 || i == data_size - 2) {
        return second_der2(i, data, data_size, h);
    }
    return (16 * (data[i + 1] + data[i - 1]) - (data[i + 2] + data[i - 2]) - 30 * data[i]) / 12 / h / h;
}


int main() {
    int n = 1001;
    double left = -5;
    double right = 5;

    double h = (right - left) / (n - 1);
    double data[n];
    for (int i = 0; i < n; i++) {
        data[i] = func(left + i * h);
    }

    const char csv_file_name[64] = "work2.csv";
    std::ofstream csv_file;
    csv_file.open(csv_file_name);
    csv_file << "x,y\n";
    for (size_t i = 0; i < n; ++i){
        csv_file << left + i * h << "," << second_der4(i, data, n, h) << "\n";
    }
    csv_file.close();

}