#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <random>
#include "windows.h"


using namespace std;


void quantization(double* data, int len_data, int levels){
    double right = INT64_MIN, left = INT64_MAX;
    for (int i = 0; i < len_data; i++) {
        if (data[i] > right) {
            right = data[i];
        }
        if (data[i] < left) {
            left = data[i];
        }
    }   

    double quants[levels];
    for (int i = 0; i < levels; i++) {
        quants[i] = left + i * (right - left) / (levels - 1); 
    }

    double up, down;
    for (int j = 0; j < len_data; j++) {
        for (int i = 1; i < levels; i++){
            if (data[j] <= quants[i]){
                down = quants[i - 1];
                up = quants[i];
                break;
            }
        }
        if (data[j] - down < up - data[j]) {
            data[j] = down;
        }
        else {
            data[j] = up;
        }
    }
}


double my_sin(double x, double period, double phi){
    return sin(x * 2 * 3.1415 / period + phi);
}


void add_impulse_noise(double P, double* data, int len_data) {
    int check = 1. / P;

    double a = rand();
    for (int i = 0; i < len_data; i++) {
        if (i % check == 0) {
            data[i] = a;
        }
    }
}


void add_Gauss_noise(double disp, double* data, int len_data) {
    double sdvig = len_data / 2;
    for (int i = 0; i < len_data; i++) {
        data[i] += 1 / (sqrt(2 * M_PI) * disp) * exp((i - sdvig) * (i - sdvig) / len_data / len_data / 2. / disp / disp);
    }
}


void make_histogram1d(double* data, int len_data, double* hist_data, int len_hist_data) {
    double right = INT64_MIN, left = INT64_MAX;
    for (int i = 0; i < len_data; i++) {
        if (data[i] > right) {
            right = data[i];
        }
        if (data[i] < left) {
            left = data[i];
        }
    }

    double quants[len_hist_data];
    for (int i = 0; i < len_hist_data; i++) {
        quants[i] = left + i * (right - left) / (len_hist_data - 1); 
    }

    for (int i = 0; i < len_data; i++) {
        for (int j = 0; j < len_hist_data; j++) {
            if (data[i] == quants[j]) {
                hist_data[j]++;
                break;
            }
        }
    }
}


void make_histogram2d(double** data, int len_data_x, int len_data_y, double* hist_data, int len_hist_data) {
    double right = INT64_MIN, left = INT64_MAX;
    for (int i = 0; i < len_data_y; i++) {
        for (int j = 0; j < len_data_x; j++) {
            if (data[i][j] > right) {
                right = data[i][j];
            }
            if (data[i][j] < left) {
                left = data[i][j];
            }
        }
    }
    double h = (right - left) / len_hist_data;
    for (int i = 0; i < len_data_y; i++) {
        for (int k = 0; k < len_data_x; k++) {
            for (int j = 0; j < len_hist_data; j++) {
                if (data[i][k] < left + (j + 1) * h) {
                    hist_data[j]++;
                    break;
                }
            }
        }
    }
}


double mean(double* hist_data, int len_hist_data, double left, double right) {
    double h = (right - left) / len_hist_data;
    double sum = 0;
    for (int i = 0; i < len_hist_data; i++) {
        sum += (left + i * h + h / 2.) * hist_data[i];
    }
    return sum;
}


double variance(double* hist_data, int len_hist_data, double left, double right) {
    double h = (right - left) / len_hist_data;
    double sum = 0;
    double me = mean(hist_data, len_hist_data, left, right);
    for (int i = 0; i < len_hist_data; i++) {
        sum += pow((left + i * h + h / 2. - me), 2) * hist_data[i];
    }
    return sum;
}


double quartile_dist(double* hist_data, int len_hist_data, double left, double right) {
    double summa = 0;
    for (int i = 0; i < len_hist_data; i++) {
        summa += hist_data[i];
    }

    double h = (right - left) / len_hist_data;
    double tmp_summa = 0;
    double left_point, right_point;
    int tmp_j = 0;
    for (int i = 0; i < len_hist_data; i++) {
        tmp_summa += hist_data[i];
        if (tmp_summa > summa / 4.) {
            left_point = left + i * h + h / 2.;
            tmp_j = i;
            break;
        }
    }
    for (int i = tmp_j + 1; i < len_hist_data; i++) {
        tmp_summa += hist_data[i];
        if (tmp_summa > 3 * summa / 4.) {
            right_point = left + i * h + h / 2.;
        }
    }
    return right_point - left_point;
}


int main() {
    double data[1000];
    double h = 0.01;
    for (int i = 0; i < 1000; i++) {
        data[i] = my_sin(i * h, 5, 0);
    }
    add_impulse_noise(0.1, data, 1000);
    const char csv_file_name[64] = "work3.csv";
    std::ofstream csv_file;
    csv_file.open(csv_file_name);
    csv_file << "time,signal\n";
    for (size_t i = 0; i < 1000; ++i){
        csv_file << (i * h) << "," << data[i] << "\n";
    }
    csv_file.close();
}

