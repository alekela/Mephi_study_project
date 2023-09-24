#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>

using namespace std;


double quantization(double value, int levels){
    bool negative = false;
    if (value < 0){
        value = -value;
        negative = true;
    }
    double up, down;
    for (int i = 1; i < levels + 1; i++){
        if (value <= (double) i / (double) levels){
            down = (double) (i - 1) / (double) levels;
            up = (double) i / (double) levels;
            break;
        }
    }
    double ans;
    if (value - down < up - value){
        ans = down;
    }
    else{
        ans = up;
    }
    if (negative){
        ans = -ans;
    }
    return ans;
}


double my_sin(double x, double period, double phi){
    return sin(x * 2 * 3.1415 / period + phi);
}


double triangle(double x, double period, double t, double phi){
    x = x + phi;
    while (x > 0){
        x = x - period;
    }
    x = x + period;
    if (x > t){
        return 0;
    }
    else if (x < t / 2.){
        return 2 * x / t;
    }
    return 2 - 2 * x / t;
}

double rectangular(double x, double period, double t, double phi){
    x = x + phi;
    while (x > 0){
        x = x - period;
    }
    x = x + period;
    if (x > t){
        return 0;
    }
    return 1;
}


double two_rects(double x, double t, double first, double second){
    if ((x > first && x < first + t)){
        return 1;
    }
    if ((x > second && x < second + t)){
        return 1;
    }
    return 0;
}


int main(){
    // без квантования значений функции
    double F = 10.;
    double N = 5;
    int N_samples = (int) N * F;
    double signal_value;
    double phi = 3.14 / 2.;
    double period = 10;
    const char csv_file_name[64] = "work1.csv";
    std::ofstream csv_file;
    csv_file.open(csv_file_name);
    csv_file << "time,signal\n";
    // для случая двух случайных прямоугольников
    /*double t = 1.;
    double right_point = N - 2 * t;
    double left_point = 0;
    double first = rand() / (double) RAND_MAX * (right_point - left_point);
    right_point = N - t;
    left_point = first + t;
    double second = rand() / (double) RAND_MAX * (right_point - left_point);*/
    for (size_t i = 0; i < N_samples; ++i){
        signal_value = my_sin(i / F, 3, 0);
        signal_value = quantization(signal_value, 256);
        csv_file << (i / F) << "," << signal_value << "\n";
    }
    csv_file.close();

}