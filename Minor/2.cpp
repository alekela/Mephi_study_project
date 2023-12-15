#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <complex.h>

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


void ComplexBitReverse(complex<double>* data, int size) {
    int middle = size / 2,
    revSize = size - 1,
    j = 0;
    for (int i = 0; i < revSize; ++i) {
        if(i < j) {
            swap(data[i], data[j]);
        }
        int k = middle;
        while (k <= j) {
            j -= k;
            k /= 2;
        }
        j += k;
    }
}


void FftDit( complex<double>* data,  int size, int sizeLog2, int dir ) {
    ComplexBitReverse(data, size);  // переставить в бит-реверсивном порядке
    int ptsInLeftDft,ptsInRightDft = 1;
    for ( int stage = 1; stage <= sizeLog2; ++stage )
    {
        ptsInLeftDft = ptsInRightDft;   // установить ptsInLeftDFT = 2**(stage-1)
        ptsInRightDft *= 2;             // установить ptsInRightDFT = 2**stage
        complex<double> twiddle = complex<double>(1.0, 0.0); // поворачивающий множ.
        double trigArg = M_PI / ptsInLeftDft;  
        // dir == 1 для прямого преобразования, dir == -1 для обратного
        complex<double> wFactor = complex<double>(cos(trigArg),-sin(trigArg)*dir);
        for( int butterflyPos = 0; butterflyPos < ptsInLeftDft; ++butterflyPos )
        {                             
            for(int topNode=butterflyPos; topNode < size; topNode+=ptsInRightDft )
            {                              
                int botNode = topNode + ptsInLeftDft;
                complex<double> temp = data[botNode] * twiddle;
                data[botNode] = data[topNode] - temp;
                data[topNode] += temp;
            }  // конец цикла по topNode

            twiddle *= wFactor;
        } // конец цикла "бабочка"
    } // конец цикла stage
}



int main() {
    double F = 16.;
    double N = 10;
    const int N_samples = pow(2, ceil(log2(N * F)));
    double signal_value;
    double phi = 3.14 / 2.;
    double period = 10;
    const char csv_file_name[64] = "work2.csv";
    complex<double> data[N_samples];
    // для случая двух случайных прямоугольников
    /*double t = 1.;
    double right_point = N - 2 * t;
    double left_point = 0;
    double first = rand() / (double) RAND_MAX * (right_point - left_point);
    right_point = N - t;
    left_point = first + t;
    double second = rand() / (double) RAND_MAX * (right_point - left_point);*/
    for (size_t i = 0; i < N_samples; ++i) {
        if (i > N * F) {
            data[i] = 0;
        }
        else {
            signal_value = my_sin(i / F, 3, 0);
            signal_value = quantization(signal_value, 256);
            data[i] = signal_value;
        }
    }
    FftDit(data, N_samples, (int) log2(N_samples), 1);

    for (int i = 0; i < N_samples / 2; i++) {
        swap(data[i], data[N_samples / 2 + i]);
    }

    std::ofstream csv_file;
    csv_file.open(csv_file_name);
    csv_file << "time,signal\n";
    for (size_t i = 0; i < N_samples; i++) {
        csv_file << (i / F) << "," << abs(data[i]) << "\n";
    }
    csv_file.close();
    return 0;
}
