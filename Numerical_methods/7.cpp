#include <iostream>
#include <cmath>

using namespace std;

double phi_func(double x) {
    return x - x * tan(x / 2.);
}

double psi_func(double x) {
    return 1 + x / (1 + cos(x));
}

double gamma_0_func(double t) {
    return 1 + tan(t / 2.);
}

double gamma_l_func(double t) {
    return t + 1 / (1 + cos(t - 1));
}

double f_func(double x, double t) {
    return (2 + 2 * cos(t - x) + x * sin(t - x)) / 2. / (1 + cos(t - x)) / (1 + cos(t - x));
}

double real_func(double x, double t) {
    return t + x + x * tan((t - x) / 2.);
}


void wave_solving_1(double time_left, double time_right, double left, double right, double** data, double a,
                    int N, int T, double* phi, double *psi, double** f, double* gamma_0, double* gamma_l, 
                    double alpha_0, double beta_0, double alpha_l, double beta_l) {
    double h = (right - left) / (N - 1);
    double tau = (time_right - time_left) / (T - 1);
    for (int j = 0; j < N; j++) {
        data[0][j] = phi[j];
        data[1][j] = phi[j] + tau * psi[j];
    }
    for (int i = 2; i < T; i++) {
        for (int j = 1; j < N - 1; j++) {
            data[i][j] = a * a * tau * tau / h / h * 
            (data[i - 1][j + 1] - 2 * data[i - 1][j] + data[i - 1][j - 1]) -
            data[i - 2][j] + 2 * data[i - 1][j] + tau * tau * f[i][j];
        }
        data[i][0] = (gamma_0[i] * h - beta_0 * data[i][1]) / (alpha_0 * h - beta_0);
        data[i][N - 1] = (gamma_l[i] * h + beta_l * data[i][N - 2]) / (alpha_l * h + beta_l);
    }
}


void wave_solving_2(double time_left, double time_right, double left, double right, double** data, double a,
                    int N, int T, double* phi, double *psi, double** f, double* gamma_0, double* gamma_l, 
                    double alpha_0, double beta_0, double alpha_l, double beta_l) {
    double h = (right - left) / (N - 1);
    double tau = (time_right - time_left) / (T - 1);
    for (int j = 0; j < N; j++) {
        data[0][j] = phi[j];
    }
    for (int j = 0; j < N; j++) {
        data[1][j] = (2. * tau * psi[j] + 2. * phi[j] + 
        a * a * tau * tau / h / h * (data[0][j + 1] - 2 * phi[j] + data[0][j - 1]) + 
        tau * tau * f[0][j]) / 2.;
    }

    for (int i = 2; i < T; i++) {
        for (int j = 1; j < N - 1; j++) {
            data[i][j] = a * a * tau * tau / h / h * 
            (data[i - 1][j + 1] - 2 * data[i - 1][j] + data[i - 1][j - 1]) -
            data[i - 2][j] + 2 * data[i - 1][j] + tau * tau * f[i][j];
        }
        data[i][0] = (gamma_0[i] - beta_0 * (4 * data[i][1] - data[i][2]) / 2. / h) / (alpha_0 - beta_0 * 3 / 2. / h);
        data[i][N - 1] = (gamma_l[i] - beta_l * (data[i][N - 3] - 4 * data[i][N - 2]) / 2. / h) / (alpha_l + 3 * beta_l / 2. / h);
    }

}


int main() {
    double time_left = 0;
    double time_right = 1;
    double left = 0;
    double right = 1;
    int N = 21;
    int T = 21;
    double h = (right - left) / (N - 1);
    double tau = (time_right - time_left) / (T - 1);
    double** function = (double**) malloc (sizeof (double*) * T);
    for (int i = 0; i < T; i++) {
        function[i] = (double*) malloc (sizeof (double) * N);
    }
    double a = 0.5;

    double *phi = (double*) malloc(sizeof(double) * N);
    double *psi = (double*) malloc(sizeof(double) * N);
    for (int j = 0; j < N; j++) {
        phi[j] = phi_func(left + h * j);
        psi[j] = psi_func(left + h * j);
    }

    double** f = (double**) malloc (sizeof (double*) * T);
    for (int i = 0; i < T; i++) {
        f[i] = (double*) malloc (sizeof (double) * N);
    }
    for (int i = 0; i < T; i++) {
        for (int j = 0; j < N; j++) {
            f[i][j] = f_func(left + h * j, time_left + tau * i);
        }
    }

    double* gamma_0 = (double*) malloc(sizeof(double) * T);
    double* gamma_l = (double*) malloc(sizeof(double) * T);
    for (int i = 0; i < T; i++) {
        gamma_0[i] = gamma_0_func(time_left + i * tau);
        gamma_l[i] = gamma_l_func(time_left + i * tau);
    }
    double alpha_0 = 0;
    double beta_0 = 1;
    double alpha_l = 1;
    double beta_l = -1;
    wave_solving_2(time_left, time_right, left, right, function, a, N, T, phi, psi, f, gamma_0, gamma_l, alpha_0, beta_0, alpha_l, beta_l);

    for (int i = 0; i < T; i++) {
        cout << function[i][1] << " ";
    }
    cout << endl;

    wave_solving_1(time_left, time_right, left, right, function, a, N, T, phi, psi, f, gamma_0, gamma_l, alpha_0, beta_0, alpha_l, beta_l);
    
    for (int i = 0; i < T; i++) {
        cout << function[i][1] << " ";
    }
    cout << endl;
    for (int i = 0; i < T; i++) {
        cout << real_func(left + 1 * h, time_left + i * tau) << " ";
    }
}