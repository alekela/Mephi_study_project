#include <iostream>
#include <cmath>
#include <vector>


using namespace std;

double phi_func(double x) {
    return cosh(x);
}

double gamma_0_func(double t) {
    return exp(t);
}

double gamma_l_func(double t) {
    return exp(1 - t);
}

double f_func(double x, double t) {
    return -sinh(x - t) - cosh(x - t);
}

double real_func(double x, double t) {
    return cosh(x - t);
}


void tridiag_matrix_algorithm(vector<double> a, vector<double> b, vector<double> c, vector<double>& d, int n) {
    c[0] /= b[0];
    d[0] /= b[0];

    for (int i = 1; i < n; i++) {
        c[i] /= b[i] - a[i]*c[i-1];
        d[i] = (d[i] - a[i]*d[i-1]) / (b[i] - a[i]*c[i-1]);
    }

    for (int i = n - 2; i >= 0; i--) {
        d[i] -= c[i]*d[i+1];
    }
}


void heat_solving_1(double time_left, double time_right, double left, double right, double** data, double a,
                    int N, int T, double* phi, double** f, double* gamma_0, double* gamma_l, 
                    double alpha_0, double beta_0, double alpha_l, double beta_l, double sigma) {
    double h = (right - left) / (N - 1);
    double tau = (time_right - time_left) / (T - 1);
    for (int j = 0; j < N; j++) {
        data[0][j] = phi[j];
    }
    vector<double> A(N), B(N), C(N);
    // fill A, B, C (they are constant on each iteration)
    for (int i = 1; i < N - 1; i++) {
        A[i] = (1 - sigma) * tau * a * a / h / h;
        B[i] = -1 - 2 * (1 - sigma) * tau * a * a / h / h;
        C[i] = (1 - sigma) * tau * a * a / h / h;
    }
    A[0] = 0;
    C[N - 1] = 0;
    B[0] = alpha_0 - beta_0 / h;
    C[0] = beta_0 / h;
    B[N - 1] = alpha_l + beta_l / h;
    A[N - 1] = -beta_l / h;
    for (int t = 1; t < T; t++) {
        vector<double> ans(N);
        for (int i = 1; i < N - 1; i++) {
            ans[i] = -data[t - 1][i] - f[t - 1][i] * tau - 
            a * a * tau * sigma / h / h * (data[t - 1][i + 1] - 2 * data[t - 1][i] + data[t - 1][i - 1]);
        }
        ans[0] = gamma_0[t];
        ans[N - 1] = gamma_l[t];
        tridiag_matrix_algorithm(A, B, C, ans, N);
        for (int i = 0; i < N; i++) {
            data[t][i] = ans[i];
        }
    }
}


void heat_solving_2(double time_left, double time_right, double left, double right, double** data, double a,
                    const int N, int T, double* phi, double** f, double* gamma_0, double* gamma_l, 
                    double alpha_0, double beta_0, double alpha_l, double beta_l, double sigma) {
    double h = (right - left) / (N - 1);
    double tau = (time_right - time_left) / (T - 1);
    for (int j = 0; j < N; j++) {
        data[0][j] = phi[j];
    }
    vector<double> A(N), B(N), C(N);
    // fill A, B, C (they are constant on each iteration)
    for (int i = 1; i < N - 1; i++) {
        A[i] = (1 - sigma) * tau * a * a / h / h;
        B[i] = -1 - 2 * (1 - sigma) * tau * a * a / h / h;
        C[i] = (1 - sigma) * tau * a * a / h / h;
    }
    A[0] = 0;
    C[N - 1] = 0;
    B[0] = alpha_0 - beta_0 * (1 / h + h / 2 / a / a / tau);
    C[0] = beta_0 / h ;
    B[N - 1] = alpha_l + beta_l * (1 / h + h / 2 / a / a / tau);
    A[N - 1] = -beta_l / h ;
    for (int t = 1; t < T; t++) {
        vector<double> ans(N);
        for (int i = 1; i < N - 1; i++) {
            ans[i] = -data[t - 1][i] - f[t - 1][i] * tau - 
            a * a * tau * sigma / h / h * (data[t - 1][i + 1] - 2 * data[t - 1][i] + data[t - 1][i - 1]);
        }
        ans[0] = gamma_0[t] - beta_0 * h / 2 / a / a / tau * (data[t - 1][0] + f[t][0] * tau);
        ans[N - 1] = gamma_l[t] + beta_l * h / 2 / a / a / tau * (data[t - 1][N - 1] + f[t][N - 1] * tau);
        tridiag_matrix_algorithm(A, B, C, ans, N);
        for (int i = 0; i < N; i++) {
            data[t][i] = ans[i];
        }
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
    double a = 1;
    double** function = (double**) malloc (sizeof (double*) * T);
    for (int i = 0; i < T; i++) {
        function[i] = (double*) malloc (sizeof (double) * N);
    }

    double *phi = (double*) malloc(sizeof(double) * N);
    for (int j = 0; j < N; j++) {
        phi[j] = phi_func(left + h * j);
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
    double alpha_0 = 1;
    double beta_0 = -1;
    double alpha_l = 1;
    double beta_l = 1;
    double sigma = 0.5;
    heat_solving_1(time_left, time_right, left, right, function, a,
    N, T, phi, f, gamma_0, gamma_l, alpha_0, beta_0, alpha_l, beta_l, sigma);

    for (int i = 0; i < T; i++) {
        cout << function[i][1] << " ";
    }
    cout << endl;

    heat_solving_2(time_left, time_right, left, right, function, 
    a, N, T, phi, f, gamma_0, gamma_l, alpha_0, beta_0, alpha_l, beta_l, sigma);
    
    for (int i = 0; i < T; i++) {
        cout << function[i][1] << " ";
    }
    cout << endl;
    for (int i = 0; i < T; i++) {
        cout << real_func(left + 1 * h, time_left + i * tau) << " ";
    }
}