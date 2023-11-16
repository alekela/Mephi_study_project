#include <iostream>
#include <cmath>

using namespace std;


double right_side_f(double x, double u, double du) {
    return (3. - 2. * x + 4. * x * x) / (1 + x * x) * exp(-2. * x) + u / (1 + x * x) - du * x / (1 + x * x);
    //return du;
}


void Euler(double* x_grid, double* f, double* df, int len, double h) {
    for (int i = 1; i < len; i++) {
        df[i] = right_side_f(x_grid[i - 1], f[i - 1], df[i - 1]) * h + df[i - 1];
        f[i] = df[i - 1] * h + f[i - 1];
    }
}


void Runge_Kutta4(double* x_grid, double* f, double* df, int len, double h) {
    double q0, q1, q2, q3;
    double k0, k1, k2, k3;

    for (int i = 1; i < len; i++) {
    
        q0 = right_side_f(x_grid[i - 1], f[i - 1], df[i - 1]);
        k0 = df[i - 1];
        q1 = right_side_f(x_grid[i - 1] + h / 2., f[i - 1] + k0 * h / 2., df[i - 1] + q0 * h / 2.);
        k1 = df[i - 1] + q0 * h / 2.;
        q2 = right_side_f(x_grid[i - 1] + h / 2., f[i - 1] + k1 * h / 2., df[i - 1] + q1 * h / 2.);
        k2 = df[i - 1] + q1 * h / 2.;
        q3 = right_side_f(x_grid[i - 1] + h, f[i - 1] + k2 * h, df[i - 1] + q2 * h);
        k3 = df[i - 1] + q2 * h;

        df[i] = df[i - 1] + h / 6. * (q0 + 2 * q1 + 2 * q2 + q3);
        f[i] = f[i - 1] + h / 6. * (k0 + 2 * k1 + 2 * k2 + k3);
    }
}


void Adams(double* x_grid, double* f, double* df, int len, double h) {
    double k1, k2, k3;

    for (int i = 3; i < len; i++) {
    
        k1 = right_side_f(x_grid[i - 1], f[i - 1], df[i - 1]);
        k2 = right_side_f(x_grid[i - 2], f[i - 2], df[i - 2]);
        k3 = right_side_f(x_grid[i - 3], f[i - 3], df[i - 3]);

        df[i] = df[i - 1] + (23. * k1 - 16. * k2 + 5. * k3) * h / 12.;
        f[i] = f[i - 1] + (23. * df[i - 1] - 16. * df[i - 2] + 5. * df[i - 3]) * h / 12.;
    }
}


int main() {
    // инициализация первой сетки
    const int n1 = 101;
    double x_grid1[n1];
    double h1 = 0.05;
    for (int i = 0; i < n1; i++) {
        x_grid1[i] = i * h1;
    }
    double f1[n1];
    double df1[n1];

    // начальные условия
    f1[0] = 2;
    df1[0] = -2;
    f1[1] = 1.90609;
    f1[2] = 1.82372;
    df1[1] = -1.75974;
    df1[2] = -1.53796;

    // инициализация второй сетки
    const int n2 = n1 / 2 + 1;
    double x_grid2[n2];
    double h2 = 0.1;
    for (int i = 0; i < n2; i++) {
        x_grid1[i] = i * h2;
    }
    double f2[n2];
    double df2[n2];

    // начальные условия
    f2[0] = 2;
    df2[0] = -2;
    f2[1] = 1.90609;
    f2[2] = 1.82372;
    df2[1] = -1.75974;
    df2[2] = -1.53796;

    // запуск метода
    Runge_Kutta4(x_grid1, f1, df1, n1, h1);
    Runge_Kutta4(x_grid2, f2, df2, n2, h2);

    // оценка погрешностей
    double max_error = 0;
    double tmp;
    for (int i = 0; i < n1; i++ ) {
        if (i % 2 == 0) {
            tmp = abs(f1[i] - f2[i / 2]);
        }
        else {
            tmp = abs(f1[i] - (f2[i / 2] + f2[i / 2 + 1]) / 2.);
        }
        if (tmp > max_error) {
            max_error = tmp;
        }
    }
    max_error /= 15.; // для Рунге-Кутта
    // max_error /= 7.; // для Адамса
    cout << "Максимальный модуль ошибки равен " << max_error << endl;
}