#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

using namespace std;


// функция, нули которой необходимо найти
long double f(long double x) {
    return exp(sin(x/2)) - atan(x) + 1;
}


// производная
long double f_deriv(long double x) {
    return exp(sin(x/2)) * cos(x/2) / 2 - 1/(pow(x,2) + 1);
}


// численное нахождение производной функции
long double centr_der(int i, long double *data, int data_size, long double h) {
    if (i == 0) {
        return (4 * data[1] - data[2] - 3 * data[0]) / 2 / h;
    }
    else if (i == data_size - 1) {
        return -(4 * data[data_size - 2] - data[data_size - 3] - 3 * data[data_size - 1]) / 2 / h;
    }
    return (data[i + 1] - data[i - 1]) / 2 / h;
}


// метод бисекции
long double bisection(long double left, long double right, long double eps) {
    long double x = (right + left) / 2., prev_x = right;
    int n = log2((right - left) / eps) + 1; 
    /* у данного метода можно рассчитать количество итераций 
    и реализовать его циклом for, без while и prev_x, получится точнее*/
    while (abs(prev_x - x) > eps){
        prev_x = x;
        if (f(x) * f(right) < 0)
            left = x;
        else
            right = x;
        x = (right + left) / 2.;
    }
    return x;
}


// метод Ньютона
long double newton_modif(long double f( long double x),
                long double f0, 
                long double x0, 
                long double eps) 
{
    /*
        f - левая часть решаемого уравнения
        f0 - производная левой части в точке x0
        x0 - начальное приближение
        eps - точность решения
    */

    long double x1 = x0 - f(x0)/f0;
    bool crit1 = abs(x1 - x0) > eps;        // флаг для отслеживания первого критерия останова
    bool crit2 = abs(f(x1) - f(x0)) > eps;  // флаг для отслеживания второго критерия останова

    while (crit1) {
        x0 = x1;
        x1 = x1 - f(x1)/f0;    // считаем новые приближенные значения корня
        crit1 = abs(x1 - x0) > eps;    // проверяем выполнение критериев останова
        crit2 = abs(f(x1) - f(x0)) > eps;
    }
    return x1;
}


// метод хорд
long double chords(long double left, long double right, long double eps) {
    long double x, prev_x = left;
    x = -f(left) / (f(right) - f(left)) * (right - left);
    while (abs(x - prev_x) > eps) {
        prev_x = x;
        if (f(x) * f(left) < 0) right = x;
        else left = x;
        x = -f(left) / (f(right) - f(left)) * (right - left);
    }
    return x;
}


int main() {
    long double left0 = 0, right0 = 10;
    const int n0 = 10000;
    long double h = (right0 - left0) / (n0 - 1);

    long double data[n0];
    for (int i = 0; i < n0; i++) {
        data[i] = f(left0 + i * h);
    }

    long double der_data[n0];
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
    vector <long double> clear_starts, clear_ends;
    for (int i = 0; i < starts.size(); i++) {
        if (data[starts[i]] * data[ends[i]] < 0) {
            clear_starts.push_back(left0 + starts[i] * h);
            clear_ends.push_back(left0 + ends[i] * h);
        }
    }

    long double eps_s[3] = {1.0E-3, 1.0E-6, 1.0E-9};
    long double eps;
    for (int j = 0; j < 3; j++) {
        eps = eps_s[j];
        cout << "Для точности " << eps << ":" << endl;

        cout << "Метод бисекции:" << endl;
        for (int i = 0; i < clear_starts.size(); i++) {
            long double left = clear_starts[i];
            long double right = clear_ends[i];
            long double x;

            x = bisection(left, right, eps);
            cout << "Корень номер " << i + 1 << ": " << setprecision(j * 3 + 4) << x << endl;
            // |x - previous_x| < eps = True, так как является критерием останова для метода
            cout << "|x - previous_x| < eps = True, " << "|f(x)| < eps = ";
            if (abs(f(x)) < eps) cout << "True" << endl;
            else cout << "False" << endl;
            cout << endl;
        }

        cout << "Модифицированный метод Ньютона:" << endl;
        for (int i = 0; i< clear_starts.size(); i++) {
            long double left = clear_starts[i];
            long double right = clear_ends[i];
            long double x;
            long double x0 = (right + left) / 2.;

            x = newton_modif(f, f_deriv(x0), x0, eps);
            cout << "Корень номер " << i + 1 << ": " << setprecision(j * 3 + 4) << x << endl;
            // |x - previous_x| < eps = True, так как является критерием останова для метода
            cout << "|x - previous_x| < eps = True, " << "|f(x)| < eps = ";
            if (abs(f(x)) < eps) cout << "True" << endl;
            else cout << "False" << endl;
            cout << endl;
        }

        cout << "Метод хорд:" << endl;
        for (int i = 0; i< clear_starts.size(); i++) {
            long double left = clear_starts[i];
            long double right = clear_ends[i];
            long double x;

            x = chords(left, right, eps);
            cout << "Корень номер " << i + 1 << ": " << setprecision(j * 3 + 4) << x << endl;
            // |x - previous_x| < eps = True, так как является критерием останова для метода
            cout << "|x - previous_x| < eps = True, " << "|f(x)| < eps = ";
            if (abs(f(x)) < eps) cout << "True" << endl;
            else cout << "False" << endl;
            cout << "\n\n";
        }
    }
}