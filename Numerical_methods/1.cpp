#include <iostream>
#include <cmath>

using namespace std;


double f(double x){
    return cos(exp(x / 3.) / 10.);
}


double Lagrange(double x, int train_grid_size, double* train_grid){
    double ans = 0;
    double tmp;
    for (int i = 0; i < train_grid_size; i++){
        tmp = 1;
        for (int j = 0; j < train_grid_size; j++){
            if (i != j){
                tmp = tmp * (x - train_grid[j]) / (train_grid[i] - train_grid[j]);
            }
        }
        ans += tmp * f(train_grid[i]);
    }
    return ans;
}


void Lagrange_equable(double left, double right, int n, double* test_values){
    double h = (double) (right - left) / (n - 1);
    double train_grid[n];
    for (int i = 0; i < n; i++){
        train_grid[i] = left + i * h;
    }
    for (int i = 0; i < n - 1; i++){
        test_values[i] = Lagrange(left + (2 * i + 1) * h / 2, n, train_grid);
    }

}


void Lagrange_Chebishev(double left, double right, int n, double* test_values){
    double train_grid[n];
    for (int i = 0; i < n; i++){
        train_grid[n - i - 1] = (left + right) / 2. + (right - left) / 2. * cos((2 * i + 1) / 2. / n * M_PI);
    }
    double test_grid[n - 1];
    for (int i = 0; i < n - 1; i++){
        test_grid[i] = (train_grid[i] + train_grid[i + 1]) / 2.;
    }
    for (int i = 0; i < n - 1; i++){
        test_values[i] = Lagrange(test_grid[i], n, train_grid);
    }
}


void Newton_train(int train_grid_size, double* train_grid, double* koeffs){
    double sep_sub;
    double tmp;
    for (int i = train_grid_size - 1; i > -1; i--){
        sep_sub = 0;
        for (int j = 0; j < i + 1; j++){
            tmp = 1;
            for (int k = 0; k < i + 1; k++){
                if (k != j){
                    tmp *= (train_grid[j] - train_grid[k]);
                }
            }
            sep_sub += f(train_grid[j]) / tmp;
        }
        koeffs[i] = sep_sub;
    }
}


double Newton(double x, int train_grid_size, double* train_grid, double* koeffs){
    double ans = 0;
    double tmp;
    double mn;
    for (int i = train_grid_size - 1; i > -1; i--){
        ans += koeffs[i];
        if (i != 0){
            ans *= (x - train_grid[i - 1]);
        }
    }
    return ans;
}


void Newton_equable(double left, double right, int n, double* test_values){
    double h = (double) (right - left) / (n - 1);
    double train_grid[n];
    for (int i = 0; i < n; i++){
        train_grid[i] = left + i * h;
    }
    double koeffs[n];
    Newton_train(n, train_grid, koeffs);
    for (int i = 0; i < n - 1; i++){
        test_values[i] = Newton(left + (2 * i + 1) * h / 2., n, train_grid, koeffs);
    }
}


void Newton_Chebishev(double left, double right, int n, double* test_values){
    double train_grid[n];
    for (int i = 0; i < n; i++){
        train_grid[n - i - 1] = (left + right) / 2. + (right - left) / 2. * cos((2 * i + 1) / 2. / n * M_PI);
    }
    double test_grid[n - 1];
    for (int i = 0; i < n - 1; i++){
        test_grid[i] = (train_grid[i] + train_grid[i + 1]) / 2.;
    }
    double koeffs[n];
    Newton_train(n, train_grid, koeffs);
    for (int i = 0; i < n - 1; i++){
        test_values[i] = Newton(test_grid[i], n, train_grid, koeffs);
    }
}


int main(){
    double left = 0;
    double right = 10;
    int n = 51;

    // Лагранж с равномерными узлами
    double test_values1[n - 1];
    Lagrange_equable(left, right, n, test_values1);

    // Лагранж с узлами Чебышева
    double test_values2[n - 1];
    Lagrange_Chebishev(left, right, n, test_values2);

    // Ньютон с равномерными узлами
    double test_values3[n - 1];
    Newton_equable(left, right, n, test_values3);
    for (int i = 0; i < n - 1; i++){
        cout << test_values3[i] << endl;
    }
    //Ньютон с узлами Чебышева
    double test_values4[n - 1];
    //Newton_Chebishev(left, right, n, test_values4);
    return 0;
}
