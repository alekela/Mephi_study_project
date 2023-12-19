import numpy as np
import matplotlib.pyplot as plt


def approximate_diff_equation(a, b, h, p, q):
    '''
        Функция аппроксимирует дифференциальное уравнение краевой задачи

        a, b - границы области изменения независимой переменной
        h - шаг сетки
        p - коэффициент при первой производной искомой функции в ДУ
        q - коэффициент при искомой функции в ДУ
    '''
    n = int((b - a) / h + 1)  # количество узлов в сетке
    x_grid = np.linspace(a, b, n)

    # создаем векторы коэффициентов в СЛАУ для исходного ДУ
    A = [1 / h ** 2 - p(x) / (2 * h) for x in x_grid[1:n - 1]]
    B = [-2 / h ** 2 + q(x) for x in x_grid[1:n - 1]]
    C = [1 / h ** 2 + p(x) / (2 * h) for x in x_grid[1:n - 1]]

    return A, B, C


def approximate_boundaries_single_prec(alpha, beta, h):
    '''
        Функция аппроксимирует граничные условия краевой задачи
        с первым порядком точности

        alpha - коэффициенты при искомой функции в граничных условиях
        beta - коэффициенты при производной искомой функции в граничных условиях
        h - шаг сетки
        precision - порядок точности (по умолчанию второй)
    '''

    b0 = alpha[0] - beta[0] / h
    c0 = beta[0] / h

    aN = -beta[1] / h
    bN = beta[1] / h + alpha[1]

    return b0, c0, aN, bN


def approximate_boundaries_double_prec(alpha, beta, h, p, q, a, b):
    '''
        Функция аппроксимирует граничные условия краевой задачи
        со вторым порядком точности
    '''
    b0 = -2 / h ** 2 + 2 * alpha[0] / (h * beta[0]) - p(a) * alpha[0] / beta[0] + q(a)
    c0 = 2 / h ** 2
    aN = 2 / h ** 2
    bN = -2 * alpha[1] / (h * beta[1]) - 2 / h ** 2 - p(b) * alpha[1] / beta[1] + q(b)

    return b0, c0, aN, bN


def solve_bvp(a, b, h, p, q, f, alpha, beta, gamma, precision=2):
    '''
        Функция решает краевую задачу методом конечных разностей

        a, b - границы области изменения независимой переменной
        h - шаг сетки
        p - коэффициент при первой производной искомой функции в ДУ
        q - коэффициент при искомой функции в ДУ
        alpha - коэффициенты при искомой функции в граничных условиях
        beta - коэффициенты при производной искомой функции в граничных условиях
        precision - порядок точности (по умолчанию второй)
    '''
    A, B, C = approximate_diff_equation(a, b, h, p, q)

    if precision == 1:
        b0, c0, aN, bN = approximate_boundaries_single_prec(alpha, beta, h)
        D = [gamma[0]]
        D += [f(x) for x in np.linspace(a + h, b - h, int((b - a) / h) - 1)]
        D.append(gamma[1])
    else:
        b0, c0, aN, bN = approximate_boundaries_double_prec(alpha, beta, h, p, q, a, b)
        D = [f(a) - p(a) * gamma[0] / beta[0] + 2 * gamma[0] / (beta[0] * h)]
        D += [f(x) for x in np.linspace(a + h, b - h, int((b - a) / h) - 1)]
        D.append(f(b) - 2 * gamma[1] / (beta[1] * h) - p(b) * gamma[1] / beta[1])

    A.append(aN)
    B.append(bN)
    B.insert(0, b0)
    C.insert(0, c0)

    solution = tridiag_matrix_algorithm(A, B, C, D)

    return solution


def tridiag_matrix_algorithm(a, b, c, d):
    n = len(d)
    a.insert(0, 0)
    c.append(0)
    c[0] /= b[0]
    d[0] /= b[0]

    # прямой ход прогонки
    for i in range(1, n):
        c[i] /= b[i] - a[i] * c[i - 1]
        d[i] = (d[i] - a[i] * d[i - 1]) / (b[i] - a[i] * c[i - 1])

    # обратный ход прогонки
    for i in range(n - 2, -1, -1):
        d[i] -= c[i] * d[i + 1]

    return d


# функции - коэффициенты в дифференциальном уравнение
def p(x):
    return np.tan(x)


def q(x):
    return -2 * x / np.cos(x)


def f(x):
    return 2 - 2 * x ** 3 / np.cos(x)


# теоретическое решение задачи
def u_theor(x):
    return np.sin(x) + x ** 2


def main():
    # параметры задачи
    a = 0
    b = 1
    h = 0.05
    n = int((b - a) / h + 1)

    # коэффиценты в граничных условиях
    # для левой границы
    alpha_a = 2
    beta_a = -1
    gamma_a = -1

    # для правой границы
    alpha_b = 3
    beta_b = 1
    gamma_b = 8.0647

    # решаем задачу методом конечных разностей с первым порядком точности
    sol = solve_bvp(a, b, h, p, q, f, (alpha_a, alpha_b),
                    (beta_a, beta_b), (gamma_a, gamma_b),
                    precision=1)

    # построение графика теоретического решения
    X = np.linspace(a, b, 1000)
    sol_theor = [u_theor(x) for x in X]
    plt.plot(X, sol_theor, label="theoretical")

    # построение графика численного решения
    grid = np.linspace(a, b, n)
    plt.plot(grid, sol, label="numerical")

    plt.legend()
    plt.show()


main()
