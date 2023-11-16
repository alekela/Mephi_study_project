from math import exp
import numpy as np
import matplotlib.pyplot as plt


def right_side_f(x, u, du):
    return (3. - 2. * x + 4. * x * x) / (1 + x * x) * exp(-2. * x) + u / (1 + x * x) - du * x / (1 + x * x)
    # return du


def Euler(x_grid, f, df, length, h):
    for i in range(1, length):
        df[i] = right_side_f(x_grid[i - 1], f[i - 1], df[i - 1]) * h + df[i - 1]
        f[i] = df[i - 1] * h + f[i - 1]


def Runge_Kutta4(x_grid, f, df, length, h):
    for i in range(1, length):
        q0 = right_side_f(x_grid[i - 1], f[i - 1], df[i - 1])
        k0 = df[i - 1]
        q1 = right_side_f(x_grid[i - 1] + h / 2., f[i - 1] + k0 * h / 2., df[i - 1] + q0 * h / 2.)
        k1 = df[i - 1] + q0 * h / 2.
        q2 = right_side_f(x_grid[i - 1] + h / 2., f[i - 1] + k1 * h / 2., df[i - 1] + q1 * h / 2.)
        k2 = df[i - 1] + q1 * h / 2.
        q3 = right_side_f(x_grid[i - 1] + h, f[i - 1] + k2 * h, df[i - 1] + q2 * h)
        k3 = df[i - 1] + q2 * h

        df[i] = df[i - 1] + h / 6. * (q0 + 2 * q1 + 2 * q2 + q3)
        f[i] = f[i - 1] + h / 6. * (k0 + 2 * k1 + 2 * k2 + k3)


def Adams(x_grid, f, df, length, h):
    for i in range(3, length):
        k1 = right_side_f(x_grid[i - 1], f[i - 1], df[i - 1])
        k2 = right_side_f(x_grid[i - 2], f[i - 2], df[i - 2])
        k3 = right_side_f(x_grid[i - 3], f[i - 3], df[i - 3])

        df[i] = df[i - 1] + (23. * k1 - 16. * k2 + 5. * k3) * h / 12.
        f[i] = f[i - 1] + (23. * df[i - 1] - 16. * df[i - 2] + 5. * df[i - 3]) * h / 12.


for q in range(3):
    # инициализация первой сетки
    n1 = 101
    h1 = 0.05
    x_grid1 = [i * h1 for i in range(n1)]
    f1 = [0] * n1
    df1 = [0] * n1

    # начальные условия
    f1[0] = 2
    df1[0] = -2
    f1[1] = 1.90609
    f1[2] = 1.82372
    df1[1] = -1.75974
    df1[2] = -1.53796

    # инициализация второй сетки
    n2 = n1 // 2 + 1
    h2 = 0.1
    x_grid2 = [i * h2 for i in range(n2)]
    f2 = [0] * n2
    df2 = [0] * n2

    # начальные условия
    f2[0] = 2
    df2[0] = -2
    f2[1] = 1.90609
    f2[2] = 1.82372
    df2[1] = -1.75974
    df2[2] = -1.53796

    # запуск метода
    if q == 0:
        Euler(x_grid1, f1, df1, n1, h1)
        Euler(x_grid2, f2, df2, n2, h2)
    elif q == 1:
        Runge_Kutta4(x_grid1, f1, df1, n1, h1)
        Runge_Kutta4(x_grid2, f2, df2, n2, h2)
    elif q == 2:
        Adams(x_grid1, f1, df1, n1, h1)
        Adams(x_grid2, f2, df2, n2, h2)

    # оценка погрешностей
    errors = []
    for i in range(n1):
        if i % 2 == 0:
            tmp = abs(f1[i] - f2[i // 2])
        else:
            tmp = abs(f1[i] - (f2[i // 2] + f2[i // 2 + 1]) / 2.)
        errors.append(tmp)

    errors = np.array(errors)
    if q == 1:
        errors /= 15.  # для Рунге-Кутта
    elif q == 2:
        errors /= 7.  # для Адамса

    plt.plot(x_grid1, errors)

plt.legend(["Эйлер", "Рунге-Кутта 4 порядка", "Адамс"])
plt.show()
