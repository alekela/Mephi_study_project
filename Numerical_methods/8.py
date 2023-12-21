from math import *
import matplotlib.pyplot as plt
import numpy as np


def phi_func(x):
    return cosh(x)


def gamma_0_func(t):
    return exp(t)


def gamma_l_func(t):
    return exp(1 - t)


def f_func(x, t):
    return -sinh(x - t) - cosh(x - t)


def real_func(x, t):
    return cosh(x - t)


def tridiag_matrix_algorithm(a, b, c, d, n):
    c[0] /= b[0]
    d[0] /= b[0]

    for i in range(1, n):
        c[i] /= b[i] - a[i] * c[i - 1]
        d[i] = (d[i] - a[i] * d[i - 1]) / (b[i] - a[i] * c[i - 1])

    for i in range(n - 2, -1, -1):
        d[i] -= c[i] * d[i + 1]


def heat_solving_1(time_left, time_right, left, right, data, a, N,
                   T, phi, f, gamma_0, gamma_l, alpha_0, beta_0, alpha_l, beta_l, sigma):
    h = (right - left) / (N - 1)
    tau = (time_right - time_left) / (T - 1)
    for j in range(N):
        data[0][j] = phi[j]

    A = [0] * N
    B = [0] * N
    for i in range(1, N - 1):
        A[i] = (1 - sigma) * tau * a * a / h / h
        B[i] = -1 - 2 * (1 - sigma) * tau * a * a / h / h
    A[0] = 0
    B[0] = alpha_0 - beta_0 / h
    B[N - 1] = alpha_l + beta_l / h
    A[N - 1] = -beta_l / h
    for t in range(1, T):
        C = [0] * N
        C[N - 1] = 0
        C[0] = beta_0 / h
        ans = [0] * N
        for i in range(1, N - 1):
            C[i] = (1 - sigma) * tau * a * a / h / h
            ans[i] = -data[t - 1][i] - f[t - 1][i] * tau - \
                     a * a * tau * sigma / h / h * (data[t - 1][i + 1] - 2 * data[t - 1][i] + data[t - 1][i - 1])

        ans[0] = gamma_0[t]
        ans[N - 1] = gamma_l[t]
        tridiag_matrix_algorithm(A, B, C, ans, N)
        for i in range(N):
            data[t][i] = ans[i]


def heat_solving_2(time_left, time_right, left, right, data, a, N,
                   T, phi, f, gamma_0, gamma_l, alpha_0, beta_0, alpha_l, beta_l, sigma):
    h = (right - left) / (N - 1)
    tau = (time_right - time_left) / (T - 1)
    for j in range(N):
        data[0][j] = phi[j]

    A = [0] * N
    B = [0] * N
    for i in range(1, N - 1):
        A[i] = (1 - sigma) * tau * a * a / h / h
        B[i] = -1 - 2 * (1 - sigma) * tau * a * a / h / h
    A[0] = 0
    B[0] = alpha_0 - beta_0 * (1 / h + h / 2 / a / a / tau)
    B[N - 1] = alpha_l + beta_l * (1 / h + h / 2 / a / a / tau)
    A[N - 1] = -beta_l / h
    for t in range(1, T):
        C = [0] * N
        C[N - 1] = 0
        C[0] = beta_0 / h
        ans = [0] * N
        for i in range(1, N - 1):
            C[i] = (1 - sigma) * tau * a * a / h / h
            ans[i] = -data[t - 1][i] - f[t - 1][i] * tau - \
                     a * a * tau * sigma / h / h * (data[t - 1][i + 1] - 2 * data[t - 1][i] + data[t - 1][i - 1])

        ans[0] = gamma_0[t] - beta_0 * h / 2 / a / a / tau * (data[t - 1][0] + f[t][0] * tau)
        ans[N - 1] = gamma_l[t] + beta_l * h / 2 / a / a / tau * (data[t - 1][N - 1] + f[t][N - 1] * tau)
        tridiag_matrix_algorithm(A, B, C, ans, N)
        for i in range(N):
            data[t][i] = ans[i]


time_left = 0
time_right = 1
left = 0
right = 1
N = 21
T = 21
h = (right - left) / (N - 1)
tau = (time_right - time_left) / (T - 1)
a = 1
function = [[0 for _ in range(N)] for _ in range(T)]

phi = [0] * N
for j in range(N):
    phi[j] = phi_func(left + h * j)

f = [[0 for _ in range(N)] for _ in range(T)]
for i in range(T):
    for j in range(N):
        f[i][j] = f_func(left + h * j, time_left + tau * i)

gamma_0 = [0] * T
gamma_l = [0] * T
for i in range(T):
    gamma_0[i] = gamma_0_func(time_left + i * tau)
    gamma_l[i] = gamma_l_func(time_left + i * tau)

alpha_0 = 1
beta_0 = -1
alpha_l = 1
beta_l = 1
sigma = 0.5

fig, ax = plt.subplots(subplot_kw=dict(projection='3d'))
x_grid = [left + i * h for i in range(N)]
y_grid = [time_left + i * tau for i in range(N)]
x_grid, y_grid = np.meshgrid(x_grid, y_grid)

heat_solving_1(time_left, time_right, left, right, function,
               a, N, T, phi, f, gamma_0, gamma_l, alpha_0, beta_0, alpha_l, beta_l, sigma)
for i in range(T):
    print(function[i][1], end=" ")
print()
function = np.array(function)
ax.plot_surface(x_grid, y_grid, function)


heat_solving_2(time_left, time_right, left, right, function,
               a, N, T, phi, f, gamma_0, gamma_l, alpha_0, beta_0, alpha_l, beta_l, sigma)
for i in range(T):
    print(function[i][1], end=" ")
print()
for i in range(T):
    print(real_func(left + 1 * h, time_left + i * tau), end=" ")
function = np.array(function)
ax.plot_surface(x_grid, y_grid, function)

ax.plot_surface(x_grid, y_grid,
                np.array([[real_func(left + j * h, time_left + i * tau) for j in range(N)] for i in range(T)]))

plt.show()
