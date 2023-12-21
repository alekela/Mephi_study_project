from math import *
import matplotlib.pyplot as plt
import numpy as np


def phi_func(x):
    return x - x * tan(x / 2.)


def psi_func(x):
    return 1 + x / (1 + cos(x))


def gamma_0_func(t):
    return 1 + tan(t / 2.)


def gamma_l_func(t):
    return t + 1 / (1 + cos(t - 1))


def f_func(x, t):
    return (2 + 2 * cos(t - x) + x * sin(t - x)) / 2. / (1 + cos(t - x)) / (1 + cos(t - x))


def real_func(x, t):
    return t + x + x * tan((t - x) / 2.)


def wave_solving_1(time_left, time_right, left, right, data, a, N,
                   T, phi, psi, f, gamma_0, gamma_l, alpha_0, beta_0, alpha_l, beta_l):
    h = (right - left) / (N - 1)

    tau = (time_right - time_left) / (T - 1)
    for j in range(N):
        data[0][j] = phi[j]
        data[1][j] = phi[j] + tau * psi[j]

    for i in range(2, T):
        for j in range(1, N - 1):
            data[i][j] = a * a * tau * tau / h / h * \
                         (data[i - 1][j + 1] - 2 * data[i - 1][j] + data[i - 1][j - 1]) - \
                         data[i - 2][j] + 2 * data[i - 1][j] + tau * tau * f[i][j]

        data[i][0] = (gamma_0[i] * h - beta_0 * data[i][1]) / (alpha_0 * h - beta_0)
        data[i][N - 1] = (gamma_l[i] * h + beta_l * data[i][N - 2]) / (alpha_l * h + beta_l)


def wave_solving_2(time_left, time_right, left, right, data,
                   a, N, T, phi, psi, f, gamma_0, gamma_l, alpha_0, beta_0, alpha_l, beta_l):
    h = (right - left) / (N - 1)

    tau = (time_right - time_left) / (T - 1)
    for j in range(N):
        data[0][j] = phi[j]

    for j in range(1, N - 1):
        data[1][j] = 0.5 * (2. * tau * psi[j] + 2. * phi[j] +
                            a * a * tau * tau / h / h * (data[0][j + 1] - 2 * phi[j] + data[0][j - 1]) +
                            tau * tau * f[0][j])
    data[1][0] = 0.5 * (2. * tau * psi[0] + 2. * phi[0] +
                        a * a * tau * tau / h / h * (2 * data[0][0] - 5 * data[0][1] + 4 * data[0][2] - data[0][3]) +
                        tau * tau * f[0][0])
    data[1][N - 1] = 0.5 * (2. * tau * psi[N - 1] + 2. * phi[N - 1] -
                            a * a * tau * tau / h / h * (
                                        2 * data[0][N - 1] - 5 * data[0][N - 2] + 4 * data[0][N - 3] - data[0][N - 4]) +
                            tau * tau * f[0][N - 1])
    for i in range(2, T):
        for j in range(1, N - 1):
            data[i][j] = a * a * tau * tau / h / h * \
                         (data[i - 1][j + 1] - 2 * data[i - 1][j] + data[i - 1][j - 1]) - \
                         data[i - 2][j] + 2 * data[i - 1][j] + tau * tau * f[i][j]

        data[i][0] = (gamma_0[i] - beta_0 * (4 * data[i][1] - data[i][2]) / 2. / h) / (
                    alpha_0 - beta_0 * 3 / 2. / h)
        data[i][N - 1] = (gamma_l[i] - beta_l * (data[i][N - 3] - 4 * data[i][N - 2]) / 2. / h) / (
                    alpha_l + 3 * beta_l / 2. / h)


time_left = 0
time_right = 1
left = 0
right = 1
N = 21
T = 21
h = (right - left) / (N - 1)
tau = (time_right - time_left) / (T - 1)
function = [[0 for _ in range(N)] for _ in range(T)]
a = 0.5

phi = [0] * N
psi = [0] * N
for j in range(N):
    phi[j] = phi_func(left + h * j)
    psi[j] = psi_func(left + h * j)

f = [[0 for _ in range(N)] for _ in range(T)]
for i in range(T):
    for j in range(N):
        f[i][j] = f_func(left + h * j, time_left + tau * i)

gamma_0 = [0] * T
gamma_l = [0] * T
for i in range(T):
    gamma_0[i] = gamma_0_func(time_left + i * tau)
    gamma_l[i] = gamma_l_func(time_left + i * tau)

alpha_0 = 0
beta_0 = 1
alpha_l = 1
beta_l = -1

fig, ax = plt.subplots(subplot_kw=dict(projection='3d'))
x_grid = [left + i * h for i in range(N)]
y_grid = [time_left + i * tau for i in range(N)]
x_grid, y_grid = np.meshgrid(x_grid, y_grid)

wave_solving_2(time_left, time_right, left, right, function,
               a, N, T, phi, psi, f, gamma_0, gamma_l, alpha_0, beta_0, alpha_l, beta_l)
plt.figure(2)
plt.plot(y_grid, [function[i][15] for i in range(T)])

plt.figure(1)
function = np.array(function)
ax.plot_surface(x_grid, y_grid, function)

function = [[0 for _ in range(N)] for _ in range(T)]
wave_solving_1(time_left, time_right, left, right, function,
               a, N, T, phi, psi, f, gamma_0, gamma_l, alpha_0, beta_0, alpha_l, beta_l)
plt.figure(2)
plt.plot(y_grid, [function[i][15] for i in range(T)])

plt.figure(1)
function = np.array(function)
ax.plot_surface(x_grid, y_grid, function)

plt.figure(2)
plt.plot(y_grid, [real_func(left + 15 * h, time_left + i * tau) for i in range(T)])
plt.legend(["Первый порядок точности", "Второй порядок точности", "Реальная функция"])

plt.figure(1)
ax.plot_surface(x_grid, y_grid,
                np.array([[real_func(left + j * h, time_left + i * tau) for j in range(N)] for i in range(T)]))
plt.show()
