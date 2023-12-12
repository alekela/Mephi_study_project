from math import *
import matplotlib.pyplot as plt


def func(x):
    return cosh(exp(-pow(x, 2)))


def right_der(i, data, data_size, h):
    if (i != data_size - 1):
        return (data[i + 1] - data[i]) / h
    return (data[i] - data[i - 1]) / h


def left_der(i, data, data_size, h):
    if (i != 0):
        return (data[i] - data[i - 1]) / h
    return (data[i + 1] - data[i]) / h


def centr_der(i, data, data_size, h):
    if (i == 0):
        return (4 * data[1] - data[2] - 3 * data[0]) / 2 / h
    elif (i == data_size - 1):
        return -(4 * data[data_size - 2] - data[data_size - 3] - 3 * data[data_size - 1]) / 2 / h
    return (data[i + 1] - data[i - 1]) / 2 / h


def second_der2(i, data, data_size, h):
    if (i == 0):
        return (data[0] - 2 * data[1] + data[2]) / h / h
    elif (i == data_size - 1):
        return (data[data_size - 1] - 2 * data[data_size - 2] + data[data_size - 3]) / h / h
    return (data[i + 1] - 2 * data[i] + data[i - 1]) / h / h


def second_der4(i, data, data_size, h):
    if (i == 0):
        return (data[1] * (-77. / 12) + data[2] * (107. / 12) + data[3] * (-13. / 2) +
                data[4] * (61. / 24) + data[5] * (-5. / 12) + data[0] * (15. / 8)) / h / h * 2.
    elif (i == 1):
        return (data[0] * 5. / 12 + data[2] * (-1. / 6) + data[3] * (7. / 12) +
                data[4] * (-0.25) + data[5] * (1. / 24) - data[1] * (5. / 8)) * 2. / h / h
    elif (i == data_size - 1):
        return (data[data_size - 2] * (-77. / 12) + data[data_size - 3] * (107. / 12) + data[data_size - 4] * (
                -13. / 2) +
                data[data_size - 5] * (61. / 24) + data[data_size - 6] * (-5. / 12) + data[data_size - 1] * (
                        15. / 8)) / h / h * 2.
    elif (i == data_size - 2):
        return (data[data_size - 1] * 5. / 12 + data[data_size - 3] * (-1. / 6) + data[data_size - 4] * (7. / 12) +
                data[data_size - 5] * (-0.25) + data[data_size - 6] * (1. / 24) - data[data_size - 2] * (
                        5. / 8)) * 2. / h / h
    return (16 * (data[i + 1] + data[i - 1]) - (data[i + 2] + data[i - 2]) - 30 * data[i]) / 12. / h / h


n = 1001
left = -5
right = 5

h = (right - left) / (n - 1)
data = [func(left + i * h) for i in range(n)]
second_der_data = []
for i in range(n):
    second_der_data.append(second_der4(i, data, n, h))

plt.plot([left + i * h for i in range(n)], second_der_data)
plt.show()
