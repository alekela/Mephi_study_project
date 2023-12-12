from math import *
import matplotlib.pyplot as plt
import numpy as np


def func(x):
    return sqrt(2 - x * x)


def rect(n, left, right):
    h = (right - left) / n
    summ = 0
    for i in range(n):
        summ += func(left + i * h + h / 2.)
    return summ * h


def trapezoid(n, left, right):
    h = (right - left) / n
    summ = 0
    a = left
    for i in range(n):
        summ += (func(left + h * (i + 1)) + func(left + i * h))
    return summ * h / 2


def simpson(n, left, right):
    h = (right - left) / n
    summ = 0
    for i in range(n - 1):
        summ += (func(left + i * h) + 4.0 * func(left + (i + 1) * h) + func(left + (i + 2) * h))
    return summ * h / 6.


left = 0
right = 1
true_value = 1.285398163397448
ns = [10, 100, 1000, 10000, 100000, 1000000]
logns = [i for i in range(1, 7)]
rect_errors = []
trapezoid_errors = []
simpson_errors = []
for i in range(6):
    n = ns[i]
    rect_errors.append(abs(rect(n, left, right) - true_value))
    trapezoid_errors.append(abs(trapezoid(n, left, right) - true_value))
    simpson_errors.append(abs(simpson(n, left, right) - true_value))

rect_errors = np.array(rect_errors)
trapezoid_errors = np.array(trapezoid_errors)
simpson_errors = np.array(simpson_errors)
ns = np.array(ns)

plt.plot(np.log(ns), np.log(rect_errors))
plt.plot(np.log(ns), np.log(trapezoid_errors))
plt.plot(np.log(ns), np.log(simpson_errors))
plt.show()
