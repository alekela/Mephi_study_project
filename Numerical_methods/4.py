from math import sin, log2


def func(x):
    return sin(x)


def centr_der(i, data, data_size, h):
    if i == 0:
        return (4 * data[1] - data[2] - 3 * data[0]) / 2 / h
    elif i == data_size - 1:
        return -(4 * data[data_size - 2] - data[data_size - 3] - 3 * data[data_size - 1]) / 2 / h
    return (data[i + 1] - data[i - 1]) / 2 / h


def bisection(left, right, eps):
    n = int(log2((right - left) / eps)) + 1
    for i in range(n):
        x = (right + left) / 2.
        if func(x) * func(right) < 0:
            left = x
        else:
            right = x
    return x


left0 = 1
right0 = 8
n0 = 8
h = (right0 - left0) / (n0 - 1)

data = [0] * n0
for i in range(n0):
    data[i] = func(left0 + i * h)

der_data = [0] * n0
for i in range(n0):
    der_data[i] = centr_der(i, data, n0, h)

starts = []
ends = []
starts.append(0)
for i in range(1, n0):
    if der_data[i - 1] * der_data[i] < 0:
        starts.append(i)
        ends.append(i - 1)

ends.append(n0 - 1)
clear_starts = []
clear_ends = []
for i in range(len(starts)):
    if data[starts[i]] * data[ends[i]] < 0:
        clear_starts.append(left0 + starts[i] * h)
        clear_ends.append(left0 + ends[i] * h)

eps = float(input("Введите точность:\n"))

for i in range(len(clear_starts)):
    left = clear_starts[i]
    right = clear_ends[i]

    x = bisection(left, right, eps)
    print(f"Найденный корень из промежутка от {clear_starts[i]} до {clear_ends[i]} равен {x}")
