from math import sin, log2, cos, exp, atan


# функция, нули которой необходимо найти
def f(x):
    return exp(sin(x / 2)) - atan(x) + 1

# производная
def f_deriv(x):
    return exp(sin(x / 2)) * cos(x / 2) / 2 - 1 / (x ** 2 + 1)


# численное нахождение производной функции
def centr_der(i, data, data_size, h):
    if i == 0:
        return (4 * data[1] - data[2] - 3 * data[0]) / 2 / h
    elif i == data_size - 1:
        return -(4 * data[data_size - 2] - data[data_size - 3] - 3 * data[data_size - 1]) / 2 / h
    return (data[i + 1] - data[i - 1]) / 2 / h


# метод бисекции
def bisection(left, right, eps):
    n = int(log2((right - left) / eps)) + 1
    # у данного метода можно рассчитать количество итераций и реализовать его циклом for, без while и prev_x
    for i in range(n):
        x = (right + left) / 2.
        if f(x) * f(right) < 0:
            left = x
        else:
            right = x
    return x


# метод Ньютона
def newton_modif(f, f0, x0, eps):
    '''
        f - левая часть решаемого уравнения
        f0 - производная левой части в точке x0
        x0 - начальное приближение
        eps - точность решения
    '''
    x1 = x0 - (f(x0)) / f0  # считаем первое приближение

    crit1 = abs(x1 - x0) > eps  # флаг для отслеживания первого критерия останова
    crit2 = abs(f(x1) - f(x0)) > eps  # флаг для отслеживания второго критерия останова

    while (crit1):
        x0, x1 = x1, x1 - (f(x1)) / f0  # считаем новые приближенные значения корня
        crit1 = abs(x1 - x0) > eps  # проверяем выполнение критериев останова
        crit2 = abs(f(x1) - f(x0)) > eps

    return x1


# метод хорд
def chords(left, right, eps):
    prev_x = left
    x = -f(left) / (f(right) - f(left)) * (right - left)
    while abs(x - prev_x) > eps:
        prev_x = x
        if f(x) * f(left) < 0:
            right = x
        else:
            left = x
        x = -f(left) / (f(right) - f(left)) * (right - left)
    return x


left0 = 0
right0 = 10
n0 = 10000
h = (right0 - left0) / (n0 - 1)

data = [f(left0 + i * h) for i in range(n0)]

der_data = [centr_der(i, data, n0, h) for i in range(n0)]

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


eps_s = [1.0E-3, 1.0E-6, 1.0E-9]
for j in range(3):
    eps = eps_s[j]
    print(f"Для точности {eps}:")

    print("Метод бисекции:")
    for i in range(len(clear_starts)):
        left = clear_starts[i]
        right = clear_ends[i]

        x = bisection(left, right, eps)
        print(f"Корень номер {i + 1}: {x:.{(3 * j + 4)}}")
        # |x - previous_x| < eps = True, так как является критерием останова для метода
        print("|x - previous_x| < eps = True, |f(x)| < eps = ", end='')
        if (abs(f(x)) < eps):
            print("True")
        else:
            print("False")
        print()

    print("Модифицированный метод Ньютона:")
    for i in range(len(clear_starts)):
        left = clear_starts[i]
        right = clear_ends[i]
        x0 = (right + left) / 2.

        x = newton_modif(f, f_deriv(x0), x0, eps)
        print(f"Корень номер {i + 1}: {x:.{(3 * j + 4)}}")
        # |x - previous_x| < eps = True, так как является критерием останова для метода
        print("|x - previous_x| < eps = True, |f(x)| < eps = ", end='')
        if (abs(f(x)) < eps):
            print("True")
        else:
            print("False")
        print()

    print("Метод хорд:")
    for i in range(len(clear_starts)):
        left = clear_starts[i]
        right = clear_ends[i]

        x = chords(left, right, eps)
        print(f"Корень номер {i + 1}: {x:.{(3 * j + 4)}}")
        # |x - previous_x| < eps = True, так как является критерием останова для метода
        print("|x - previous_x| < eps = True, |f(x)| < eps = ", end='')
        if (abs(f(x)) < eps):
            print("True")
        else:
            print("False")
        print("\n\n")
