from numpy import cos, exp, pi
import matplotlib.pyplot as plt


def f(x):
    return cos(exp(x / 3.) / 10.)


def Lagrange(x, train_grid_size, train_grid):
    ans = 0
    for i in range(train_grid_size):
        tmp = 1
        for j in range(train_grid_size):
            if (i != j):
                tmp *= (x - train_grid[j]) / (train_grid[i] - train_grid[j])
        ans += tmp * f(train_grid[i])
    return ans


def Lagrange_equable(left, right, n, test_values):
    h = (right - left) / (n - 1)
    train_grid = [0] * n
    for i in range(n):
        train_grid[i] = left + i * h
    for i in range(n - 1):
        test_values[i] = Lagrange(left + (2 * i + 1) * h / 2, n, train_grid)
    ax1.scatter([left + (2 * i + 1) * h / 2 for i in range(n - 1)], test_values)


def Lagrange_Chebishev(left, right, n, test_values):
    train_grid = [0] * n
    for i in range(n):
        train_grid[n - i - 1] = (left + right) / 2. + (right - left) / 2. * cos((2 * i + 1) / 2. / n * pi)
    test_grid = [0] * (n - 1)
    for i in range(n - 1):
        test_grid[i] = (train_grid[i] + train_grid[i + 1]) / 2.
    for i in range(n - 1):
        test_values[i] = Lagrange(test_grid[i], n, train_grid)
    ax1.scatter(test_grid, test_values)


def Newton_train(train_grid_size, train_grid, koeffs):
    for i in range(train_grid_size):
        sep_sub = 0
        for j in range(i + 1):
            tmp = 1
            for k in range(i + 1):
                if (k != j):
                    tmp *= (train_grid[j] - train_grid[k])
            sep_sub += f(train_grid[j]) / tmp
        koeffs[i] = sep_sub


def Newton(x, train_grid_size, train_grid, koeffs):
    ans = 0
    for i in range(train_grid_size - 1, -1, -1):
        ans += koeffs[i]
        if (i != 0):
            ans *= (x - train_grid[i - 1])
    return ans


def Newton_equable(left, right, n, test_values):
    h = (right - left) / (n - 1)
    train_grid = [0] * n
    for i in range(n):
        train_grid[i] = left + i * h
    koeffs = [0] * n
    Newton_train(n, train_grid, koeffs)
    for i in range(n - 1):
        test_values[i] = Newton(left + (2 * i + 1) * h / 2., n, train_grid, koeffs)
    ax2.scatter([left + (2 * i + 1) * h / 2 for i in range(n - 1)], test_values)


def Newton_Chebishev(left, right, n, test_values):
    train_grid = [0] * n
    for i in range(n):
        train_grid[n - i - 1] = (left + right) / 2. + (right - left) / 2. * cos((2 * i + 1) / 2. / n * pi)
    test_grid = [0] * (n - 1)
    for i in range(n - 1):
        test_grid[i] = (train_grid[i] + train_grid[i + 1]) / 2.
    koeffs = [0] * n
    Newton_train(n, train_grid, koeffs)
    for i in range(n - 1):
        test_values[i] = Newton(test_grid[i], n, train_grid, koeffs)
    ax2.scatter(test_grid, test_values)


left = 0
right = 10
n = 31
grid = [left + i * (right - left) / (n - 1) for i in range(n)]
func = [f(grid[i]) for i in range(n)]

fig1, ax1 = plt.subplots()
ax1.plot(grid, func)

fig2, ax2 = plt.subplots()
ax2.plot(grid, func)

# Лагранж с равномерными узлами
test_values1 = [0] * (n - 1)
Lagrange_equable(left, right, n, test_values1)

# Лагранж с узлами Чебышева
test_values2 = [0] * (n - 1)
Lagrange_Chebishev(left, right, n, test_values2)

# Ньютон с равномерными узлами
test_values3 = [0] * (n - 1)
Newton_equable(left, right, n, test_values3)

# Ньютон с узлами Чебышева
test_values4 = [0] * (n - 1)
Newton_Chebishev(left, right, n, test_values4)

ax1.legend(["True function", "Lagrange_equable", "Lagrange_Chebishev"])
ax2.legend(["True function", "Newton_equable", "Newton_Chebishev"])
plt.show()
