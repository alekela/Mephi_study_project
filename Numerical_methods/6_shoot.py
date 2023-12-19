from math import sin, log2, cos, exp, atan, tanh, cosh
import matplotlib.pyplot as plt


# Ответ
def answer(x):
    return sin(x) + tanh(x / 2)


# Производная ответа
def answer_deriv(x):
    return cos(x) + 1 / (2 * cosh(0.5 * x) * cosh(0.5 * x))


def W(x, u, v):
    return -sin(x) - 0.5 * sin(2 * x) - v * tanh(x / 2.) + u * cos(x)


def V(v):
    return v


# метод хорд
def chords(A, B, realB):
    return A[1] - ((B[1] - realB) * (A[1] - A[0])) / ((B[1] - realB) - (B[0] - realB))


def Runge(x, u, v, u0, v0, h):
    u.append(u0)
    v.append(v0)
    for i in range(len(x) - 1):
        q0 = W(x[i], u[i], v[i])
        k0 = V(v[i])
        q1 = W(x[i] + h / 2, u[i] + k0 * h / 2, v[i] + q0 * h / 2)
        k1 = V(v[i] + q0 * h / 2)
        q2 = W(x[i] + h / 2, u[i] + k1 * h / 2, v[i] + q1 * h / 2)
        k2 = V(v[i] + q1 * h / 2)
        q3 = W(x[i] + h, u[i] + k2 * h, v[i] + q2 * h)
        k3 = V(v[i] + q2 * h)
        u.append(u[i] + h / 6 * (k0 + 2 * k1 + 2 * k2 + k3))
        v.append(v[i] + h / 6 * (q0 + 2 * q1 + 2 * q2 + q3))


def main():
    a = 0.
    b = 1.
    h = 0.05
    right = 0.9335
    n = int((b - a) / h) + 1

    A = [9876543456789, 23456789876543]
    B = [0, 0]

    u0 = -1.5 + A[0]
    du0 = A[0]

    x_grid = [(a + i * h) for i in range(n)]
    u_runge = []
    v_runge = []

    Runge(x_grid, u_runge, v_runge, u0, du0, h)

    B[0] = v_runge[int((b - a) / h)]

    u_runge = []
    v_runge = []

    eps = 1.0E-9

    while True:

        u0 = -1.5 + A[1]
        du0 = A[1]

        u_runge = []
        v_runge = []

        Runge(x_grid, u_runge, v_runge, u0, du0, h)

        B[1] = v_runge[int((b - a) / h)]

        if (abs(B[1] - right) < eps): break

        tmp = chords(A, B, right)  # Ak+1

        A[0] = A[1]
        A[1] = tmp
        B[0] = B[1]

        u_runge = []
        v_runge = []

    plt.figure(1)
    plt.plot(x_grid, u_runge, label='Numerical')
    plt.plot(x_grid, list(map(answer, x_grid)), label='Theoretical')
    plt.legend()
    plt.figure(2)
    plt.plot(x_grid, v_runge, label='Numerical')
    plt.plot(x_grid, list(map(answer_deriv, x_grid)), label='Theoretical')
    plt.legend()
    plt.show()

    print("Max function error:", round(max(abs(u_runge[i] - answer(x_grid[i])) for i in range(n)), 8))
    print("Max function derivative error:", round(max(abs(v_runge[i] - answer_deriv(x_grid[i])) for i in range(n)), 8))


main()
