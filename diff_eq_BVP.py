from linear_algebra import *
from diff_eq_IVP import *
from functions import *
from math import *
from numpy.linalg import *


def lin_fin_diff(p, q, r, a, b, alpha, beta, N):
    """returns a list of N+2 points (including (a, alpha), (b, beta)) (x, y) which approximate the solution to the BVP y'' + p(x)y' + q(x)y = r(x) with y(a) = alpha and y(b) = beta using method of finite differences for linear BVPs"""

    h = (b - a)/(N+1)

    B = vector([-(h**2)*r(a+h) + (1 + (h)*p(a+h)/2)*alpha]+[-(h**2)*r(a + (i+2)*h) for i in range(N-2)] + [-(h**2)*r(b-h) + (1 - (h)*p(b-h)/2)*beta])

    A = [[2 + (h**2)*q(a+h), -1 + (h)*p(a+h)/2] + [0 for _ in range(N-2)]]

    for i in range(2, N):
        A.append([0 for _ in range(i-2)] + [-1 + (h)*p(a + i*h)/2, 2 + (h**2)*q(a + i*h), -1 - (h)*p(a + i*h)/2] + [0 for _ in range(N-i-1)])

    A.append([0 for _ in range(N-2)] + [-1 - (h)*p(b-h)/2, 2 + (h**2)*q(b-h)])

    A = matrix(A)
    W = gauss_elim(A, B)

    print(solve(A.rows, B.coordinates))

    return [(a, alpha)] + [(a + (i+1)*h, w) for i, w in enumerate(B.coordinates)] + [(b, beta)]
 
def nonlin_fin_diff(f, a: float, b: float, alpha: float, beta: float, W_0: vector, N: int)->list[tuple]:
    """returns a list of N+2 points (including (a, alpha), (b, beta)) (x, y) which approximate the solution to the BVP y'' = f(x, y, y') with y(a) = alpha and y(b) = beta using method of finite differences for nonlinear BVPs"""

    h = (b - a)/(N+1)

    f_y = lambda x, y, y1: Func(lambda y: f(x, y, y1)).Der().f(y)
    f_y1 = lambda x, y, y1: Func(lambda y1: f(x, y, y1)).Der().f(y1)

    f_yi = lambda i: f_y(a + i*h, W_0.coordinates[i], (W_0.coordinates[i+1] - W_0.coordinates[i-1])/(2*h))
    f_y1i = lambda i: f_y1(a + i*h, W_0.coordinates[i], (W_0.coordinates[i+1] - W_0.coordinates[i-1])/(2*h))

    J = [[0 for _ in range(N)] for _ in range(N)]

    for k in range(100):
        for i in range(N):
            for j in range(i+1):
                J[i][j] = 0
            
            J[i][i+1] = -1 + (h/2)*f_y1i(i)
            J[i][i] = 2 + (h**2)*f_yi(i)
            J[i][i-1] = -1 - (h/2)*f_y1i(i)

            for j in range(i+2, N):
                J[i][j] = 0

        W_0 -= matrix(J).inverse()(vector([f(w) for w in W_0.coordinates]))

    return [(a, alpha)] + [(a + (i+1)*h, w) for i, w in enumerate(W_0.coordinates)] + [(b, beta)]

print(lin_fin_diff(lambda x: -2/x, lambda x: 2/x**2, lambda x: sin(log(x))/x**2, 1, 2, 1, 2, 9))