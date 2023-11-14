from linear_algebra import *
from diff_eq_IVP import *
from functions import *
from math import *
import numpy


def lin_fin_diff(p, q, r, a, b, alpha, beta, N):
    """returns a list of N+2 points (including (a, alpha), (b, beta)) (x, y) which approximate the solution to the BVP y'' + p(x)y' + q(x)y = r(x) with y(a) = alpha and y(b) = beta using method of finite differences for linear BVPs"""

    h = (b - a)/(N+1)

    B = vector([-(h**2)*r(a+h) + alpha*(1 + (h)*p(a+h)/2)]+[-(h**2)*r(a + (i+2)*h) for i in range(N-2)] + [-(h**2)*r(b-h) + beta*(1 - (h)*p(b-h)/2)])

    A = [[2 + (h**2)*q(a+h), -1 + (h)*p(a+h)/2] + [0 for _ in range(N-2)]]
    
    for i in range(2, N):
        A.append([0 for _ in range(i-2)] + [-1 - (h)*p(a + i*h)/2, 2 + (h**2)*q(a + i*h), -1 + (h)*p(a + i*h)/2] + [0 for _ in range(N-i-1)])

    A.append([0 for _ in range(N-2)] + [-1 - (h)*p(b-h)/2, 2 + (h**2)*q(b-h)])

    A = matrix(A)
    W = gauss_elim(A, B)

    return [(a, alpha)] + [(a + (i+1)*h, w) for i, w in enumerate(W.coordinates)] + [(b, beta)]
 
def nonlin_fin_diff(f, a: float, b: float, alpha: float, beta: float, W_0: vector, N: int, M: int)->list[tuple]:
    """returns a list of N+2 points (including (a, alpha), (b, beta)) (x, y) which approximate the solution to the BVP y'' = f(x, y, y') with y(a) = alpha and y(b) = beta using method of finite differences for nonlinear BVPs"""

    h = (b - a)/(N+1)

    f_y = lambda x, y, y1: Func(lambda z: f(x, z, y1)).Der()(y)
    f_y1 = lambda x, y, y1: Func(lambda z: f(x, y, z)).Der()(y1)

    for _ in range(M):
        J = [[2 + (h**2)*f_y(a+h, W_0[0], (W_0[1] - alpha)/h), -1 + (h)*f_y1(a+h, W_0[0], (W_0[1] - alpha)/h)/2] + [0 for _ in range(N-2)]]
        
        for i in range(1, N-1):
            J.append([0 for _ in range(i-1)] + [-1 - (h)*f_y1(a + (i+1)*h, W_0[i], (W_0[i+1] - W_0[i-1])/(2*h))/2, 2 + (h**2)*f_y(a + (i+1)*h, W_0[i], (W_0[i+1] - W_0[i-1])/(2*h)), -1 + (h)*f_y1(a + (i+1)*h, W_0[i], (W_0[i+1] - W_0[i-1])/(2*h))/2] + [0 for _ in range(N-i-2)])

        J.append([0 for _ in range(N-2)] + [-1 - (h)*f_y1(b-h, W_0[N-1], (beta - W_0[N-2])/h)/2, 2 + (h**2)*f_y(b-h, W_0[N-1], (beta - W_0[N-2])/h)])

        W_0 -= matrix(J).inverse()*(vector([f(a+h, W_0[0], (W_0[1] - alpha)/(2*h))] + [f(a + (i+1)*h, W_0[i], (W_0[i+1] - W_0[i-1])/(2*h)) for i in range(1, N-1)] + [f(b-h, W_0[N-1], (beta - W_0[N-2])/(2*h))]))

    return [(a, alpha)] + [(a + (i+1)*h, w) for i, w in enumerate(W_0.coordinates)] + [(b, beta)]



# print(lin_fin_diff(lambda x: -2/x, lambda y: 2/y**2, lambda z: sin(log(z))/z**2, 1, 2, 1, 2, 9))
print(nonlin_fin_diff(lambda x, y, y1: 4 + 0.25*x**3 - 0.125*y*y1, 1, 3, 17, 43/3, vector([1 + i/21 for i in range(1, 21)]), 20, 10))