from polynomials import *
from differentiation import *
from integration import *
from functions import *
from scipy.misc import derivative
from roots import *
from math import *
from linear_algebra import *

def euler(f, a, b, alpha, h):
    """returns a list of points (x, y) which approximate the solution to the IVP y' = f(x, y) with y(a) = alpha on the interval [a, b] with step size h"""

    x = a
    y = alpha

    result = [(x, y)]

    while x < b:
        y += h*f(x, y)
        x += h
        result.append((x, y))

    return result

def T(f, n, t, w, h):
    return sum([(h**i/factorial(i))*derivative(lambda t: f(t, w), t, n = i, dx=0.00001, order=n+1) for i in range(n)])

def taylor(n, f, a, b, alpha, h):
    """returns a list of points (x, y) which approximate the solution to the IVP y' = f(x, y) with y(a) = alpha on the interval [a, b] with step size h using n'th order taylor method"""

    x = a
    y = alpha

    result = [(x, y)]

    while x < b:
        y += T(f, n, x, y, h)
        x += h
        result.append((x, y))

    return result

def RK2(f, a, b, alpha, h):     #seems correct
    """returns a list of points (x, y) which approximate the solution to the IVP y' = f(x, y) with y(a) = alpha on the interval [a, b] with step size h using RK-2 method"""

    x = a
    y = alpha

    result = [(x, y)]

    while x < b:
        k1 = f(x, y)
        k2 = f(x + h/2, y + h*k1/2)
        y += h*k2
        x += h
        result.append((x, y))

    return result

def RK4(f, a, b, alpha, h):     # hoping its correct 
    """returns a list of points (x, y) which approximate the solution to the IVP y' = f(x, y) with y(a) = alpha on the interval [a, b] with step size h using RK-4 method"""

    x = a
    y = alpha

    result = [(x, y)]

    while x < b:
        k1 = f(x, y)
        k2 = f(x + h/2, y + h*k1/2)
        k3 = f(x + h/2, y + h*k2/2)
        k4 = f(x + h, y + h*k3)
        y += h*(k1 + 2*k2 + 2*k3 + k4)/6
        x += h
        result.append((x, y))

    return result

def trap_IVP(f, a, b, alpha, h):
    """returns a list of points (x, y) which approximate the solution to the IVP y' = f(x, y) with y(a) = alpha on the interval [a, b] with step size h using trapezoidal method"""

    x = a
    y = alpha

    result = [(x, y)]

    while x < b:
        y = fixed_point(lambda t: y + h*(f(x, y) + f(x + h, t))/2, y, 0.00001, 100)
        x += h
        result.append((x, y))

    return result

def neg_bin_int(k):
    """returns the value (-1)^k times the integral of comb(-s, k) from 0 to 1)"""
    P = Poly([1])

    for i in range(k):
        P = Poly_mul(P, Poly([i, 1]))
    P = P.Int()

    result = P(1)

    for i in range(1, k+1):
        result /= i

def adam_bash(f, a, b, alpha, h, m, W):
    """uses adam-bashforth m-step method to solve the IVP y' = f(x, y) with y(a) = alpha on the interval [a, b] with step size h, given initial mesh points W = [(a, alpha), (a + h, y1), ..., (a + (m-1)*h, y(m-1))]"""

    F = [f(W[i][0], W[i][1]) for i in range(m)]

    def f_lol(i):
        return F[i]
    
    f_lel = Func(f_lol)
    
    x = a + (m-1)*h
    y = W[-1]

    result = W[:]

    i = m
    while x < b:
        y = y + h*sum([neg_bin_int(k)*f_lel.backward_diff_undiv(i, 1, k) for k in range(m)])
        x += h
        result.append((x, y))
        F.append(f(x, y))
        i += 1
    
    return result

def adam_mutton(f, a, b, alpha, h, m, W):
    """uses adam-moulton m-step method to solve the IVP y' = f(x, y) with y(a) = alpha on the interval [a, b] with step size h, given initial mesh points W = [(a, alpha), (a + h, y1), ..., (a + (m-1)*h, y(m-1))]"""

    F = [f(W[i][0], W[i][1]) for i in range(m)]

    def f_lol(i):
        return F[i]
    
    f_lel = Func(f_lol)
    
    x = a + (m-1)*h
    y = W[-1]

    result = W[:]

    i = m
    while x < b:
        y = fixed_point(lambda t: y + h*(sum([neg_bin_int(k)*f_lel.backward_diff_undiv(i, 1, k) for k in range(m)])) + neg_bin_int(m)*f(x + h, t), y, 0.00001, 100)
        x += h
        result.append((x, y))
        F.append(f(x, y))
        i += 1
    
    return result

def pred_corr(f, a, b, alpha, h, m, W):
    """uses m step predictor corrector method for solving the IVP y' = f(x, y) with y(a) = alpha on the interval [a, b] with step size h, given initial mesh points W = [(a, alpha), (a + h, y1), ..., (a + (m-1)*h, y(m-1))]"""
    
    F = [f(W[i][0], W[i][1]) for i in range(m)]

    def f_lol(i):
        return F[i]
    
    f_lel = Func(f_lol)
    
    x = a + (m-1)*h
    y = W[-1]

    result = W[:]

    i = m
    while x < b:
        y = y + h*sum([neg_bin_int(k)*f_lel.backward_diff_undiv(i, 1, k) for k in range(m)])
        y = y + h*(sum([neg_bin_int(k)*f_lel.backward_diff_undiv(i, 1, k) for k in range(m)])) + neg_bin_int(m)*f(x + h, y)
        x += h
        result.append((x, y))
        F.append(f(x, y))
        i += 1
    
    return result


