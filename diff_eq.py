from polynomials import *
from differentiation import *
from integration import *
from functions import *
from scipy.misc import derivative
from roots import *


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

def RK4(f, a, b, alpha, h):     # hoping its correct (mostly not ig)
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



