from polynomials import *
from differentiation import *
from integration import *
from functions import *
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


# def taylor(f, a, b, alpha, h, n):
#     """returns a list of points (x, y) which approximate the solution to the IVP y' = f(x, y) with y(a) = alpha on the interval [a, b] with step size h using n'th order taylor method"""

#     x = a
#     y = alpha

#     result = [(x, y)]

#     while x < b:
#         y += h*sum([Func(f).mult_Der(i)(x, y)/(factorial(i)) for i in range(1, n+1)])
#         x += h
#         result.append((x, y))

#     return result

# def runge_kutta(f, a, b, alpha, h, n):
