from math import *
from functions import *
from polynomials import *

def forward_der(f, x, h, n):
    """returns the derivative of f at x using forward difference metthod of order n with step size h"""

    p = newt_interpol(f, x, x + n*h, n)

    return p.Der()(x)

def backward_der(f, x, h, n):
    """returns the derivative of f at x using backward difference metthod of order n with step size h"""

    p = newt_interpol(f, x - n*h, x, n)

    return p.Der()(x)

def central_der(f, x, h, n):
    """returns the derivative of f at x using central difference method of order n with step size h"""

    p = newt_interpol(f, x - n*floor(h/2), x + n*ceil(h/2), n)

    return p.Der()(x)

print(forward_der(lambda x : x**2, 1, 1, 10))
