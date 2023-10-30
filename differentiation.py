from math import *
from functions import *
from polynomials import *

def forward_der(f, x, h, n):
    """returns the derivative of f at x using forward difference method of order n with step size h"""

    p = newt_interpol(f, x, x + n*h, n)

    return p.Der()(x)

def backward_der(f, x, h, n):
    """returns the derivative of f at x using backward difference method of order n with step size h"""

    p = newt_interpol(f, x - n*h, x, n)

    return p.Der()(x)

def central_der(f, x, h, n):
    """returns the derivative of f at x using central difference method of order n with step size h"""

    p = newt_interpol(f, x - (h/2)*floor(n/2), x + (h/2)*ceil(n/2), n)

    return p.Der()(x)

def mult_der(f, x, n, h = 0.00001, N = 100):
    """returns the n'th derivative of f at x using central differences with step size h"""

    p = newt_interpol(f, x - (h/2)*floor(N/2), x + (h/2)*ceil(N/2), n)

    return p.mult_Der(n)(x)
