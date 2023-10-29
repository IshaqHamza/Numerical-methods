from polynomials import *

def closed_cotes(f, a, b, n):
    """returns the integral of f from a to b using closed cotes method of order n"""

    p = newt_interpol(f, a, b, n)

    return p.Int()(b) - p.Int()(a)

def rect(f, a, b):
    """uses rectangle method to return the integral of f from a to b"""
    return closed_cotes(f, a, b, 0)

def trap(f, a, b):
    """uses trapezoidal method to return the integral of f from a to b"""
    return closed_cotes(f, a, b, 1)

def simp(f, a, b):
    """uses simpson's method to return the integral of f from a to b"""
    return closed_cotes(f, a, b, 2)

def simp_3_8(f, a, b):
    """uses simpson's 3/8 method to return the integral of f from a to b"""
    return closed_cotes(f, a, b, 3)

def boole(f, a, b):
    """uses boole's method to return the integral of f from a to b"""
    return closed_cotes(f, a, b, 4)

def open_cotes(f, a, b, n):
    """returns the integral of f from a to b using open cotes method of order n"""

    h = (b - a)/(n + 2)

    p = newt_interpol(f, a+h, b-h, n)

    return p.Int()(b) - p.Int()(a)

def mid_point(f, a, b):
    """uses mid point method to return the integral of f from a to b"""
    return open_cotes(f, a, b, 0)

def int_comp(f, a, b, n, rule):
    """returns the integral of f from a to b using composite rule with n intervals"""
    h = (b - a)/n
    return sum([rule(f, a + i*h, a + (i+1)*h) for i in range(n)])



def gauss_quadrature(f, a, b, n):
    """uses n'th order gaussian quadrature to evaluate integral of f from a to b"""

    #first we convert the integral from (a, b) to (-1, 1)   

    def g(x):
        return f((b - a)*x/2 + (b + a)/2)
    
    #then we use roots of n'th legendre polynomial to find the nodes and weights

    nodes = [legendre_root(n, k) for k in range(1, n+1)]
    weights = [2/((1 - nodes[k-1]**2)*(legendre_monic(n).Der()(nodes[k-1]))**2) for k in range(1, n+1)]

    #finally we use the nodes and weights to evaluate the integral

    return sum([weights[k-1]*g(nodes[k-1]) for k in range(1, n+1)])*(b - a)/2


def double_int(f, a, b, c, d, n, m, x_rule, y_rule):
    """returns the integral of f from a to b and c to d using n intervals in x and m intervals in y"""
    return int_comp(lambda y : int_comp(lambda x : f(x, y), a, b, n, x_rule), c, d, m, y_rule)


