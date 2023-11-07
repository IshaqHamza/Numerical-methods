from math import *
import typing
import matplotlib.pyplot as plt
from roots import *
import numpy.polynomial as P
class Poly():
    def __init__(self, coefficients: list[float]): # coefficients are in the form [a_0, a_1, ..., a_n] for a_n*x^n + ... + a_1*x + a_0
        self.coefficients = coefficients
        self.degree = len(coefficients) - 1

    def __call__(self, x : [float]) -> float:  
        """returns the value of a polynomial at a point x, using Horner's method"""
        value = 0
        for i in range(len(self.coefficients)-1, -1, -1):
            value = self.coefficients[i] + x*value
        return value
    
    def copy(self):
        return Poly(self.coefficients.copy())
    
    def print(self):
        """prints the polynomial"""
        print(self.coefficients)
    
    def plot(self, a, b):
        x = [a + (b - a)*i/1000 for i in range(1001)]
        y = [self(point) for point in x]
        plt.plot(x, y)
        plt.show()

    def __add__(self, other):
        """returns the sum of two polynomials"""
        if self.degree > other.degree:
            return Poly([self.coefficients[i] + other.coefficients[i] if i <= other.degree else self.coefficients[i] for i in range(self.degree + 1)])
        return Poly([self.coefficients[i] + other.coefficients[i] if i <= self.degree else other.coefficients[i] for i in range(other.degree + 1)])

    def __sub__(self, other):
        """returns the difference of two polynomials"""
        if self.degree > other.degree:
            return Poly([self.coefficients[i] - other.coefficients[i] if i <= other.degree else self.coefficients[i] for i in range(self.degree + 1)])
        return Poly([self.coefficients[i] - other.coefficients[i] if i <= self.degree else -other.coefficients[i] for i in range(other.degree + 1)])

    def __mul__(self, other):
        """returns the product of two polynomials"""
        coefficients = [0 for _ in range(self.degree + other.degree + 1)]
        for i in range(self.degree + 1):
            for j in range(other.degree + 1):
                coefficients[i + j] += self.coefficients[i]*other.coefficients[j]
        return Poly(coefficients)

    def __truediv__(self, other):
        """returns the quotient of when self is divided by other"""


    def div_diff(self, xs: list[float]) -> float:
        """returns the divided difference of a polynomial at the points in xs"""
        if len(xs) == 1:
            return self(xs[0])
        return (self.div_diff(xs[1:]) - self.div_diff(xs[:-1])) / (xs[-1] - xs[0])

    def Der(self):
        """returns the derivative of a polynomial"""
        return Poly([self.coefficients[i]*i for i in range(1, len(self.coefficients))]) if self.degree > 0 else Poly([0])
    
    def Int(self):
        """returns the integral of a polynomial"""
        return Poly([0] + [self.coefficients[i]/(i + 1) for i in range(len(self.coefficients))])
    
    def mult_Der(self, n):
        """returns the nth derivative of a polynomial"""
        if n == 0:
            return self
        return self.Der().mult_Der(n-1)
    
    def root(self, a : float, b : float, tol = 0.00001, N = 10000):
        """returns a root of a polynomial in the interval (a, b) if it exists"""
        return regula_falsi(self, a, b, tol, N)


def Poly_add(f : Poly, g : Poly) -> Poly:
    """returns the sum of two polynomials"""
    if f.degree > g.degree:
        return Poly([f.coefficients[i] + g.coefficients[i] if i <= g.degree else f.coefficients[i] for i in range(f.degree + 1)])
    return Poly([f.coefficients[i] + g.coefficients[i] if i <= f.degree else g.coefficients[i] for i in range(g.degree + 1)])

def Poly_sub(f : Poly, g : Poly) -> Poly:
    """returns the difference of two polynomials"""
    if f.degree > g.degree:
        return Poly([f.coefficients[i] - g.coefficients[i] if i <= g.degree else f.coefficients[i] for i in range(f.degree + 1)])
    return Poly([f.coefficients[i] - g.coefficients[i] if i <= f.degree else -g.coefficients[i] for i in range(g.degree + 1)])

def Poly_mul(f : Poly, g : Poly) -> Poly:
    """returns the product of two polynomials"""
    coefficients = [0 for _ in range(f.degree + g.degree + 1)]
    for i in range(f.degree + 1):
        for j in range(g.degree + 1):
            coefficients[i + j] += f.coefficients[i]*g.coefficients[j]
    return Poly(coefficients)

def Poly_div(f : Poly, g : Poly) -> Poly:
    """returns the quotient of when f is divided by g"""
    
    if g.degree == 0:
        return Poly([f.coefficients[i]/g.coefficients[0] for i in range(f.degree + 1)])
    
    if f.degree < g.degree:
        return Poly([0])
    
    coefficients = [0 for _ in range(f.degree - g.degree + 1)]
    h = f.copy()
    for i in range(f.degree - g.degree, -1, -1):
        coefficients[i] = h.coefficients[i + g.degree]/g.coefficients[g.degree]
        for j in range(g.degree + 1):
            h.coefficients[i + j] -= coefficients[i]*g.coefficients[j]
    return Poly(coefficients)

    

def Poly_rem(f : Poly, g : Poly) -> Poly:
    """returns the remainder of when f is divided by g using Euclidean division"""
    if f.degree < g.degree:
        return f
    
    h = f.copy()
    q = Poly([0])

    while h.degree >= g.degree:
        q = Poly_add(q, Poly([h.coefficients[0]/g.coefficients[0]] + [0 for _ in range(h.degree - g.degree)]))
        h = Poly_sub(h, Poly_mul(q, g))

    return h

def Poly_gcd(f : Poly, g : Poly) -> Poly:
    """returns the greatest common divisor of two polynomials using Euclidean division"""
    if f.degree < g.degree:
        return Poly_gcd(g, f)
    
    h = f.copy()

    while h.degree >= g.degree:
        h = Poly_rem(h, g)

    return h

def Poly_lcm(f : Poly, g : Poly) -> Poly:
    """returns the least common multiple of two polynomials"""
    return Poly_mul(f, g)/Poly_gcd(f, g)

def Poly_eq(f : Poly, g : Poly) -> bool:
    """returns whether two polynomials are equal"""
    return f.coefficients == g.coefficients


def legendre_poly(n):
    """returns the legendre polynomial of degree n"""
    if n == 0:
        return Poly([1])
    if n == 1:
        return Poly([0, 1])

    P_n_1 = legendre_poly(n-1)
    P_n_2 = legendre_poly(n-2)

    A = Poly_mul(P_n_1, Poly([2*n - 1]))
    A = Poly([c/n for c in A.coefficients])
    A = Poly_mul(A, Poly([0, 1]))

    B = Poly_mul(P_n_2, Poly([n - 1]))
    B = Poly([c/n for c in B.coefficients])

    return Poly_sub(A, B)



def legendre_monic(n):
    P = legendre_poly(n)
    return Poly_div(P, Poly([P.coefficients[-1]]))

def legendre_root(n, k):
    """returns the kth root of the legendre polynomial of degree n"""
    if not n:
        if k != 1:
            return None
        return 0
    if k > n:
        return None
    e_nk = (1 - 1/(8*n**2) + 1/(8*n**3))*cos(pi*(4*k - 1)/(4*n + 2))

    i = 0
    p = e_nk

    f = legendre_poly(n)
    df = f.Der()

    while i < 6:
        p1 = p - f(p)/df(p)
        if abs(p1 - p) < 0.00001:
            return p1
        p = p1
        i += 1


def legendre_roots(n):
    """returns the roots of the legendre polynomial of degree n"""
    return [legendre_root(n, k) for k in range(1, n+1)]
