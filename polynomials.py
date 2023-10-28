from math import *
from functions import *
import typing
import matplotlib.pyplot as plt
from roots import *
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
    return Poly([-f.coefficients[i] + g.coefficients[i] if i <= f.degree else -g.coefficients[i] for i in range(g.degree + 1)])

def Poly_mul(f : Poly, g : Poly) -> Poly:
    """returns the product of two polynomials"""
    coefficients = [0 for _ in range(f.degree + g.degree + 1)]
    for i in range(f.degree + 1):
        for j in range(g.degree + 1):
            coefficients[i + j] += f.coefficients[i]*g.coefficients[j]
    return Poly(coefficients)

def Poly_div(f : Poly, g : Poly) -> Poly:
    """returns the quotient of when f is divided by g"""
    if f.degree < g.degree:
        return Poly([0])
    
    h = f.copy()
    q = Poly([0])

    while h.degree >= g.degree:
        q = Poly_add(q, Poly([h.coefficients[0]/g.coefficients[0]] + [0 for _ in range(h.degree - g.degree)]))
        h = Poly_sub(h, Poly_mul(q, g))

    return q

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

    

def laggy_interpol(points : list) -> Poly:
    """uses lagrange interpolation to return an interpolating polynomial whose graph passes through all the points in the list"""
    
    n = len(points)

    result = Poly([0])

    for i in range(len(points)):
        yL_nk = Poly([points[i][1]])

        x_i = points[i][0]

        for j in range(i):
            yL_nk = Poly_mul(yL_nk, Poly([-points[j][0]/(x_i - points[j][0]), 1/(x_i - points[j][0])]))

        for j in range(i+1, n):
            yL_nk = Poly_mul(yL_nk, Poly([-points[j][0]/(x_i - points[j][0]), 1/(x_i - points[j][0])]))

        result = Poly_add(result, yL_nk)

    return result


def interpol(f : Func, points : list[float]):
    """uses newton interpolation to return an n-degree polynomial which interpolates f at all the values in points"""

    F = Func(f)
    result = Poly([f(points[0])])
    P = Poly([1])

    for k in range(1, len(points)):
        P = Poly_mul(P, Poly([-points[k-1], 1]))
        result = Poly_add(result, Poly_mul(Poly([F.div_diff(points[:k+1])]), P))

    return result

    

def newt_interpol(f, a, b, n):
    """returns an n-degree polynomial which interpolates f at n + 1 equispaced points starting at a and ending at b"""

    if n == 0:
        return Poly([f(a)])

    h = (b - a)/n

    return interpol(f, [a + i*h for i in range(n+1)])


def legendre_poly(n):
    """returns the nth legendre polynomial in monic form (coefficient of x^n is 1)"""
    if n == 0:
        return Poly([1])
    if n == 1:
        return Poly([0, 1])
    
    return Poly_div(Poly_mul((Poly_mul(Poly([2*n - 1]), Poly([0, 1])), legendre_poly(n-1)) - Poly_mul(Poly([n - 1]), legendre_poly(n-2))), Poly([n]))

def legendre_root(n, k):
    """returns the kth root of the legendre polynomial of degree n"""
    if not n:
        if k != 1:
            return None
        return 0
    if k > n:
        return None
    e_nk = (1 - 1/(8*n**2) + 1/(8*n**3))*cos(pi*(4*k - 1)/(4*n + 2))

    return newton(legendre_poly(n), e_nk, 0.00001, 10)

def legendre_roots(n):
    """returns the roots of the legendre polynomial of degree n"""
    return [legendre_root(n, k) for k in range(1, n+1)]

# print(legendre_poly(3))
# print(legendre_roots(3))
