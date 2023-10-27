import math
import functions
import typing
import matplotlib.pyplot as plt

class Poly(functions.Func):
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


def interpol(f : functions, points : list[float]):
    """uses newton interpolation to return an n-degree polynomial which interpolates f at all the values in points"""

    F = functions.Func(f)
    result = Poly([f(points[0])])
    P = Poly([1])

    for k in range(1, len(points)):
        P = Poly_mul(P, Poly([-points[k-1], 1]))
        result = Poly_add(result, Poly_mul(Poly([F.div_diff(points[:k+1])]), P))

    return result

    

def newt_interpol(f, a, b, n):
    """returns an n-degree polynomial which interpolates f at n + 1 equispaced points starting at a and ending at b"""

    h = (b - a)/n

    return interpol(f, [a + i*h for i in range(n+1)])





(newt_interpol(lambda x: x**2, 0, 10, 10)).print()
(newt_interpol(lambda x: x**2, 0, 10, 10)).plot(0, 10)
