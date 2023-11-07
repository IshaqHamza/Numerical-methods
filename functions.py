import matplotlib.pyplot as plt
from roots import *
from polynomials import *
from scipy.misc import *
from math import *

class Func:
    def __init__(self, f):
        self.f = f

    def __call__(self, x : float) -> float:
        return self.f(x)

    def div_diff(self, xs : list[float]) -> float:
        
        if len(xs) == 1:
            return self(xs[0])
        
        return (self.div_diff(xs[1:]) - self.div_diff(xs[:-1])) / (xs[-1] - xs[0])

    def forward_diff_undiv(self, x : float, h : float, k) -> float:
        if k == 0:
            return self(x)
        return self.forward_diff_undiv(x + h, h, k-1) - self.forward_diff_undiv(x, h, k-1)

    def backward_diff_undiv(self, x : float, h : float, k) -> float:
        if k == 0:
            return self(x)
        return self.backward_diff_undiv(x, h, k-1) - self.backward_diff_undiv(x - h, h, k-1)

    def __add__(self, other):
        return Func(lambda x : self(x) + other(x))
    
    def copy(self):
        return Func(self.f)
    

    def plot(self, a : float, b : float):
        x = [a + (b - a)*i/1000 for i in range(1001)]
        y = [self(point) for point in x]
        plt.plot(x, y)
        plt.show()


    def Der(self):
        return Func(lambda x: derivative(self.f, x))
    
    def mult_Der(self, n,):
        if n == 0:
            return self
        return Func(lambda x : derivative(self.f, x, dx = 0.000001, n = n, order = n+1))

    def root(self, a : float, b : float, tol = 0.00001, N = 10000):
        return regula_falsi(self.f, a, b, tol, N)


    


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


def interpol(f, points : list[float]):
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