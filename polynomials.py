import math

import matplotlib.pyplot as plt

# def poly_add(poly1, poly2):
#     """returns the sum of two polynomials"""
#     if len(poly1) < len(poly2):
#         poly1, poly2 = poly2, poly1
#     poly3 = poly1.copy()
#     for i in range(len(poly2)):
#         poly3[i] += poly2[i]
#     return poly3

# def poly_sub(poly1, poly2):
#     """returns the difference of two polynomials"""
#     if len(poly1) < len(poly2):
#         poly1, poly2 = poly2, poly1
#     poly3 = poly1.copy()
#     for i in range(len(poly2)):
#         poly3[i] -= poly2[i]
#     return poly3

# def poly_mult(poly1, poly2):
#     """returns the product of two polynomials"""
#     poly3 = [0]*(len(poly1) + len(poly2) - 1)
#     for i in range(len(poly1)):
#         for j in range(len(poly2)):
#             poly3[i + j] += poly1[i]*poly2[j]
#     return poly3

# def poly_div(poly1, poly2):
#     """returns the quotient of two polynomials"""
#     poly3 = poly1.copy()
#     poly4 = [0]*(len(poly1) - len(poly2) + 1)
#     for i in range(len(poly1) - len(poly2) + 1):
#         poly4[i] = poly3[0]/poly2[0]
#         for j in range(len(poly2)):
#             poly3[i + j] -= poly4[i]*poly2[j]
#         poly3.pop(0)
#     return poly4

# def poly_rem(poly1, poly2):
#     """returns the remainder of two polynomials"""
#     poly3 = poly1.copy()
#     poly4 = [0]*(len(poly1) - len(poly2) + 1)
#     for i in range(len(poly1) - len(poly2) + 1):
#         poly4[i] = poly3[0]/poly2[0]
#         for j in range(len(poly2)):
#             poly3[i + j] -= poly4[i]*poly2[j]
#         poly3.pop(0)
#     return poly3

# def poly_eval(poly, x):
#     """returns the value of a polynomial at a point"""
#     value = 0
#     for i in range(len(poly)):
#         value += poly[i]*(x**i)
#     return value

# def poly_deriv(poly):
#     """returns the derivative of a polynomial"""
#     poly2 = poly.copy()
#     for i in range(len(poly2)):
#         poly2[i] *= i
#     poly2.pop(0)
#     return poly2

# def poly_integ(poly):
#     """returns the integral of a polynomial"""
#     poly2 = poly.copy()
#     for i in range(len(poly2)):
#         poly2[i] /= (i + 1)
#     poly2.insert(0, 0)
#     return poly2

# def lagrange_poly(points):
#     """returns lagrange interpolating polynomial for a list of points in the form of a list of tuples"""
#     poly = [0]*len(points)
#     for i in range(len(points)):
#         poly2 = [points[i][1]]
#         for j in range(len(points)):
#             if i != j:
#                 poly2 = poly_mult(poly2, [-points[j][0], 1])
#                 poly2 = [term/(points[i][0] - points[j][0]) for term in poly2]
#         poly = poly_add(poly, poly2)
#     return poly

# def poly_approx_lagrange(f, a, b, n):
#     """returns the lagrange interpolating polynomial for a function f on the interval [a, b] with n points"""
#     points = [(a + (b - a)*i/n, f(a + (b - a)*i/n)) for i in range(n + 1)]
#     return lagrange_poly(points)

# def poly_approx_lagrange_pts(f, points: list):
#     return lagrange_poly([(point, f(point)) for point in points])

class Poly:
    def __init__(self, coefficients: list): # coefficients are in the form [a_0, a_1, ..., a_n] for a_n*x^n + ... + a_1*x + a_0
        self.coefficients = coefficients
        self.degree = len(coefficients) - 1

    def __call__(self, x):  
        """returns the value of a polynomial at a point x, using Horner's method"""
        value = 0
        for i in range(len(self.coefficients)-1, -1, -1):
            value = self.coefficients[i] + x*value
        return value
    
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
    
p = Poly([1, 2, 3])
print(p(2))
p.plot(-1, 1)