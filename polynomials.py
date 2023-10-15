import math

import matplotlib.pyplot as plt

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
    

def Poly_add(f, g):
    """returns the sum of two polynomials"""
    if f.degree > g.degree:
        return Poly([f.coefficients[i] + g.coefficients[i] if i <= g.degree else f.coefficients[i] for i in range(f.degree + 1)])
    return Poly([f.coefficients[i] + g.coefficients[i] if i <= f.degree else g.coefficients[i] for i in range(g.degree + 1)])

def Poly_sub(f, g):
    """returns the difference of two polynomials"""
    if f.degree > g.degree:
        return Poly([f.coefficients[i] - g.coefficients[i] if i <= g.degree else f.coefficients[i] for i in range(f.degree + 1)])
    return Poly([-f.coefficients[i] + g.coefficients[i] if i <= f.degree else -g.coefficients[i] for i in range(g.degree + 1)])

def Poly_mul(f, g):
    """returns the product of two polynomials"""
    coefficients = [0 for _ in range(f.degree + g.degree + 1)]
    for i in range(f.degree + 1):
        for j in range(g.degree + 1):
            coefficients[i + j] += f.coefficients[i]*g.coefficients[j]
    return Poly(coefficients)

