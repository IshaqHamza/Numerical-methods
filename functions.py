import matplotlib.pyplot as plt
from differentiation import *
from roots import *
from polynomials import *
class Func:
    def __init__(self, f):
        self.f = f

    def __call__(self, x : float) -> float:
        return self.f(x)

    def div_diff(self, xs : list[float]) -> float:
        
        if len(xs) == 1:
            return self(xs[0])
        
        return (self.div_diff(xs[1:]) - self.div_diff(xs[:-1])) / (xs[-1] - xs[0])
    
    def copy(self):
        return Func(self.f)
    

    def plot(self, a : float, b : float):
        x = [a + (b - a)*i/1000 for i in range(1001)]
        y = [self(point) for point in x]
        plt.plot(x, y)
        plt.show()


    # def bisec(self, a : float, b : float):
    #     """returns a root of f in the interval (a, b) if it exists"""


    

    def Der(self):
        return Func(lambda x : central_der(self.f, x, 0.00001, 100))
    
    def mult_Der(self, n, h = 0.00001, N = 100):
        if n == 0:
            return self
        return Func(lambda x: mult_der(self.f, x, n, h, N))

    def root(self, a : float, b : float, tol = 0.00001, N = 10000):
        return regula_falsi(self.f, a, b, tol, N)
    

    # def Int(self):


# Func(lambda x : x**2 + 3).plot(-2, 5)

# print((Func(lambda x : x**2)).div_diff([1, 2, 4]))