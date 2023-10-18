import matplotlib.pyplot as plt

class Func:
    def __init__(self, f):
        self.f = f

    def __call__(self, x : float) -> float:
        return self.f(x)

    def div_diff(self, xs : list[float]) -> float:
        
        if len(xs) == 1:
            return self(xs[0])
        
        return (self.div_diff(xs[1:]) - self.div_diff(xs[:-1])) / (xs[-1] - xs[0])
    

    def plot(self, a : float, b : float):
        x = [a + (b - a)*i/1000 for i in range(1001)]
        y = [self(point) for point in x]
        plt.plot(x, y)
        plt.show()

    # def bisec(self, a : float, b : float):
    #     """returns a root of f in the interval (a, b) if it exists"""


    

    # def Der(self):

    # def Int(self):


# Func(lambda x : x**2 + 3).plot(-2, 5)

# print((Func(lambda x : x**2)).div_diff([1, 2, 4]))