def box(f, a, b):
    return (b-a)*f(a)

def naiive_int(f, a, b, tol):
    h = b - a
    int1 = box(f, a, b)
    int2 = box(f, a, a+h/2) + box(f, a+h/2, b)
    while abs(int1 - int2) > tol:
        h = h/2
        int1 = int2
        int2 = sum([box(f, a+i*h, a+(i+1)*h) for i in range(int((b-a)/h) + 1)])
    return int2

print(naiive_int(lambda x: x**2, 0, 1, 0.0001))