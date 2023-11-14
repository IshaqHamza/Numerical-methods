import math
import numpy as np
def linearfinitediff(a,b,y_a,y_b,n):
    h=(b-a)/(n+1)
    A=np.zeros((n,n))
    b=np.zeros(n)
    for i in range(n): #fill the diagonal
        x=a+(i+1)*h
        A[i][i]=2+h**2*q(x)
    for i in range(n-1): #fill the upper diagonal
        x=a+(i+1)*h
        A[i][i+1]=-1+h/2*p(x)
        A[i+1][i]=-1-h/2*p(x)
    print(A)
    for i in range(n): #fill the b vector
        x=a+(i+1)*h
        print(x)
        if i==0:
            b[i]=-h**2*r(x)+y_a*(1+h/2*p(x))
        elif i==n-1:
            b[i]=-h**2*r(x)+y_b*(1-h/2*p(x))
        else:
            b[i]=-h**2*r(x)
        print(b[i])
    # print(b[0])
    # print(b)
    y=np.linalg.solve(A,b)
    return y

def p(x):
    return -2/x

def q(x):
    return 2/x**2

def r(x):
    return math.sin(math.log(x))/x**2

def error(a,b,y_a,y_b,n):
    h=(b-a)/(n+1)
    y=linearfinitediff(a,b,y_a,y_b,n)
    print("approx:", y)
    error=[]
    for i in range(n):
        x=a+(i+1)*h
        error.append(abs(y[i]-u(x)))
        print(u(x))
    return error

def u(x):
    return 1.139*x - 0.039/x**2 -0.3*math.sin(math.log(x)) - 0.1*math.cos(math.log(x))
