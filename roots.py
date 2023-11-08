import math
from scipy.misc import derivative
from linear_algebra import *



def bisection(f, a, b, tol, N):
    a1 = a
    b1 = b
    i = 0
    E = (b1 - a1)/2
    while E > tol and i < N:
        p =  (a1 + b1)/2
        if f(p) == 0:
            return (p, i)
        elif f(a1) < 0 and f(p) < 0 or f(a1) > 0 and f(p) > 0:
            a1 = p
        else:
            b1 = p
        E = (b1 - a1)/2
        i += 1
    if E < tol:
        return (p, i)
    return "method failed"

def fixed_point(f, p0, tol, N):
    i = 0
    p = p0
    while i < N:
        p1 = f(p)
        if abs(p1 - p) < tol:
            return p1
        p = f(p)
        i += 1
    return p1

def newton(f, p0, tol, N):
    i = 0
    p = p0

    df = derivative(f, p)

    while i < N:
        p1 = p - f(p)/df(p)
        if abs(p1 - p) < tol:
            return (p1, i)
        p = p1
        i += 1
    return "method failed"

def secant(f, p0, p1, tol, N):
    i = 0
    q0 = f(p0)
    q1 = f(p1)

    while i < N:
        p = p1 - q1*(p1 - p0)/(q1 - q0)
        if abs(p - p1) < tol:
            return (p, i)
        p0 = p1
        q0 = q1
        p1 = p
        q1 = f(p)
        i += 1
    return "method failed"

def regula_falsi(f, a, b, tol, N):
    i = 0
    fa = f(a)
    fb = f(b)
    p = (a*fb - b*fa)/(fb - fa)

    while i < N:
        fp = f(p)
        if abs(fp) < tol:
            return (p, i)
        if fa*fp < 0:
            b = p
            fb = fp
        else:
            a = p
            fa = fp
        p = (a*fb - b*fa)/(fb - fa)
        i += 1
    return "method failed"

