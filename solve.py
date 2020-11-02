# !/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Author:  Qingchun Wang @ NJU
# E-Mail:  qingchun720@foxmail.com
#

from sympy import *
from scipy import optimize

o1,o2 = 0.5*pi,1.5*pi
e1,e2,e3 = 0.0169,3.34,81.5

n = Symbol('n', real=True)

A = Matrix([
[   cos(n*o1),     sin(n*o1),    -cos(n*o1),    -sin(n*o1),           0,                  0],
[e1*sin(n*o1), -e1*cos(n*o1), -e2*sin(n*o1),  e2*cos(n*o1),           0,                  0],
[           0,             0,     cos(n*o2),     sin(n*o2),     -cos(n*o2),      -sin(n*o2)],
[           0,             0,  e2*sin(n*o2), -e2*cos(n*o2),  -e3*sin(n*o2),    e3*cos(n*o2)],
[          -1,             0,             0,             0,    cos(2*n*pi),     sin(2*n*pi)],
[           0,            e1,             0,             0, e3*sin(2*n*pi), -e3*cos(2*n*pi)]])

fn = A.det()
# fn = -8*e1*e2*e3 \
#      +(e1+e2)*(e2+e3)*(e3+e1)*cos(2*n*pi) \
#      -(e1+e2)*(e2-e3)*(e3-e1)*cos(2*n*o1) \
#      +(e1-e2)*(e2+e3)*(e3-e1)*cos(2*n*o1) \
#      +(e1-e2)*(e2-e3)*(e3+e1)

plotting.plot(fn, (n, -pi, pi))
# fn.evalf(subs={n: 0})


def bisect(f, x, x0, x1):
    '''
    bisect for solve equation
    
    :param f:
    :param x:
    :param x0:
    :param x1:
    :return:
    '''
    
    while True:
        x3 = (x0+x1)/2
        f3 = f.evalf(subs={x: x3})
        if abs(f3) < 1.0e-12: break
        else:
            f0 = f.evalf(subs={x: x0})
            f1 = f.evalf(subs={x: x1})
            if f0*f3<0: x0,x1 = x0,x3
            elif f1*f3: x0,x1 = x3,x1
            else:
                raise RuntimeError(F'f({x0}), f({x3}), f({x1}) are same sign')

    return x3


x = Symbol('x')
fx = 2*x**2-1
plotting.plot(fx, (x, -1, 1))
# x,*_ = solve(fx, x)
# x = bisect(fx, x, 0.5, 1)
fxpy = lambdify(x, fx)
# x = optimize.brentq(fxpy, 0.5, 1)
# x = optimize.newton(fxpy, 0.6)
x,*_ = optimize.fsolve(fxpy, 0.6)
print(F'x = {x}')

# n,*_ = solve(fn,n)
# n = bisect(fn, n, 0.6, 1.0)
fnpy = lambdify(n, fn)
# n = optimize.brentq(fnpy, 0.5, 1)
# n = optimize.newton(fnpy, 0.7)
n,*_ = optimize.fsolve(fnpy, 0.7)
print(F'n = {n}')


