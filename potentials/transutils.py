#!/usr/bin/env python
# coding: utf-8

import sys
import os
from gravipy.tensorial import *
from sympy import init_printing
import sympy as sym

import numpy as np
import matplotlib.pyplot as plt

#from tqdm.contrib.telegram import tqdm, trange
from scipy.interpolate import Akima1DInterpolator as AkimaSpline
from scipy.interpolate import CubicSpline

import matplotlib.colors as colors
import matplotlib.cm as cmx

def set_functions(potential):
    location = f"{os.getcwd()}/potentials/{potential}"
    sys.path.append(location)
    setup = __import__(f"setup_{potential}")
    return setup.G_sympy, setup.V_sympy



def linecolors(nvals, mapname):
    map = cm = plt.get_cmap(mapname)
    cNorm  = colors.Normalize(vmin=0, vmax=nvals-1)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=map)
    line_colors = []
    for i in range(nvals):
        line_colors.append(scalarMap.to_rgba(i))
    return line_colors

def det(mat):
    """Square root of determinant of metric"""
    return np.sqrt(np.linalg.det(mat))


def V_NN(fields, dV, ddV, time=0):
    r, t = symbols('f1, f2')
    x = Coordinates('chi', [r,t])
    G = MetricTensor('G', x, diag(fields.G_sympy([r,t],fields.params)))
    Ga = Christoffel('Ga', G)
    #print(Ga(All,All,All))

    sum = 0
    for a in range(fields.nF):
        for b in range(fields.nF):
            gamma_sum = 0
            for c in range(fields.nF):
                print('(',time,a,b,c,')')
                gamma_sum += Ga(-(c+1),a+1,b+1)*dV[time,c]
            gamma_sum = sym.lambdify([r, t], gamma_sum)
            sum += fields.N[time,a]*fields.N[time,b]*(ddV[time,a,b]-gamma_sum(*fields.value[time]))
    return sum


def ricci(fields, time=0):
    r, t = symbols('f1, f2')
    x = Coordinates('chi', [r,t])
    G = MetricTensor('G', x, diag(fields.G_sympy([r,t],fields.params)))
    Ga = Christoffel('Ga', G)
    Ri = Ricci('Ri', G).scalar()
    Ri = sym.lambdify([r, t], Ri)
    return Ri(*fields.value[time])





class Cubic_spline(object):
    def __init__(self, x, y):
        self.derivative()
        self.N = len(x)
        self.coeff = self.spline(x, y)
        self.x = x.copy()

    def __call__(self, x):
        m = self.order
        interpol = np.zeros(x.shape)
        for i in range(x.shape[0]):
            for j in range(self.N - 1):
                if x[i] >= self.x[j] and x[i] <= self.x[j + 1]:
                    interpol[i] = (
                        self.diff_factor(0)
                        * self.coeff[0][j]
                        * (x[i] - self.x[j]) ** (0 - m)
                        + self.diff_factor(1)
                        * self.coeff[1][j]
                        * (x[i] - self.x[j]) ** (1 - m)
                        + self.diff_factor(2)
                        * self.coeff[2][j]
                        * (x[i] - self.x[j]) ** (2 - m)
                        + self.diff_factor(3)
                        * self.coeff[3][j]
                        * (x[i] - self.x[j]) ** (3 - m)
                    )
        return interpol

    def derivative(self, order=0):
        self.order = order

    def spline(self, x, y):
        a = y.copy()
        b = np.zeros(self.N - 1)
        d = np.zeros(self.N - 1)
        c = np.zeros(self.N)

        l = c.copy()
        mu = c.copy()
        z = c.copy()
        h = x[1:] - x[:-1]

        for i in range(1, self.N - 1):
            l = 2 * (x[i + 1] - x[i - 1]) - h[i - 1] * mu[i - 1]
            mu[i] = h[i] / l
            z[i] = (
                3 * ((a[i + 1] - a[i]) / h[i] - (a[i] - a[i - 1]) / h[i - 1])
                - h[i - 1] * z[i - 1]
            ) / l

        for i in reversed(range(self.N - 1)):
            c[i] = z[i] - mu[i] * c[i + 1]
            b[i] = (a[i + 1] - a[i]) / h[i] - h[i] * (c[i + 1] + 2 * c[i]) / 3.0
            d[i] = (c[i + 1] - c[i]) / (3 * h[i])
        return a[:-1], b, c[:-1], d

    def diff_factor(self, exponent):
        factor = 1.0
        for i in range(self.order):
            factor *= exponent - i
        return factor


def eff_V(field, p):
    return (
        p[0] * field[:,1]
        * (
            1.
            + p[2] / 2.0 * (field[:,0] - p[1]) ** 2
            + p[3] / 6.0 * (field[:,0] - p[1]) ** 3
        )
    )



"""
# test data
x = np.array([1.0000,4.3333,7.6667,11.000,14.333,17.667,21.000])
y = np.array([1.4925,15.323,3.2356,-29.472,-22.396,24.019,36.863])

x_interp = np.linspace(x.min(), x.max(), 1000)
plt.scatter(x,y)

spline = CubicSpline(x,y)
interp = spline(x_interp)
plt.plot(x_interp, interp)

spline_diff = spline.derivative(1)
interp = spline_diff(x_interp)
plt.plot(x_interp, interp)

spline = Cubic_spline(x,y)
interp = spline(x_interp)
plt.plot(x_interp, interp)

spline.derivative(2)
interp = spline(x_interp)
plt.plot(x_interp, interp)

plt.show()
#"""
