#import sympy as sym
from sympy import *
import math
import sys
import os
import gravipy
import numpy

Mp = 1.0
'''
def kahler(sufi, p, S=1):
    sum = sufi[0] + conjugate(sufi[0])
    interm0 = Symbol('S*_0')
    interm1 = Symbol('S*_1')
    func = (
        sum
        - p[2]
        * ((sum - 1) * cos(p[1]) - 1j * (sufi[0] - conjugate(sufi[0])) * sin(p[1])) ** 4
        #* ((sum - 1) ) ** 4
    )
    if S == 0:
        superterm = 0
    else:
        superterm = conjugate(sufi[1])* sufi[1]/ (sum ** 3)

    potential = -3 * p[0] * Mp * Mp * ln(func) + superterm
    return potential.subs(conjugate(sufi[0]), interm0).subs(conjugate(sufi[1]), interm1)


def F(sufi, p):
    return (3.0 / 4) ** 0.5 * p[4] / p[3] * (sufi[0] - p[3])


def V_sympy(field, p):
    sufi = symarray("S", 2)
    interm = Symbol('S*_1')
    sufi[0] = field[0] + 1j * field[1]

    dKdSS = kahler(sufi, p).diff(sufi[1]).diff(interm)


    FF = F(sufi, p)*conjugate(F(sufi, p))
    print(FF)
    return exp(kahler(sufi, p, S=0) / (Mp * Mp)) * FF / dKdSS


def G_sympy(field, p=None):
    sufi = symarray("S", 2)
    interm = Symbol('S*_0')
    sum = field[0] + 1j * field[1]
    factor = 2 * kahler(sufi, p).diff(sufi[0]).diff(interm)#.simplify()

    factor = factor.subs(interm, conjugate(sum))
    factor = factor.subs(sufi[0], sum)#.simplify()
    return factor * Matrix([[1, 0], [0, 1]])
#'''

#'''
def V_sympy(field, p):
    upper = 6*p[1]**2*field[0]**3*((p[3]-field[0])**2+field[1]**2)
    lower = (2*field[0] - p[2]*(1-2*field[0])**4)**(3*p[0])*p[3]**2
    #upper = Mp**(3*p[0]) * (p[1]**2*field[1]**2 + (p[2]+p[1]*field[0])**2)
    #lower = (2*field[0])**(3*p[0])
    return upper/lower


def G_sympy(field, p=None):
    c = p[2]
    r = 2*field[0]
    F = (r-1)
    factor = 3*p[0]*(1 + 4*c*F**2 * (3*r - F*2 + c*F**3))/(r-c*F**4)**2
    #factor = 3*p[2]/(2*field[0**2])
    return 2*factor * Matrix([[1, 0], [0, 1]])
#'''

def setup():
    try:
        location = "/net/vdesk/data2/boxelaar/Inflation/PyTransport-master/PyTransport"
        sys.path.append(location)
        import PyTransSetup
    except:
        location = "/Users/jortboxelaar/Documents/Study/Masters/Project_2/Inflation/PyTransport-master/PyTransport"
        sys.path.append(location)
        import PyTransSetup
    import PyTransSetup

    nF = 2
    nP = 4
    f = symarray("f", nF, real=True)
    p = symarray("p", nP, real=True)

    #print(G_sympy(f,p))
    #print(V_sympy(f,p))
    #sys.exit()

    pot_name = str(os.getcwd().split("/")[-1])
    PyTransSetup.potential(V_sympy(f, p), nF, nP, False, G_sympy(f, p))
    PyTransSetup.compileName3(pot_name, True)

    os.system(
        f"cp {location}/PyTrans/lib/python/PyTrans{pot_name}.cpython-37m-darwin.so\
         {location}/PyTrans/lib64/python/PyTrans{pot_name}.cpython-37m-darwin.so"
    )
    os.system(
        f"cp {location}/PyTrans/lib/python/PyTrans{pot_name}-1.0-py3.7.egg-info\
         {location}/PyTrans/lib64/python/PyTrans{pot_name}-1.0-py3.7.egg-info"
    )


if __name__ == "__main__":
    setup()
