import sympy as sym
import math
import sys
import os
import gravipy
import numpy

Mp = 1.0


def V_sympy(field, p):
    upper = 6*p[1]**2*field[0]**3*((p[3]-field[0])**2+field[1]**2)
    lower = (2*field[0] - p[2]*(1-2*field[0])**4)**(3*p[0])*p[3]**2
    return upper/lower


def G_sympy(field, p=None):
    c = p[2]
    r = 2*field[0]
    F = (r-1)
    factor = 3*p[0]*(1 + 4*c*F**2 * (3*r - F*2 + c*F**3))/(r-c*F**4)**2
    return 2*factor * sym.Matrix([[1, 0], [0, 1]])


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
    f = sym.symarray("f", nF)
    p = sym.symarray("p", nP)

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
