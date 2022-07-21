import sympy as sym
import math
import sys
import os
import gravipy
import numpy

Mp = 1.0


def func(field, p=None):
    return field[0] ** 2
    #return sym.exp(2*field[0]/p[1])


def eff_V(field, p):
    return (
        p[0] * field[1]
        * (
            1.
            + p[2] / 2.0 * (field[0] - p[1]) ** 2
            + p[3] / 6.0 * (field[0] - p[1]) ** 3
        )
    )


def V_sympy(field, p):
    H = (
        p[0] * field[1]
        * (
            1.
            + p[2] / math.factorial(2) * (field[0] - p[1]) ** 2
            + p[3] / math.factorial(3) * (field[0] - p[1]) ** 3
            #+ p[4] / math.factorial(4) * (field[0] - p[1]) ** 4
            #+ p[5] / math.factorial(5) * (field[0] - p[1]) ** 5
        )
    )

    H_rho = (
        p[0] * field[1]
        * (
            p[2] * (field[0] - p[1])
            + p[3] / 2 * (field[0] - p[1]) ** 2
        )
    )
    return 3 * H ** 2 - 2 * ((H / field[1]) ** 2 / func(field) +  H_rho**2)

def G_sympy(field, p=None):
    if p is None:
        p = sym.symarray("p", nP)
    return sym.Matrix([[1, 0], [0, sym.exp(2*field[0]/p[1])]])


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
