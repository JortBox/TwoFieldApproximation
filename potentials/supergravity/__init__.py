from potentials.classes import *
from potentials.classes import set_potential

import sympy as sym
import math
import sys
import os
import gravipy

from .setup_supergravity import *

path = os.path.dirname(os.path.realpath(__file__))
name = str(path.split('/')[-1])

try:
    set_potential(name)
except:
    print(name)
    setup(name)
