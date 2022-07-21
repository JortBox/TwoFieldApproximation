from potentials.classes import *
from potentials.classes import set_potential

import sympy as sym
import math
import sys
import os
import gravipy

from .setup import *


path = os.path.dirname(os.path.realpath(__file__))
try:
    name = str(path.split('/')[-1])
    set_potential(name)
except:
    print(str(path.split('/')[-1]))
    setup(str(path.split('/')[-1]))
