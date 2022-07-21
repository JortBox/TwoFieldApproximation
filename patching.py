#!/usr/bin/env python
# coding: utf-8

#!/usr/bin/env python
# coding: utf-8

import time
import math
import sys
import os
import pickle

import sympy as sym
import numpy as np
import matplotlib.pyplot as plt

from pylab import *
from matplotlib import cm
from matplotlib import style
from matplotlib.ticker import LinearLocator

from potentials import transutils as utils
#from potentials import orbital as pot
from potentials.supergravity import setup_supergravity
#from potentials.orbital import setup_orbital

from tqdm.contrib.telegram import tqdm, trange


Mp = 1

# plt.style.use('classic')
plt.rcParams.update(
    {
        "text.usetex": True,
        "font.family": "serif",
        "font.serif": ["Computer Modern Roman"],
    }
)
plt.rcParams.update({"axes.labelsize": 15, "axes.titlesize": 25, "legend.fontsize": 10})



def run_bispectrum(bispec):
    print("start bispectrum")
    bispec()
    bispec.plot_fnl(save=True, show=False)
    bispec.plot_bispectrum(save=True, show=False)
    bispec.plot_three_pt_correlations(save=True, show=False)
    return bispec


def run_powerspectrum(twopt):
    print("start spectrum")
    twopt()
    twopt.plot_sigma_evolution(save=True, show=False)
    twopt.plot_zeta_evolution(save=True, show=False)
    #pz, k_out = twopt.p_spectra(pool=True)
    #twopt.plot_p_spectra(pz, k_out, save=True, show=False)

    twopt.spectral_index()
    print(twopt.n_s, twopt.r)
    return twopt

# EGNO
alpha = 1
M = 1e-3
c = 1e3
a = 0.5


def evolution(N_start, N_end, n_steps, R0=1, A0=100, V0=1):
    fields = np.array([R0, A0], dtype=np.float64)
    p = np.array([V0, R0, 0, 0], dtype=np.float64)
    print(f"python3 orbital_tr_analyser.py orbital -i {N_start} -f {N_end} -steps {n_steps} -R {R0} -V {V0} -A {A0}")
    os.system(f"python3 orbital_tr_analyser.py orbital -i {N_start} -f {N_end} -steps {n_steps} -R {R0} -V {V0} -A {A0}")

    with open(f"./data/orbital.obj", "rb") as file:
        orb = pickle.load(file)
    return orb


def fit_trajectory(sfield, N_start, N_end, steps, A0=100):
    patch = (sfield.efolds >= N_start)&(sfield.efolds < N_end)

    kappa_EGNO = sfield.radius_of_curvature[patch][len(patch)//2]
    epsilon_EGNO_half = sfield.slow_roll[patch][len(patch)//2]
    epsilon_EGNO_init = sfield.slow_roll[patch][0]
    mass_EGNO = sfield.entropy_mass[patch][len(patch)//2]

    '''
    This function needs to have the implementation of the equations of the main text.
    Now, only the radius of curvature is fitted as a test. Also, the coordinate
    transformation is not properly implemented yet.
    '''
    R0 = float(
        sfield.value[sfield.efolds==N_end,0]
        + sfield.value[sfield.efolds==N_start,0]
    )/2
    print('R0',R0)
    orb = evolution(N_start, N_end, steps, R0, A0=A0)
    return orb

def patching(sfield, patch_init=6, treshold=3e-5):
    patches = [sfield.efolds[sfield.efolds>=patch_init][0]]
    turningrate = sfield.radius_of_curvature
    reference = turningrate[sfield.efolds>=patches[-1]][0]
    sel = (sfield.efolds>=patches[-1])&(sfield.slow_roll<1)
    patch_ref = np.zeros(turningrate.shape)

    average = np.mean(turningrate[sel])
    for i, rate in enumerate(turningrate[sel]):
        if abs(rate - reference)/average > treshold:
            patches.append(sfield.efolds[sel][i])
            reference = rate

        patch_ref[len(patch_ref[sfield.efolds<=patches[0]])+i] += float(len(patches))
    patches.append(sfield.efolds[sel][-1])
    return patches, patch_ref






s = open('./data/supergravity.obj', 'rb')
sugra = pickle.load(s)
sfield = sugra.fields
patches, patch_ref = patching(sfield, treshold=1e-1)
print(patches)

orb_fields = list()
back=None

for i in range(1,len(patches)):
    if i==1:
        A0 = 100
    else:
        A0 = orb_fields[-1].value[-1,1]

    cosmo = fit_trajectory(
        sfield,
        float(patches[i-1]),
        float(patches[i]),
        int(np.shape(np.where(patch_ref==i)[0])[0]),
        A0
    )
    if back is None:
        back = cosmo
    else:
        back = back + cosmo

    fields = cosmo.fields

with open(f"./data/patching.obj", "wb") as file:
    pickle.load(back, file)
