#!/usr/bin/env python
# coding: utf-8

import time
import sys
import os
import pickle
import dill

dill.settings["recurse"] = True

import sympy as sym
import numpy as np
import matplotlib.pyplot as plt

from pylab import *
from matplotlib import cm
from matplotlib import style
from matplotlib.ticker import LinearLocator

from potentials import transutils as utils

# from potentials import supergravity as pot

Mp = 1

import argparse

parser = argparse.ArgumentParser(description="set potential")
parser.add_argument("potential", help="potential name", type=str)
parser.add_argument(
    "-i", help="(float) initial e-folding. Default: 0", default=0, type=float
)
parser.add_argument(
    "-f", help="(float) final e-folding. Default: 60", default=60, type=float
)
parser.add_argument(
    "-steps", help="(int) number of steps. Default: 8000", default=8000, type=int
)
parser.add_argument(
    "-twopt",
    help="(bool) whether to perform twopt function. Default: False",
    default=False,
    type=bool,
)
parser.add_argument(
    "-bispec",
    help="(bool) whether to perform bispec function. Default: False",
    default=False,
    type=bool,
)

parser.add_argument(
    "-f1", help="(float) initial value field 1. Default: 1", default=1.0, type=float
)
parser.add_argument(
    "-f2", help="(float) initial value field 2. Default: 10", default=10.0, type=float
)
parser.add_argument(
    "-V",
    help="(float) V parameter Orbital Inflation. Default: 1",
    default=1.0,
    type=float,
)
parser.add_argument(
    "-lab",
    help="(float) lambda parameter Orbital Inflation. Default: 0",
    default=0,
    type=float,
)
parser.add_argument(
    "-beta",
    help="(float) beta parameter Orbital Inflation. Default: 0",
    default=0,
    type=float,
)

parser.add_argument(
    "-alp",
    help="(float) alpha parameter supergravity Inflation. Default: 1",
    default=1.0,
    type=float,
)
parser.add_argument(
    "-M",
    help="(float) M parameter supergravity Inflation. Default: 10",
    default=1.0e-3,
    type=float,
)
parser.add_argument(
    "-c",
    help="(float) c parameter supergravity Inflation. Default: 1",
    default=1.0e3,
    type=float,
)
parser.add_argument(
    "-a",
    help="(float)  a parameter supergravity Inflation.Default: 0",
    default=0.5,
    type=float,
)

args = parser.parse_args()


location = f"{os.getcwd()}/potentials/"
sys.path.append(location)
pot = __import__(args.potential)

plt.rcParams.update(
    {
        "text.usetex": True,
        "font.family": "serif",
        "font.serif": ["Computer Modern Roman"],
    }
)
plt.rcParams.update({"axes.labelsize": 15, "axes.titlesize": 25, "legend.fontsize": 10})


def main():
    if pot.name == "orbital":
        fields = np.array([args.f1, args.f2])  # 2.4*np.pi
        p = np.asarray([args.V, args.f1, args.lab, args.alp])

        fsym = sym.symarray("f", len(fields))
        psym = sym.symarray("p", len(p))
        dHdf = sym.lambdify([fsym, psym], pot.eff_V(fsym, psym).diff(fsym[1]))

        Cosmo = pot.BackgroundEvolution(fields, p)
        Cosmo.initial = np.concatenate(
            (fields, [0.0, -2.0 * dHdf(fields, p) / pot.func(fields, p)])
        )
        Cosmo(N_start=args.i, N_end=args.f, n_steps=args.steps, slow_roll=True)

    elif pot.name == "supergravity":
        if args.f1 == 1.0 or args.f2 == 10.0:
            fields = np.array([a, 3.7 * a * np.sqrt(2.0 / 3)])
        else:
            fields = np.array([args.f1, args.f2])
        p = np.asarray([args.alp, args.M, args.c, args.a])

        Cosmo = pot.BackgroundEvolution(fields, p)
        vel = np.array([-0.22, -0.67]) * 1e-3
        Cosmo.initial = np.concatenate((fields, vel))
        # Cosmo.transform_coord()
        Cosmo(N_start=args.i, N_end=args.f, n_steps=args.steps, slow_roll=False)

    with open(f"./data/{pot.name}.obj", "wb") as file:
        pickle.dump(Cosmo, file)

    Cosmo.plot_background(save=False, show=True)
    Cosmo.plot_turningrate(save=False)
    pot.plots_for_background(Cosmo, save=False, show=True, to_xy=True)

    if args.twopt:
        twopt = pot.PowerSpectrum(Cosmo, NB=4, N_exit=5)
        twopt()
        twopt.plot_sigma_evolution(save=True, show=False)
        twopt.plot_zeta_evolution(save=True, show=False)
        pz, k_out = twopt.p_spectra(pool=True)
        twopt.plot_p_spectra(pz, k_out, save=True, show=False)

        with open(f"./data/{pot.name}_twopt.obj", "wb") as file:
            pickle.dump(twopt, file)

    if args.bispec:
        bispec = pot.Bispectrum(Cosmo, alpha=0.0, beta=1.0 / 3, NB=6, N_exit=10)
        bispec()
        bispec.plot_fnl(save=True, show=False)
        bispec.plot_bispectrum(save=True, show=False)
        bispec.plot_three_pt_correlations(save=True, show=False)

        with open(f"./data/{pot.name}_bispec.obj", "wb") as file:
            pickle.dump(bispec, file)


if __name__ == "__main__":
    sys.exit(main())
