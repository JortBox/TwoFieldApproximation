#!/usr/bin/env python
# coding: utf-8

import time
import math
import sys
import os
import pickle
import dill

dill.settings["recurse"] = True
from copy import deepcopy

import numpy as np
import sympy as sym
from gravipy.tensorial import *

from scipy import interpolate
import timeit
from functools import partial

from pylab import *
from matplotlib import pyplot as plt
from scipy.misc import derivative

import potentials.transutils as utils

Mp = 1.0  # TEST


def set_potential(potential):
    global PyT
    global PyS
    global potential_name
    global plotpath
    global potential_setup

    potential_name = copy(potential)
    location = "./PyTransport-master/PyTransport"
    sys.path.append(location)
    path = os.path.dirname(os.path.realpath(__file__))
    plotpath = "/".join(path.split("/")) + "/" + str(potential_name) + "/plots/"
    module = "PyTrans" + str(potential_name)

    if not os.path.isdir(plotpath):
        os.makedirs(plotpath)

    import PyTransSetup

    PyTransSetup.pathSet()
    PyT = __import__(module)

    import PyTransScripts as PyS

    location = f"{os.getcwd()}/potentials/{potential}"
    sys.path.append(location)
    potential_setup = __import__(f"setup_{potential}")


class Inflaton:
    def __init__(self, evolution):
        self.G_sympy = evolution.G_sympy
        self.V_sympy = evolution.V_sympy

        self.params = evolution.params.copy()
        self.nF = evolution.nF
        try:
            self.transformation = evolution.integral
        except:
            self.transformation = lambda f, p: f / f

        self.steps = evolution.background.shape[0]
        self.efolds = evolution.background[:, 0]
        self.value = evolution.background[:, 1 : 1 + self.nF]
        self.dot = evolution.background[:, -self.nF :]
        self.field = self.value.copy()

        self.T = self.tangent_vector()
        self.N = self.normal_vector()

    def __len__(self):
        return len(self.efolds)

    def abs(self, mat=None):
        if mat is None:
            mat = self.value
        return np.sqrt(self.prod(mat))

    @property
    def metric(self):
        metric = [self.G_sympy(val, self.params) for val in self.value]
        return np.array(metric, dtype=np.float64)

    @property
    def radius_of_curvature(self, mode="1"):
        return self.abs(self.dot) / self.abs(self.DTdt())

    @property
    def turning_rate(self):
        return self.abs(self.dot) / self.radius_of_curvature

    @property
    def slow_roll(self):
        return -self.dHdt() / (self.H() ** 2)

    @property
    def entropy_mass(self):
        try:
            self.V_NN *= 1
        except TypeError:
            self.V_NN()
        try:
            self.ricci *= 1
        except TypeError:
            self.ricci()

        return self.V_NN + self.H() ** 2 * self.slow_roll * (
            self.ricci + 6 / self.radius_of_curvature ** 2
        )

    def tangent_vector(self):
        tangent = self.dot / (np.tile(self.abs(self.dot), (2, 1)).T)
        self.T = tangent
        return tangent

    def normal_vector(self, mode="simple"):
        if mode == "general":
            dTdt = utils.Cubic_spline(self.efolds, self.T[:, 0])
            dTdt.derivative(order=1)
            DTdt = dTdt(self.efolds)
            # normal = DTdt/self.vector_norm(DTdt)

        elif mode == "simple":
            normal = np.zeros((self.steps, 2))
            norm = self.abs(self.dot)
            G = self.metric
            for i in range(self.steps):
                matrix = np.array(
                    [[-G[i, 0, 1], -G[i, 1, 1]], [G[i, 0, 0], G[i, 1, 0]]]
                )
                normal[i] = np.dot(matrix, self.dot[i]) / (norm[i] * utils.det(G[i]))

        self.N = normal
        return normal

    def prod(self, vector1, vector2=None, metric=None):
        if vector2 is None:
            vector2 = np.copy(vector1)
        if metric is None:
            metric = self.metric

        if len(vector1.shape) == 1:
            prod = np.dot(np.dot(metric[-1], vector1), vector2)
        else:
            prod = np.array(
                [
                    np.dot(np.dot(metric[i], vector1[i]), vector2[i])
                    for i in range(self.steps)
                ]
            )
        return prod

    def DTdt(self):
        diff = self.dVdf()
        V_N = np.array([np.dot(diff[i], self.N[i]) for i in range(self.steps)])
        return -self.N * np.tile(V_N / self.abs(self.dot), (2, 1)).T

    def H(self):
        print(self.V())
        print("value", self.value[:, 0])
        return np.sqrt((0.5 * self.prod(self.dot) + self.V()) / (3 * Mp))

    def dHdt(self):
        return -0.5 * self.prod(self.dot) / (Mp * Mp)

    def dHdN(self):
        return self.dHdt() / self.H()

    def dHdf(self):
        return np.array(
            [
                -np.dot(self.G_sympy(val, self.params), dotval)
                for val, dotval in zip(self.field, self.dot)
            ]
        ) / (2.0 * Mp ** 2)

    def V(self):
        # field = self.field.copy().reshape((-1, 2))
        potential = np.array(
            [self.V_sympy(val, self.params) for val in self.value], dtype=np.float64
        )
        return potential

    def dVdf(self):
        dvdf = np.zeros(self.field.shape)
        f = sym.symarray("f", self.nF)
        p = sym.symarray("p", len(self.params))

        for a in range(self.nF):
            differential = sym.lambdify([f, p], self.V_sympy(f, p).diff(f[a]))
            for i in range(len(dvdf)):
                dvdf[i, a] = differential(self.field[i], self.params)
        return dvdf

    def ddV(self):
        ddvdf = np.zeros((len(self), self.nF, self.nF))
        f = sym.symarray("f", self.nF)
        p = sym.symarray("p", len(self.params))

        for a in range(self.nF):
            for b in range(self.nF):
                differential = sym.lambdify(
                    [f, p], self.V_sympy(f, p).diff(f[a]).diff(f[b])
                )
                ddvdf[:, a, b] = np.array(
                    [differential(self.field[i], self.params) for i in range(len(self))]
                )

        return ddvdf

    def V_NN(self):
        V_NN = np.zeros(len(self))
        dV = self.dVdf()
        ddV = self.ddV()
        self.V_NN = np.array([utils.V_NN(self, dV, ddV, i) for i in range(len(self))])
        # self.V_NN = np.array([self.gravipy_V_NN(dV, ddV, i) for i in range(len(self))])

    def ricci(self):
        # self.ricci = np.array([self.gravipy_ricci(i) for i in range(len(self))])
        self.ricci = np.array([utils.ricci(self, i) for i in range(len(self))])

    """
    def gravipy_V_NN(self, dV, ddV, time=0):
        r, t = symbols('f1, f2')
        x = Coordinates('chi', [r,t])
        G = MetricTensor('G', x, diag(self.G_sympy([r,t],self.params)))
        Ga = Christoffel('Ga', G)

        sum = 0
        for a in range(self.nF):
            for b in range(self.nF):
                gamma_sum = 0
                for c in range(self.nF):
                    print('(',time,a,b,c,')')
                    gamma_sum += Ga(-(c+1),a+1,b+1)*dV[time,c]
                gamma_sum = sym.lambdify([r, t], gamma_sum)
                sum += self.N[time,a]*self.N[time,b]*(ddV[time,a,b]-gamma_sum(*self.value[time]))
        return sum


    def gravipy_ricci(self, time=0):
        r, t = symbols('f1, f2')
        x = Coordinates('chi', [r,t])
        G = MetricTensor('G', x, diag(self.G_sympy([r,t],self.params)))
        Ga = Christoffel('Ga', G)
        Ri = Ricci('Ri', G).scalar()
        Ri = sym.lambdify([r, t], Ri)
        scalar = Ri(*self.value[time])
        print(scalar)
    """


class BackgroundEvolution:
    def __init__(self, fields, params, names=["\\rho", "\\theta"]):
        fields = np.asarray(fields, dtype=np.float64)
        params = np.asarray(params, dtype=np.float64)
        self.params = self.set_p_value(params)
        self.nP = PyT.nP()
        self.nF = PyT.nF()
        self.names = names

        self.initial = np.concatenate(
            (
                fields,
                -1.0
                * PyT.dV(fields, self.params)
                / np.sqrt(3 * PyT.V(fields, self.params)),
            )
        )

        # print("V", PyT.V(fields, self.params))
        # print("dV", PyT.dV(fields, self.params))

        self.G_sympy = deepcopy(potential_setup.G_sympy)
        self.V_sympy = deepcopy(potential_setup.V_sympy)

        f = sym.symarray("f", 2)
        p = sym.symarray("p", len(self.params))
        # print(self.G_sympy(f,p))

    def __call__(self, N_start=0, N_end=60, n_steps=1000, slow_roll=True):
        """Background evolutuion for a number of e folds"""
        tols = np.array([10 ** -15, 10 ** -15])

        self.steps = n_steps
        self.efolds = np.linspace(N_start, N_end, n_steps)
        self.background = PyT.backEvolve(
            self.efolds, self.initial, self.params, tols, slow_roll
        )

        self.fields = Inflaton(self)

    def __len__(self):
        return self.steps

    def __add__(self, background):
        new = BackgroundEvolution(self.fields.value[0], self.params)
        new.steps = self.steps + background.steps
        new.initial = np.concatenate((self.fields.value[0], self.fields.dot[0]))
        new.efolds = np.concatenate((self.efolds, self.efolds))
        new.background = np.concatenate(
            (self.background, background.background), axis=0
        )
        new.fields = Inflaton(new)
        return new

    def set_p_value(self, params):
        # global p_value
        p_value = np.asarray(params.copy(), dtype=np.float64)
        return p_value

    def set_exit_mode(self, N_exit):
        self.N_exit = N_exit
        self.k = PyS.kexitN(N_exit, self.background, self.params, PyT)
        return self.k

    def transform_coord(self):
        """Function to transform radial coordinate from conformal
        to one of the form 1/2(r'^2 + g(r)t'^2)"""
        try:
            integral = dill.load(open("integral", "rb"))
            print("loaded coordinate transform from file")
        except:
            print("doing coordinate transformation using Sympy. This may take a while.")
            fsym = sym.symarray("f", self.nF)
            psym = sym.symarray("p", len(self.params))

            integral = sym.integrate(self.G_sympy(fsym, psym), fsym[0])
            dill.dump(sym.lambdify([fsym, psym], integral), open("integral", "wb"))

        self.integral = integral

    def plot_background(self, save=False, show=True):
        for i in range(self.nF):
            fig, ax = plt.subplots(
                2, 1, figsize=(8, 5), sharex=True, gridspec_kw={"height_ratios": [2, 1]}
            )
            fig.subplots_adjust(hspace=0.0)
            ax[0].plot(
                self.fields.efolds,
                self.fields.value[:, i],
                label=f"${self.names[i]}(N)$",
                color="black",
            )
            ax[1].plot(
                self.fields.efolds,
                self.fields.dot[:, i],
                label=f"${self.names[i]}'(N)$",
                color="black",
            )

            try:
                ax[0].axvline(self.fields.efolds[self.fields.slow_roll > 1][0])
                ax[1].axvline(self.fields.efolds[self.fields.slow_roll > 1][0])
            except:
                pass

            ax[1].set_xlabel("e-folds ($N$)")
            ax[0].set_ylabel(f"${self.names[i]}(N)$")
            ax[1].set_ylabel(f"${self.names[i]}'(N)$")
            ax[0].minorticks_on()
            ax[1].minorticks_on()
            ax[0].set_xlim(-0.75, self.fields.efolds.max())
            ax[1].set_xlim(-0.75, self.fields.efolds.max())

            if save:
                plt.savefig("%sbackground_%s_%s.pdf" % (plotpath, potential_name, i))
            if show:
                plt.show()
            plt.clf()
            plt.close()

    def plot_turningrate(self, save=False, show=True):
        fig, ax = plt.subplots()
        ax.plot(self.fields.efolds, self.fields.radius_of_curvature, color="black")
        ax.plot(self.fields.efolds, self.fields.value[:, 0], color="green")

        # ax[0].axvline(self.fields.efolds[self.slow_roll1>1][0])
        try:
            ax.axvline(self.fields.efolds[self.fields.slow_roll > 1][0])
        except:
            pass
        # ax[0].legend(frameon=False)
        # ax[1].legend(frameon=False)
        ax.set_xlabel("e-folds ($N$)")
        ax.set_ylabel(f"$\\Omega$")
        ax.minorticks_on()
        ax.set_xlim(-0.75, self.fields.efolds.max())
        # plt.tight_layout()

        if save:
            plt.savefig("%sturningrate_%s.pdf" % (plotpath, potential_name))
        if show:
            plt.show()
        plt.clf()
        plt.close()


# =========================================================================================


class PowerSpectrum(BackgroundEvolution):
    def __init__(self, Bg, NB=6, N_exit=15.0):
        # super().__init__(fields, p)
        # super().__call__(N_start, N_end, n_steps)
        self.tols = np.array([10 ** -15, 10 ** -15])
        self.Bg = Bg

        self.k = Bg.set_exit_mode(N_exit)
        self.N_exit = N_exit
        self.n_steps = Bg.steps

        self.initial_conditions(NB)

    def __call__(self):
        tols = np.array([1.0e-17, 1.0e-17])
        t_sigma = np.linspace(self.N_start, self.Bg.background[-1, 0], self.n_steps)
        self.two_pt = PyT.sigEvolve(
            t_sigma, self.k, self.ext_min, self.Bg.params, self.tols, True
        )

        # the second column is the 2pt of zeta
        zz1 = self.two_pt[:, 1]
        # the last 2nF* 2nF columns correspond to the evolution of the sigma matrix
        sigma = self.two_pt[:, 1 + 1 + 2 * self.Bg.nF :]

        self.t_sig = t_sigma
        self.sigma = sigma
        self.zz1 = zz1

        # self.spectral_index()

    def initial_conditions(self, NB, k=None):
        if k is None:
            k = self.k

        # print("initial input", NB, k)
        self.N_start, self.ext_min = PyS.ICsBE(
            NB, k, self.Bg.background, self.Bg.params, PyT
        )

    def disect_perturbations(self):
        N = self.two_pt.shape[0]
        P_zeta_k = self.two_pt[:, 1]
        Sig_ab_k = self.two_pt[:, -4 * self.Bg.nF * self.Bg.nF :].reshape(
            (N, 2 * self.Bg.nF, 2 * self.Bg.nF)
        )

        return P_zeta_k, Sig_ab_k

    def spectral_index(self):
        tols = np.array([1.0e-17, 1.0e-17])
        two_pt = PyT.sigEvolve(
            self.t_sig, self.k * 1.1, self.ext_min, self.Bg.params, self.tols, True
        )
        power_spec2 = two_pt[:, 1][-1]
        self.n_s = (np.log(power_spec2) - np.log(self.zz1[-1])) / (
            np.log(1.1 * self.k) - np.log(self.k)
        ) + 4.0

        twoG = PyT.gamEvolve(
            self.t_sig, self.k, self.ext_min, self.Bg.params, tols, True
        )
        self.r = (
            8 * twoG[-1, 1 + 1 + 2 * self.Bg.nF + 0 + 2 * self.Bg.nF * 0]
        ) / self.zz1[-1]

    def p_spectra(self, subhor=None, pool=False):
        kOut = np.array([])
        if subhor is None:
            subhor = self.N_exit
        for ii in range(1, 30):
            NExit = subhor + ii
            k = PyS.kexitN(NExit, self.Bg.background, self.Bg.params, PyT)
            kOut = np.append(kOut, k)

        if pool:
            Pz = pSpectra_pool(
                kOut, self.Bg.background, self.Bg.params, subhor, self.tols
            )
        else:
            Pz, __ = PyS.pSpectra(
                kOut, self.Bg.background, self.Bg.params, subhor, self.tols, PyT
            )  # this cacalcute P_z for this range of ks
        Pz = kOut ** 3 * Pz / (2 * np.pi)
        return Pz, kOut

    def plot_p_spectra(self, Pz, kOut, save=False, show=True):
        fig = plt.figure()
        plt.plot(np.log(kOut / kOut[0]), Pz, linewidth=2)
        plt.semilogy()
        plt.title("Power spectrum")
        plt.ylabel(r"$\cal P$ [Mpc$^3$]")
        plt.xlabel(r"$\log(k/k_{\rm pivot})$")
        plt.xlim(min(np.log(kOut / kOut[0])), max(np.log(kOut / kOut[0])))

        if save:
            plt.savefig(plotpath + "Pz.pdf")
        if show:
            plt.show()
        plt.clf()

    def plot_sigma_evolution(self, save=False, show=True):
        for ii in range(0, 2):
            for jj in range(0, 2):
                func = np.abs(self.sigma[:, ii + 2 * self.Bg.nF * jj])
                plt.plot(self.two_pt[:, 0], func)

        plt.title(r"$\Sigma$ evolution", fontsize=15)
        plt.ylabel(r"Aboslute 2pt field correlations")
        plt.xlabel(r"$N$")
        plt.grid()
        plt.yscale("log")
        # plt.legend(fontsize=15)
        if save:
            plt.savefig(plotpath + "2pt_sigma_evolve.pdf")
        if show:
            plt.show()
        plt.clf()

    def plot_zeta_evolution(self, save=False, show=True):
        plt.plot(self.t_sig, self.zz1)

        plt.title(r"$P_\zeta$ evolution")
        plt.ylabel(r"$P_\zeta(k)$")
        plt.xlabel(r"$N$")
        plt.grid(True)
        plt.yscale("log")
        # plt.legend(fontsize=15)
        if save:
            plt.savefig(plotpath + "2pt_zeta_evolve.pdf")
        if show:
            plt.show()
        plt.clf()

    def plot_background(self, save=False, show=True):
        plot_background(super(), save, show)


class Bispectrum(BackgroundEvolution):
    def __init__(
        self,
        background,
        alpha,
        beta,
        NB=6.0,
        N_exit=15.0,
        n_steps=1000,
    ):
        self.tols = np.array([10 ** -15, 10 ** -15])
        self.alpha = alpha
        self.beta = beta
        self.Bg = background
        self.n_steps = background.steps
        # super().__init__(fields, p)
        # super().__call__(N_start, N_end, n_steps)
        self.k = self.Bg.set_exit_mode(N_exit)

        self.k1 = self.k / 2 * (1 - beta)
        self.k2 = self.k / 4 * (1 + alpha + beta)
        self.k3 = self.k / 4 * (1 - alpha + beta)
        self.initial_conditions(NB)
        # print("start, exit_min", self.N_start, self.ext_min)
        self.n_steps = n_steps

    def initial_conditions(self, NB):
        self.NB = NB
        k_M = np.min(np.array([self.k1, self.k2, self.k3]))
        self.N_start, self.ext_min = PyS.ICsBM(
            NB, k_M, self.Bg.background, self.Bg.p_value, PyT
        )

    def __call__(self):
        t_alpha = np.linspace(self.N_start, self.Bg.background[-1, 0], self.n_steps)
        self.three_pt = PyT.alphaEvolve(
            t_alpha,
            self.k1,
            self.k2,
            self.k3,
            self.ext_min,
            self.Bg.p_value,
            self.tols,
            True,
        )

        # this contains the evolution of two point of zeta for each k mode involved and the 3pt of zeta
        zzz = self.three_pt[:, 1:5]
        self.three_pt_pert = self.three_pt[
            :, 1 + 4 + 2 * self.Bg.nF + 6 * 2 * self.Bg.nF * 2 * self.Bg.nF :
        ]
        self.fnl = (
            5.0
            / 6.0
            * zzz[:, 3]
            / (zzz[:, 1] * zzz[:, 2] + zzz[:, 0] * zzz[:, 1] + zzz[:, 0] * zzz[:, 2])
        )
        self.zzz = zzz
        self.t_alpha = t_alpha

    def disect_perturbations(self):
        N = self.three_pt.shape[0]

        i = 3
        P_zeta_k = self.three_pt[:, [1, 2, i]]
        i += 1
        B_zeta_kkk = self.three_pt[:, i]
        i += 4  # 2nF fields and field velocities
        Sig_ab_k1_Re = (self.three_pt[:, i : i + 16]).reshape(N, 4, 4)  # 16 = 2nF*2nF
        i += 16
        Sig_ab_k2_Re = (self.three_pt[:, i : i + 16]).reshape(N, 4, 4)
        i += 16
        Sig_ab_k3_Re = (self.three_pt[:, i : i + 16]).reshape(N, 4, 4)
        i += 16

        Sig_ab_k1_Im = (self.three_pt[:, i : i + 16]).reshape(N, 4, 4)  # 16 = 2nF*2nF
        i += 16
        Sig_ab_k2_Im = (self.three_pt[:, i : i + 16]).reshape(N, 4, 4)
        i += 16
        Sig_ab_k3_Im = (self.three_pt[:, i : i + 16]).reshape(N, 4, 4)
        i += 16

        B_abc_kkk = self.three_pt[:, i : i + 64]  # .reshape(N,4,4,4) # 64 = 2nF*2nF*2nF

        Sig_ab_k1 = Sig_ab_k1_Re + 1.0j * Sig_ab_k1_Im
        Sig_ab_k2 = Sig_ab_k2_Re + 1.0j * Sig_ab_k2_Im
        Sig_ab_k3 = Sig_ab_k3_Re + 1.0j * Sig_ab_k3_Im

        self.fnl = (
            5.0
            / 6.0
            * B_zeta_kkk
            / (
                P_zeta_k[:, 1] * P_zeta_k[:, 2]
                + P_zeta_k[:, 0] * P_zeta_k[:, 1]
                + P_zeta_k[:, 0] * P_zeta_k[:, 2]
            )
        )
        return P_zeta_k, B_zeta_kkk, B_abc_kkk, [Sig_ab_k1, Sig_ab_k2, Sig_ab_k3]

    def plot_fnl(self, save=False, show=True):
        plt.plot(self.t_alpha, self.fnl, "r")

        plt.title(r"$f_{NL}$ evolution")
        plt.ylabel(r"$f_{NL}$")
        plt.xlabel(r"$N$")
        # plt.legend(fontsize=15)
        plt.grid()
        if save:
            plt.savefig(plotpath + "f_nl.pdf")
        if show:
            plt.show()
        plt.clf()

    def plot_three_pt_correlations(self, save=False, show=True):
        for ii in range(0, 2):
            for jj in range(0, 2):
                for kk in range(0, 2):
                    func = np.abs(
                        self.three_pt_pert[
                            :, int(ii + 2 * self.Bg.nF * (jj + 2 * self.Bg.nF * kk))
                        ]
                    )
                    plt.plot(self.t_alpha, func)

        plt.title(r"$\alpha$ evolution")
        plt.ylabel(r"Absolute 3pt field correlations")
        plt.xlabel(r"$N$")
        plt.yscale("log")
        # plt.legend(fontsize=15)
        plt.grid()
        if save:
            plt.savefig(plotpath + "3pt_correlation_abs.pdf")
        if show:
            plt.show()
        plt.clf()

    def plot_bispectrum(self, save=False, show=True):
        P_zeta_k, B_zeta_kkk, B_abc_kkk, Sig_ab = self.disect_perturbations()
        for i in range(B_abc_kkk.shape[1]):
            plt.plot(np.abs(B_abc_kkk[:, i]), label=str(i), lw=1, alpha=0.5)

        plt.semilogy()
        plt.legend()
        plt.savefig(plotpath + "B_abc_test.pdf")
        plt.clf()

        plt.plot(P_zeta_k[:, 0])
        plt.plot(P_zeta_k[:, 1])
        plt.plot(P_zeta_k[:, 2])
        plt.semilogy()
        plt.savefig(plotpath + "powerspec_test.pdf")
        plt.clf()

    def plot_background(self, save=False, show=True):
        super().plot_background(save, show)


def plots_for_background(
    Bg, save=False, show=True, to_xy=True, normalisation=1, box_aspect=[1, 1, 1]
):
    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    ax.set_box_aspect(box_aspect)
    # ax.set_axis_off()
    # normalisation = 1.0 / Bg.p_value[1]

    # Make data.
    X = (
        np.linspace(
            0.15 * Bg.fields.value[:, 0].min(), 1.8 * Bg.fields.value[:, 0].max(), 3000
        )
        * normalisation
    )
    Y = np.linspace(
        0.15 * Bg.fields.value[:, 1].min(), 1.8 * Bg.fields.value[:, 1].max(), 3000
    )
    space = np.asarray(np.meshgrid(X, Y))
    X, Y = space[0], space[1]

    if to_xy:
        X, Y = space[0] * np.cos(space[1]), space[0] * np.sin(space[1])

    pot = Bg.V_sympy(space.reshape((space.shape[0], -1)), Bg.params).reshape(
        (space.shape[1:])
    )  # /Bg.params[1]**2

    # pot[pot>1e2] = np.nan
    # pot[pot<0] = np.nan
    # pot[X>.65] = np.nan
    # pot[X<.35] = np.nan

    # print(pot)
    surf = ax.plot_surface(
        X,
        Y,
        pot,
        cmap=cm.gist_heat,
        linewidth=1,
        antialiased=True,
        alpha=0.6,
        zorder=1,
        # vmax=70, vmin=0,
    )
    fig.colorbar(surf, shrink=0.5, aspect=5, label="$V(\\phi)$")

    X = Bg.fields.value[:, 0] * normalisation
    Y = Bg.fields.value[:, 1]
    field = Bg.fields.value
    field[:, 0] *= normalisation
    print(X.shape, Y.shape, field.shape)
    Z = Bg.V_sympy(field.T, Bg.params)  # /Bg.params[1]**2
    if to_xy:
        X, Y = X * np.cos(Y), X * np.sin(Y)
    ax.scatter(X, Y, Z, s=10, color="black", zorder=16)

    # Customize the z axis.
    # ax.set_zlim(0, 100)
    # ax.set_xlim(.35, .65)

    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter("{x:.02f}")
    ax.set_xlabel("$\\phi^1$")
    ax.set_ylabel("$\\phi^2$")

    if save:
        plt.savefig(plotpath + "3dplot.pdf", transparent=True)
    if show:
        plt.show()
    plt.clf()
    plt.close()


############################################################################################################################################
def pSpectra_loop(params, ii):
    kA, back, params, NB, tols = params
    zzOut = np.array([])
    # times = np.array([])
    num = np.size(kA)
    print("\n \n \n performing " + str(ii + 1) + " of " + str(num) + "\n \n \n")
    k = kA[ii]
    Nstart, backExitMinus = PyS.ICs(NB, k, back, params, PyT)
    # start_time = timeit.default_timer()

    if Nstart == np.nan:
        twoPt = numpy.empty((2, 2))
        twoPt[:] = np.nan
    else:
        t = np.linspace(Nstart, back[-1, 0], 10)
        # run solver for this triangle
        twoPt = PyT.sigEvolve(
            t, k, backExitMinus, params, tols, True
        )  # all data from three point run goes into threePt array
    # zzOut=np.append(zzOut, twoPt[-1,1])
    # print(zzOut)
    # times = np.append(times, timeit.default_timer()-start_time)
    return twoPt[-1, 1]


def pSpectra_pool(kA, back, params, NB, tols):
    num = np.size(kA)
    from multiprocessing import Pool

    with Pool(processes=4) as pool:
        return pool.map(
            partial(pSpectra_loop, (kA, back, params, NB, tols)), range(0, num)
        )
