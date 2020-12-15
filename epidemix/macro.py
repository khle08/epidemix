########################################################################

'''

Reference:
+ scipy求解微分方程: https://www.twblogs.net/a/5b7afb882b7177539c248a16

Author: Kuo Chun-Lin
Date  : 2020.05.01
'''

########################################################################

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.integrate import odeint

########################################################################


class MacroODE(object):
    N = None
    initial = None

    def derivative(self, z, t):
        assert self.N is not None, 'N is not defined.'
        assert self.initial is not None, 'States are not initiated.'
        raise NotImplementedError()


class SI(MacroODE):
    def __init__(self, N, I0, beta):
        self.N = N
        self.I0 = I0
        self.S0 = N - I0
        self.initial = (self.S0, I0)

        # Beta: transmission / infection Rate
        self.beta = beta

    def derivative(self, y, t):
        S, I = y
        dSdt = -self.beta * S * I / self.N
        dIdt = self.beta * S * I / self.N
        return dSdt, dIdt


class SIS(MacroODE):
    def __init__(self, N, I0, beta, gamma):
        self.N = N
        self.I0 = I0
        self.S0 = N - I0
        self.initial = (self.S0, I0)

        self.beta = beta
        self.gamma = gamma
        self.reproduction_num = beta / gamma    # Definition of "R_0".

    def derivative(self, y, t):
        S, I = y
        dSdt = -self.beta * S * I / self.N + self.gamma * I
        dIdt = self.beta * S * I / self.N - self.gamma * I
        return dSdt, dIdt


class SIR(MacroODE):
    """docstring for SIR"""

    def __init__(self, N, I0, R0, beta, gamma):
        self.N = N
        self.I0 = I0
        self.R0 = R0
        self.S0 = N - I0 - R0
        self.initial = (self.S0, I0, R0)

        self.beta = beta
        self.gamma = gamma
        self.reproduction_num = beta / gamma    # Definition of "R_0".

    def derivative(self, y, t):
        S, I, R = y
        dSdt = -self.beta * S * I / self.N
        dIdt = self.beta * S * I / self.N - self.gamma * I
        dRdt = self.gamma * I
        return dSdt, dIdt, dRdt


class SIRS(MacroODE):
    def __init__(self, N, I0, R0, beta, gamma, delta):
        self.N = N
        self.I0 = I0
        self.R0 = R0
        self.S0 = N - I0 - R0
        self.initial = (self.S0, I0, R0)

        self.beta = beta
        self.gamma = gamma
        self.delta = delta
        self.reproduction_num = None    # Definition of "R_0".

    def derivative(self, y, t):
        S, I, R = y
        dSdt = -self.beta * S * I / self.N + self.delta * R
        dIdt = self.beta * S * I / self.N - self.gamma * I
        dRdt = self.gamma * I - self.delta * R
        return dSdt, dIdt, dRdt


class SEIR(MacroODE):
    def __init__(self, N, E0, I0, R0, sigma, beta, gamma):
        self.N = N
        self.E0 = E0
        self.I0 = I0
        self.R0 = R0
        self.S0 = N - E0 - I0 - R0
        self.initial = (self.S0, E0, I0, R0)

        self.sigma = sigma
        self.beta = beta
        self.gamma = gamma

    def derivative(self, y, t):
        S, E, I, R = y
        dSdt = -self.beta * S * I / self.N
        dEdt = self.beta * S * I / self.N - self.sigma * E
        dIdt = self.sigma * E - self.gamma * I
        dRdt = self.gamma * I
        return dSdt, dEdt, dIdt, dRdt


class SEIRD(MacroODE):
    def __init__(self, N, E0, I0, R0, D0, sigma, beta, gamma, p, rho):
        self.N = N
        self.E0 = E0
        self.I0 = I0
        self.R0 = R0
        self.D0 = D0
        self.S0 = N - E0 - I0 - R0 - D0
        self.initial = (self.S0, E0, I0, R0, D0)

        self.sigma = sigma
        self.beta = beta
        self.gamma = gamma
        self.p = p
        self.rho = rho

    def derivative(self, y, t):
        S, E, I, R, D = y
        dSdt = - self.beta * S * I / self.N
        dEdt = self.beta * S * I / self.N - self.sigma * E
        dIdt = self.sigma * E - (1 - self.p) * self.gamma * I \
                              - self.p * self.rho * I
        dRdt = (1 - self.p) * self.gamma * I
        dDdt = self.p * self.rho * I
        return dSdt, dEdt, dIdt, dRdt, dDdt


class EpiModel(object):
    """docstring for EpiModel"""

    def __init__(self, model, data_dir=None):
        self.model = model
        self.data_dir = data_dir

    def simulate(self, t):
        results = odeint(self.model.derivative, self.model.initial, t)
        return results.T


########################################################################


if __name__ == '__main__':
    # Seq : N, I0, beta
    si = SI(1000, 10, 0.2)
    # Seq : N, I0, beta, gamma
    sis = SIS(1000, 10, 0.2, 0.1)
    # Seq : N, I0, R0, beta, gamma
    sir = SIR(1000, 50, 0, 0.4, 0.1)
    # Seq : N, E0, I0, R0, D0, beta, gamma, delta, alpha, rho
    # seird = SEIRD(1000, 1, 0, 0, 0, )

    days = np.linspace(0, 80, 80)

    # ------------------------------------------------------------------

    # epi = EpiModel(si)
    # s, i = epi.simulate(days)

    # # Plot the data on three separate curves for S(t), I(t) and R(t)
    # fig = plt.figure(facecolor='w')
    # plt.plot(days, s / si.N, 'b', alpha=0.5, lw=2, label='Susceptible')
    # plt.plot(days, i / si.N, 'r', alpha=0.5, lw=2, label='Infected')
    # plt.xlabel('Time /days')
    # plt.ylabel('Number (1000s)')
    # plt.grid(True)
    # plt.legend()

    # plt.show()

    # ------------------------------------------------------------------

    epi = EpiModel(sir)
    s, i, r = epi.simulate(days)

    # Plot the data on three separate curves for S(t), I(t) and R(t)
    fig = plt.figure(facecolor='w')
    plt.plot(days, s / sis.N, 'b', alpha=0.5, lw=2, label='Susceptible')
    plt.plot(days, i / sis.N, 'r', alpha=0.5, lw=2, label='Infected')
    plt.plot(days, r / sis.N, 'g', alpha=0.5, lw=2, label='Recovered')
    plt.xlabel('Time /days')
    plt.ylabel('Number (1000s)')
    plt.grid(True)
    plt.legend()

    plt.show()

    # ------------------------------------------------------------------

    # epi = EpiModel(sir)
    # s, i, r = epi.simulate(days)

    # # Plot the data on three separate curves for S(t), I(t) and R(t)
    # fig = plt.figure(facecolor='w')
    # plt.plot(days, s / sir.N, 'b', alpha=0.5, lw=2, label='Susceptible')
    # plt.plot(days, i / sir.N, 'r', alpha=0.5, lw=2, label='Infected')
    # plt.plot(days, r / sir.N, 'g', alpha=0.5, lw=2,
    #          label='Recovered with immunity')
    # plt.xlabel('Time /days')
    # plt.ylabel('Number (1000s)')
    # plt.grid(True)
    # plt.legend()

    # plt.show()
