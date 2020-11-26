########################################################################

'''

Reference:
+ scipy求解微分方程: https://www.twblogs.net/a/5b7afb882b7177539c248a16

Author: Kuo Chun-Lin
Date  : 2020.06.06 - 2020.08.11
'''

########################################################################

import copy
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

from scipy.integrate import odeint

########################################################################


class DifferentialEquation(object):
    # The 2 necessary parameters that are easily being ignored.
    A = None
    N = None
    initial = None

    def derivative(self, z, t):
        raise NotImplementedError()

    @property
    def model_name(self):
        return self.__class__.__name__

    @property
    def degree_base(self):
        assert self.N is not None, 'N is not defined.'
        assert self.A is not None, 'A is not defined.'
        vec = np.sum(self.A, axis=0)
        return self.N * vec / np.sum(vec)

    @property
    def acq_base(self):
        assert self.N is not None, 'N is not defined.'
        assert self.A is not None, 'A is not defined.'
        vec = np.sum(self.A.sum(axis=0) * self.A, axis=1)
        return self.N * vec / np.sum(vec)

    def first_nei(self, state):
        assert self.N is not None, 'N is not defined.'
        assert self.A is not None, 'A is not defined.'
        state = self.N * state / np.sum(state)
        vec = np.sum(self.A, axis=0) * state
        return self.N * vec / np.sum(vec)

    def second_nei(self, state):
        assert self.N is not None, 'N is not defined.'
        assert self.A is not None, 'A is not defined.'
        state = self.N * state / np.sum(state)
        vec = np.sum(self.A.sum(axis=0) * self.A, axis=1) * state
        return self.N * vec / np.sum(vec)

    def inv_idx(self, i, total):
        idx = np.ones(total, dtype=np.bool)
        idx[i] = 0
        return idx


class SI(DifferentialEquation):
    def __init__(self, A, I0, beta):
        # numpy 2D Adjacent matrix
        self.A = A
        self.N = len(A)

        # Randomly assign the non-repeated infected and recovered nodes.
        idx = np.random.choice(np.arange(self.N), I0, replace=False)
        self.I0 = np.zeros((self.N,))
        self.I0[idx] = 1

        # Init matrix should be stacked into a 1D array.
        self.initial = np.hstack([1 - self.I0, self.I0])
        self.beta = beta

    def derivative(self, z, t):
        # The initial "z" starts from "self.initial".
        b = self.beta * z[0:self.N] * np.dot(self.A, z[self.N:2 * self.N])
        return np.hstack([-b, b])


class SIS(DifferentialEquation):
    def __init__(self, A, I0, beta, gamma):
        # numpy 2D Adjacent matrix
        self.A = A
        self.N = len(A)

        # Randomly assign the non-repeated infected and recovered nodes.
        idx = np.random.choice(np.arange(self.N), I0, replace=False)
        self.I0 = np.zeros((self.N,))
        self.I0[idx] = 1

        # Init state array is the combination of "s(t)" and "x(t)".
        self.initial = np.hstack([1 - self.I0, self.I0])

        self.beta = beta
        self.gamma = gamma
        self.reproduction_num = beta / gamma    # Definition of "R_0".

    def derivative(self, z, t):
        # The initial "z" starts from "self.initial".
        b = self.beta * z[0:self.N] * np.dot(self.A, z[self.N:2 * self.N])
        r = self.gamma * z[self.N:2 * self.N]
        return np.concatenate([-b + r, b - r])


class SIR(DifferentialEquation):
    def __init__(self, A, I0, R0, beta, gamma):
        # numpy 2D Adjacent matrix
        self.A = A
        self.N = len(A)

        # Randomly assign the non-repeated infected and recovered nodes.
        idx = np.random.choice(np.arange(self.N), I0 + R0, replace=False)
        self.I0 = np.zeros((self.N,))
        self.R0 = np.zeros((self.N,))
        self.I0[idx[:I0]] = 1
        self.R0[idx[I0:I0 + R0]] = 1

        # Init matrix should be stacked into a 1D array.
        self.initial = np.hstack([1 - self.I0 - self.R0,    # s(t)
                                  self.I0,                  # i(t)
                                  self.R0])                 # r(t)
        self.beta = beta
        self.gamma = gamma
        self.reproduction_num = beta / gamma    # Definition of "R_0".

    def derivative(self, z, t):
        # The initial "z" starts from "self.initial".
        b = self.beta * z[0:self.N] * np.dot(self.A, z[self.N:2 * self.N])
        r = self.gamma * z[self.N:2 * self.N]
        return np.hstack([-b, b - r, r])


class SIRV(DifferentialEquation):
    """docstring for SIRV"""

    def __init__(self, A, I0, R0, V0, beta, gamma, eps):
        # numpy 2D Adjacent matrix
        self.A = A
        self.N = len(A)

        # Randomly assign the non-repeated infected, recovered, and vacci nodes.
        idx = np.random.choice(np.arange(self.N), I0 + R0 + V0, replace=False)
        self.I0 = np.zeros((self.N,))
        self.R0 = np.zeros((self.N,))
        self.V0 = np.zeros((self.N,))
        self.I0[idx[:I0]] = 1
        self.R0[idx[I0:I0 + R0]] = 1
        self.V0[idx[I0 + R0:I0 + R0 + V0]] = 1

        # Init matrix should be stacked into a 1D array.
        self.initial = np.hstack([1 - self.I0 - self.R0 - self.V0,
                                  self.I0, self.R0, self.V0])
        self.beta = beta
        self.gamma = gamma
        self.eps = eps
        self.reproduction_num = beta / gamma

    def derivative(self, z, t):
        Pk = 1    # default.
        # Pk = self.degree_base
        # Pk = self.acq_base
        # Pk = self.first_nei(z[self.N * 1:2 * self.N])
        # Pk = self.second_nei(z[self.N * 1:2 * self.N])

        b = self.beta * z[0:self.N] * np.dot(self.A, z[self.N:2 * self.N])
        r = self.gamma * z[self.N:2 * self.N]
        v = (1 - self.beta) * self.eps * z[0:self.N] * Pk
        return np.hstack([-b - v, b - r, r, v])


class SIRVir1(DifferentialEquation):
    def __init__(self, A, I0, R0, V0, i0, r0, beta1, gamma1, alpha1,
                 beta2, gamma2, alpha2, eps):
        # numpy 2D Adjacent matrix
        self.A = A
        self.N = len(A)

        idx = np.random.choice(np.arange(self.N),
                               I0 + R0 + V0 + i0 + r0,
                               replace=False)
        self.I0 = np.zeros((self.N,))
        self.R0 = np.zeros((self.N,))
        self.V0 = np.zeros((self.N,))
        self.i0 = np.zeros((self.N,))
        self.r0 = np.zeros((self.N,))
        self.I0[idx[:I0]] = 1
        self.R0[idx[I0:I0 + R0]] = 1
        self.V0[idx[I0 + R0:I0 + R0 + V0]] = 1
        self.i0[idx[I0 + R0 + V0:I0 + R0 + V0 + i0]] = 1
        self.r0[idx[I0 + R0 + V0 + i0:I0 + R0 + V0 + i0 + r0]] = 1

        self.initial = np.hstack(
            [1 - self.I0 - self.R0 - self.V0 - self.i0 - self.r0,
             self.I0, self.R0, self.V0, self.i0, self.r0]
        )

        self.beta1 = beta1
        self.beta2 = beta2
        self.gamma1 = gamma1
        self.gamma2 = gamma2
        self.alpha1 = alpha1
        self.alpha2 = alpha2
        self.eps = eps

    def derivative(self, z, t):
        Pk = 1    # default.

        # # Degree-based Target Immunization
        # Pk = np.zeros(self.N)
        # thres = 0.1
        # Pk[self.get_nodes(thres)] = 1

        # # Acquaintance Immunization
        # Pk = np.zeros(self.N)
        # thres = 0.1
        # Pk[self.get_neighbors(thres)] = 1

        b1 = self.beta1 * z[:1 * self.N] * \
            np.dot(self.A, z[self.N * 1:2 * self.N])
        r1 = self.gamma1 * z[self.N * 1:2 * self.N]
        a1 = self.alpha1 * z[self.N * 2:3 * self.N]
        v = self.eps * z[:1 * self.N] * Pk

        b2 = self.beta2 * z[self.N * 3:4 * self.N] * \
            np.dot(self.A, z[self.N * 1:2 * self.N])
        r2 = self.gamma2 * z[self.N * 4:5 * self.N]
        a2 = self.alpha2 * z[self.N * 5:]
        return np.hstack([-b1 - v + a1 + a2,
                          b1 - r1,
                          r1 - a1,
                          v - b2,
                          b2 - r2,
                          r2 - a2])


class SIRVir2(DifferentialEquation):
    def __init__(self, A, I0, R0, V0, i0, r0, beta1, gamma1, alpha1,
                 beta2, gamma2, alpha2, eps, eps1v, eps2v):
        # numpy 2D Adjacent matrix
        self.A = A
        self.N = len(A)

        idx = np.random.choice(np.arange(self.N),
                               I0 + R0 + V0 + i0 + r0,
                               replace=False)
        self.I0 = np.zeros((self.N,))
        self.R0 = np.zeros((self.N,))
        self.V0 = np.zeros((self.N,))
        self.i0 = np.zeros((self.N,))
        self.r0 = np.zeros((self.N,))
        self.I0[idx[:I0]] = 1
        self.R0[idx[I0:I0 + R0]] = 1
        self.V0[idx[I0 + R0:I0 + R0 + V0]] = 1
        self.i0[idx[I0 + R0 + V0:I0 + R0 + V0 + i0]] = 1
        self.r0[idx[I0 + R0 + V0 + i0:I0 + R0 + V0 + i0 + r0]] = 1

        self.initial = np.hstack(
            [1 - self.I0 - self.R0 - self.V0 - self.i0 - self.r0,
             self.I0, self.R0, self.V0, self.i0, self.r0]
        )

        self.beta1 = beta1
        self.beta2 = beta2
        self.gamma1 = gamma1
        self.gamma2 = gamma2
        self.alpha1 = alpha1
        self.alpha2 = alpha2
        self.eps = eps
        self.eps1v = eps1v
        self.eps2v = eps2v

    def derivative(self, z, t):
        Pk = 1    # default.
        P1k = 1
        P2k = 1

        b1 = self.beta1 * z[:1 * self.N] * \
            np.dot(self.A, z[self.N * 1:2 * self.N])
        r1 = self.gamma1 * z[self.N * 1:2 * self.N]
        a1 = self.alpha1 * z[self.N * 2:3 * self.N]
        e1 = self.eps1v * z[self.N * 2:3 * self.N] * P1k
        v = self.eps * z[:1 * self.N] * Pk

        b2 = self.beta2 * z[self.N * 3:4 * self.N] * \
            np.dot(self.A, z[self.N * 1:2 * self.N])
        r2 = self.gamma2 * z[self.N * 4:5 * self.N]
        a2 = self.alpha2 * z[self.N * 5:]
        e2 = self.eps2v * z[self.N * 5:] * P2k
        return np.hstack([-b1 - v + a1 + a2,
                          b1 - r1,
                          r1 - a1 - e1,
                          v - b2 + e1 + e2,
                          b2 - r2,
                          r2 - a2 - e2])


class SEINRVeinr(DifferentialEquation):
    def __init__(self, A, E0, I0, N0, R0, V0, e0, i0, n0, r0,
                 beta1, theta1, eta1, gamma1, pi1, alpha1, eps,
                 beta2, theta2, eta2, gamma2, pi2, alpha2):
        # numpy 2D Adjacent matrix
        self.A = A
        self.N = len(A)

        idx = np.random.choice(np.arange(self.N),
                               I0 + R0 + V0 + i0 + r0,
                               replace=False)
        self.E0 = np.zeros((self.N,))
        self.I0 = np.zeros((self.N,))
        self.N0 = np.zeros((self.N,))
        self.R0 = np.zeros((self.N,))
        self.V0 = np.zeros((self.N,))
        self.e0 = np.zeros((self.N,))
        self.i0 = np.zeros((self.N,))
        self.n0 = np.zeros((self.N,))
        self.r0 = np.zeros((self.N,))

        num = 0
        self.E0[idx[num:num + E0]] = 1
        num += E0
        self.I0[idx[num:num + I0]] = 1
        num += I0
        self.N0[idx[num:num + N0]] = 1
        num += N0
        self.R0[idx[num:num + R0]] = 1
        num += R0
        self.V0[idx[num:num + V0]] = 1
        num += V0
        self.e0[idx[num:num + e0]] = 1
        num += e0
        self.i0[idx[num:num + i0]] = 1
        num += i0
        self.n0[idx[num:num + n0]] = 1
        num += n0
        self.r0[idx[num:num + r0]] = 1

        self.initial = np.hstack(
            [1 - self.E0 - self.I0 - self.N0 - self.R0 - self.V0 -
             self.e0 - self.i0 - self.n0 - self.r0,
             self.E0, self.I0, self.N0, self.R0, self.V0,
             self.e0, self.i0, self.n0, self.r0]
        )

        self.beta1, self.beta2 = beta1, beta2
        self.theta1, self.theta2 = theta1, theta2
        self.eta1, self.eta2 = eta1, eta2
        self.gamma1, self.gamma2 = gamma1, gamma2
        self.pi1, self.pi2 = pi1, pi2
        self.alpha1, self.alpha2 = alpha1, alpha2
        self.eps = eps

    def derivative(self, z, t):
        # 0, 1, 2, 3, 4, 5, 6, 7, 8, 9
        # S, E, I, N, R, V, e, i, n, r
        Pk = 1    # default.

        # Degree-based Target Immunization
        Pk = np.zeros(self.N)
        thres = 10
        Pk[self.get_nodes(thres)] = 1

        # # Acquaintance Immunization
        # Pk = np.zeros(self.N)
        # thres = 0.1
        # Pk[self.get_neighbors(thres)] = 1

        b1 = self.beta1 * z[:1 * self.N] * \
            np.dot(self.A, z[self.N * 2:3 * self.N])
        th1 = self.theta1 * z[self.N * 1:2 * self.N]
        e1 = (1 - self.theta1) * self.eta1 * z[self.N * 1:2 * self.N]
        r1 = self.gamma1 * z[self.N * 2:3 * self.N]
        p1 = self.pi1 * z[self.N * 3:4 * self.N]
        a1 = self.alpha1 * z[self.N * 4:5 * self.N]
        v = (1 - self.beta1) * self.eps * z[:1 * self.N] * Pk

        b2 = self.beta2 * z[self.N * 5:6 * self.N] * \
            np.dot(self.A, z[self.N * 2:3 * self.N])
        th2 = self.theta2 * z[self.N * 6:7 * self.N]
        e2 = (1 - self.theta2) * self.eta2 * z[self.N * 6:7 * self.N]
        r2 = self.gamma2 * z[self.N * 7:8 * self.N]
        p2 = self.pi2 * z[self.N * 8:9 * self.N]
        a2 = self.alpha2 * z[self.N * 9:]
        return np.hstack([-b1 - v + a1 + a2,  # S
                          b1 - th1 - e1,      # E
                          th1 - r1,           # I
                          e1 - p1,            # N
                          r1 + p1 - a1,       # R
                          v - b2,             # V
                          b2 - th2 - e2,      # e
                          th2 - r2,           # i
                          e2 - p2,            # n
                          r2 + p2 - a2])      # r


class SIRS(DifferentialEquation):
    def __init__(self, A, I0, R0, beta, gamma, delta):
        # numpy 2D Adjacent matrix
        self.A = A
        self.N = len(A)

    def derivative(self, z, t):
        pass


class SEIR(DifferentialEquation):
    def __init__(self, A, E0, I0, R0, sigma, beta, gamma):
        # numpy 2D Adjacent matrix
        self.A = A
        self.N = len(A)

    def derivative(self, z, t):
        pass


class SEIRD(DifferentialEquation):
    def __init__(self, A, E0, I0, R0, D0, sigma, beta, gamma, p, rho):
        # numpy 2D Adjacent matrix
        self.A = A
        self.N = len(A)

    def derivative(self, z, t):
        pass


########################################################################


class EpiModel(object):
    """ Take the instance of a Epidemic model as the input for a
    simulation within a period. """

    def __init__(self, G, eq, num_state, params, data_dir=None,
                 state_colors=['blue', 'red', 'green']):
        self.G = G
        self.A = np.array(nx.adj_matrix(G).todense())
        self.N = len(self.A)
        self.model = eq(self.A, *params)

        self.num_state = num_state
        self.data_dir = data_dir

        assert self.model.initial is not None, 'Initial vector is not defined.'

    def simulate(self, t):
        # Put the model into a differential equation solver.
        results = odeint(self.model.derivative, self.model.initial, t)
        # Split the results into a proper number of part.
        return np.split(results, self.num_state, axis=1)


########################################################################


if __name__ == '__main__':
    # Seq: (The number of nodes,
    #       The number of random edges to add for each new node,
    #       Probability of adding a triangle after adding a random edge)
    G_ws = nx.watts_strogatz_graph(40, 5, 0.4)
    # G_pl = nx.powerlaw_cluster_graph(50, 5, 0.4)
    A = np.array(nx.adj_matrix(G_ws).todense())
    days = np.arange(0, 10, 0.1)

    # ------------------------------------------------------------------

    # # Seq : Graph, I0, beta
    # si = SI(A, 4, 0.6)
    # # Seq : Graph, I0, beta, gamma
    # sis = SIS(A, 4, 0.6, 0.1)
    # # Seq : Graph, I0, R0, beta, gamma
    # sir = SIR(A, 4, 2, 0.4, 0.6)
    # # Seq: Graph, I0, R0, V0, beta, gamma, eps
    # sirv = SIRV(A, 4, 2, 2, 0.6, 0.2, 0.05)

    # model = sir
    # num_state = 3

    # simulate = odeint(model.derivative, model.initial, days)
    # s, i, r = np.array(np.split(simulate, num_state, axis=1))

    # ------------------------------------------------------------------

    epi = EpiModel(G_ws, SIR, num_state=3, params=[4, 2, 0.4, 0.6])
    s, i, r = epi.simulate(days)
    # epi.step(1)

    # Plot the data on three separate curves for S(t) and I(t).
    fig = plt.figure(facecolor='w')
    plt.plot(s, color='blue')
    plt.plot(i, color='red')
    plt.plot(r, color='green')
    # plt.plot(v, color='orange')
    plt.xlabel('Time unit', fontsize=15)
    plt.ylabel('Probability (%)', fontsize=15)
    plt.xlim((0, 100))
    plt.ylim((0, 1))
    plt.grid(True)
    plt.legend()
    plt.show()
