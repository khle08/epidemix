########################################################################

'''
say something ....
'''

########################################################################

'''

Reference:
+ scipy求解微分方程: https://www.twblogs.net/a/5b7afb882b7177539c248a16

Author: Kuo Chun-Lin
Date  : 2020.06.06 - 2020.06.17
'''

########################################################################

import os
import copy
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

from tqdm import tqdm
from epidemic import EpiModel
from vaccination import VacciModel
from equations import SI, SIS, SIR, SIRV

from utils.plot import draw_probs

########################################################################

# SI   params: I0, beta
# SIS  params: I0, beta, gamma
# SIR  params: I0, R0, beta, gamma
# SIRV params: I0, R0, V0, beta, gamma, eps

# ----------------------------------------------------------------------


def main():
    # Seq: (The number of nodes,
    #       The number of random edges to add for each new node,
    #       Probability of adding a triangle after adding a random edge)
    G = nx.watts_strogatz_graph(50, 5, 0.2)
    # G = nx.powerlaw_cluster_graph(50, 5, 0.2)
    # G = nx.random_lobster(10, 0.6, 0.6, 0)  # 0: seed

    days = np.arange(0, 100, 1)
    colors = ['blue', 'red', 'green', 'orange']

    # epi = EpiModel(G, SIR, num_state=3, params=[4, 2, 0.3, 0.1])
    epi = VacciModel(G, SIR, num_state=3, params=[4, 2, 0.3, 0.1],
                     vacci_strategy='random', state_colors=colors)

    probs = epi.simulate(days)  # , param=0.5)
    # s, i, r = epi.reconstruct(days, 20, param=0.5)
    # s, i, r = epi.reconstruct(days, 40, param=0.5)
    # s, i, r = epi.reconstruct(days, 60, param=0.5)

    epi.set_propagate(0, 1, neighbor=1, update_state=False)
    epi.set_propagate(1, 2, neighbor=None, update_state=False)

    # status, _ = epi.propagation()
    status, probs = epi.propagation(reconstruct_index=[3, 7, 9],
                                    reconstruct_param=0.5)

    epi.visualize(status, np.arange(16), figsize=(15, 15), n_row=4, n_col=4)
    draw_probs(probs, colors)


########################################################################


if __name__ == '__main__':

    # ------------------------------------------------------------------
    main()    # For testing the propagation process.
    # ------------------------------------------------------------------

