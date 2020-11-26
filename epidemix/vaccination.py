########################################################################

'''
A place to define vaccination strategies.

Author: Kuo Chun-Lin
Date  : 2020.06.06 - 2020.08.11
'''

########################################################################

import copy
import numpy as np
import networkx as nx

from math import ceil
from itertools import combinations
from cdlib import algorithms as algo

from epidemic import EpiModel
from utils.partition import *

########################################################################


class VacciModel(EpiModel):
    """docstring for VacciModel"""

    def __init__(self, G, eq, num_state, params,
                 data_dir=None, vacci_strategy=None,
                 state_colors=['blue', 'red', 'green', 'orange']):
        super(VacciModel, self).__init__(
            G, eq, num_state, params, data_dir, state_colors)

        if vacci_strategy == 'random':
            self.vacci_func = self.random_vaccination
        elif vacci_strategy == 'target':
            self.vacci_func = self.degree_based_target
        elif vacci_strategy == 'acquaintance':
            self.vacci_func = self.acquaintance
        elif vacci_strategy == 'hhi':
            self.vacci_func = self.hhi

    # initialize a certain state to probability = 1.
    def initialize_nodes(self, idx_nodes, state=0):
        # Listing the total states included in the epidemic model and ...
        # ... loop through each state.
        for s in np.arange(self.num_state):
            if s == state:
                self.model.initial[s * self.N + idx_nodes] = 1
            else:
                self.model.initial[s * self.N + idx_nodes] = 0

    def community_vaccination(self, A, param=[8, 0.5]):
        # Use an algo to detect communities.
        communities = algo.louvain(nx.Graph(A)).communities
        # Get the nodes with highest "out_community" links.
        opt1 = target_nodes(communities, param[0], out_community=True)[0]
        opt2 = self.select_nodes(A, state=0, nei_thres=0)
        opt = np.intersect1d(opt1, opt2)

        idx = np.random.choice(opt, ceil(len(opt) * param[1]), replace=False)
        A[idx, :] = 0
        A[:, idx] = 0

        self.initialize_nodes(idx, state=0)
        self.vacci_num = len(idx)
        return A

    def random_vaccination(self, A, ratio):
        # opt = self.select_nodes(A, state=0, nei_thres=0)
        opt = self.select_neighbors(A, state=1, nei_state=0, nei_thres=0)

        # Randomly pick up certain selected nodes.
        idx = np.random.choice(opt, ceil(len(opt) * ratio), replace=False)
        A[idx, :] = 0
        A[:, idx] = 0

        # When the "s" nodes are separated, set the prob of "s" to 1 again.
        self.initialize_nodes(idx, state=0)
        self.vacci_num = len(idx)
        return A

    def degree_based_target(self, A, param):    # param = (kc, ratio)
        opt = self.select_nodes(A, state=0, nei_thres=param[0])
        # opt = self.select_neighbors(A, state=1, nei_state=0, nei_thres=param[0])

        idx = np.random.choice(opt, ceil(len(opt) * param[1]), replace=False)
        A[idx, :] = 0
        A[:, idx] = 0

        # When the "s" nodes are separated, set the prob of "s" to 1 again.
        self.initialize_nodes(idx, state=0)
        self.vacci_num = len(idx)
        return A

    def acquaintance(self, A, ratio):
        opt = self.select_nodes(A, state=0, nei_thres=0)
        # opt = self.select_neighbors(A, state=1, nei_state=0, nei_thres=0)

        # Randomly pick up certain selected nodes.
        idx = np.random.choice(opt, ceil(len(opt) * ratio), replace=False)
        # Pick up oonly those nodes with neighbors.
        neighbors = np.sum(A[idx], axis=0) > 0
        A[neighbors, :] = 0
        A[:, neighbors] = 0

        # When the "s" nodes are separated, set the prob of "s" to 1 again.
        self.initialize_nodes(idx, state=0)
        self.vacci_num = np.sum(neighbors)
        return A

    def hhi(self, A, ratio):
        opt = self.select_nodes(A, state=0, nei_thres=0)
        # opt = self.select_neighbors(A, state=1, nei_state=0, nei_thres=0)

        # Determine all possible combination of the selected nodes.
        subset = np.array(list(combinations(opt, ceil(len(opt) * ratio))))
        # Create an array to save the HHI value.
        H_list = np.zeros(len(subset))

        for i, idx_nodes in enumerate(subset):
            # Try each combination by changing its A matrix.
            A_temp = copy.deepcopy(A)
            A_temp[idx_nodes.astype(np.int), :] = 0
            A_temp[:, idx_nodes.astype(np.int)] = 0
            # Get the corresponding graph using the A matrix.
            G_temp = self.set_graph(self.state_list,
                                    self.color_list,
                                    G=nx.Graph(A_temp))
            # Find the total nodes of each subgraphs.
            subgraphs = list(nx.connected_components(G_temp))
            # Count how many nodes there are for each subgraph.
            node_num_list = np.array([len(z) for z in subgraphs])
            # Calculate HHI value using the following formula.
            H = 1 - (np.sum(node_num_list) ** 2) / (len(opt) ** 2)
            # Save the value to the array.
            H_list[i] = H

        # Find the one combination with minimum HHI value.
        idx = subset[np.argmin(H_list)].astype(np.int)
        # Change the A matrix according to the selected nodes.
        A[idx, :] = 0
        A[:, idx] = 0

        # When the "s" nodes are separated, set the prob of "s" to 1 again.
        self.initialize_nodes(idx, state=0)
        self.vacci_num = len(idx)
        return A


########################################################################


if __name__ == '__main__':

    from utils.plot import draw_probs
    from equations import SI, SIS, SIR, SIRV

    G = nx.watts_strogatz_graph(40, 5, 0.4)
    days = np.arange(0, 10, 0.1)
    colors = ['blue', 'red', 'green', 'orange']

    epi = VacciModel(G, SIR, num_state=3, params=[4, 2, 0.3, 0.1],
                     vacci_strategy='random', state_colors=colors)

    probs = epi.simulate(days)
    # probs = epi.reconstruct(days, 20, param=0.5)

    epi.set_propagate(0, 1, neighbor=1, update_state=False)
    epi.set_propagate(1, 2, neighbor=None, update_state=False)

    # status, _ = epi.propagation()
    status, probs = epi.propagation(reconstruct_index=[2, 4, 7, 20],
                                    reconstruct_param=0.5)

    epi.visualize(status, np.arange(16), figsize=(15, 15), n_row=4, n_col=4)
    draw_probs(probs, colors)
