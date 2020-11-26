########################################################################

'''

Reference:
+ scipy求解微分方程: https://www.twblogs.net/a/5b7afb882b7177539c248a16

Author: Kuo Chun-Lin
Date  : 2020.06.06 - 2020.06.17
'''

########################################################################

import copy
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

from scipy.integrate import odeint

########################################################################


class EpiModel(object):
    """ Take the instance of a Epidemic model as the input for a
    simulation within a period. """

    def __init__(self, G, eq, num_state, params, data_dir=None,
                 state_colors=['blue', 'red', 'green', 'orange']):
        self.G = G
        self.A = np.array(nx.adj_matrix(G).todense())
        self.N = len(self.A)
        self.model = eq(self.A, *params)
        # Get the node position so that the drawn graph can be fixed.
        self.pos = nx.spring_layout(self.G)
        self.num_state = num_state
        self.state_queue = np.arange(num_state)
        self.state_colors = np.array(state_colors)

        # Get the initial states of all nodes.
        self.initial = np.array(np.split(self.model.initial, num_state))
        # Create an array to record the state transition.
        states = np.zeros(self.N)
        # Assign each node to ...
        for i, state in enumerate(self.initial):
            # ... diff state represented by a number.
            states[state.astype(np.bool)] = i

        # Summarize the state and color of each node to an independent arr.
        self.state_list = self.state_queue[states.astype(np.int)]
        self.color_list = self.state_colors[states.astype(np.int)]

        # Assign the initial information to the node attribute.
        self.G = self.set_graph(self.state_list,
                                self.color_list,
                                G=self.G)

        self.model_name = self.model.model_name
        self.data_dir = data_dir
        self.sim = None
        self.sims = None
        self.vacci_num = 0
        self.vacci_sum = 0
        self.vacci_list = []
        self.vacci_func = None
        self.t = None
        self.propagating = []

        assert self.model.initial is not None, '[!] Initial vector is not defined.'

    def set_graph(self, state_list, color_list, G):
        # This func take a graph as input and append attributes accordingly.
        G_ = copy.deepcopy(G)
        states_ = zip(np.arange(len(state_list)), state_list)
        colors_ = zip(np.arange(len(color_list)), color_list)
        nx.set_node_attributes(G_, dict(states_), 'state')
        nx.set_node_attributes(G_, dict(colors_), 'color')
        return G_

    def select_nodes(self, A, state=0, nei_thres=0):
        # Determine those nodes in "0" state.
        opt_1 = np.where(self.state_list == state)[0]
        # Determine those nodes with neighbors.
        opt_2 = np.where(np.sum(A, axis=1) > nei_thres)[0]
        # Find only those "0" states while also connecting with the others.
        return np.intersect1d(opt_1, opt_2)

    def select_neighbors(self, A, state, nei_state, nei_thres=0):
        # Determine which nodes to pick up neighbors.
        opt = self.select_nodes(A, state)
        # Count the total neighbor numbers of the selected nodes.
        nei_nums = np.sum(A[opt, :], axis=0)
        # Determine which nodes are the neighbors of the selected nodes.
        idx_nei = np.where(nei_nums > nei_thres)[0]
        # Find out the nodes with selected state.
        idx_in = np.where(self.state_list == nei_state)[0]
        return np.intersect1d(idx_nei, idx_in)

    def disconnect_links(self, state, ratio=1):
        opt = self.select_nodes(self.A, state=3, nei_thres=0)
        # Randomly pick up certain selected nodes.
        idx = np.random.choice(opt, round(len(opt) * ratio), replace=False)
        self.A[idx, :] = 0
        self.A[:, idx] = 0

        self.G = self.set_graph(self.state_list,
                                self.color_list,
                                G=nx.Graph(self.A))

    def simulate(self, t, param=None):
        # Only the first simulation will assign self.t attribute.
        if self.t is None:
            self.t = t

        # Change the adjacent matrix using a pre-defined strategy.
        if param is not None:
            assert self.vacci_func is not None, '[!] Vacci func is not defined.'
            self.A = self.vacci_func(self.A, param)
            # Update graph if the A matrix is modified.
            self.G = self.set_graph(self.state_list,
                                    self.color_list,
                                    G=nx.Graph(self.A))

            self.vacci_sum += self.vacci_num

        # Put the model into a differential equation solver.
        self.sim = odeint(self.model.derivative, self.model.initial, t)
        # Split the results into a proper number of part.
        self.sims = np.array(np.split(self.sim, self.num_state, axis=1))
        # shape 3D: (num_state, num_period, num_node)
        return self.sims

    def reconstruct(self, t, idx, param=None):
        assert self.sims is not None, '[!] simulation is not done yet.'
        # Take a copy of "self.sim" and find the idx-th row as ...
        sim_ori = copy.deepcopy(self.sim)
        # ... the new initial probability of each node.
        self.model.initial = sim_ori[idx]

        # Change the adjacent matrix using a pre-defined strategy.
        if param is not None:
            assert self.vacci_func is not None, '[!] Vacci func is not defined.'
            self.A = self.vacci_func(self.A, param)
            # Update graph if the A matrix is modified.
            self.G = self.set_graph(self.state_list,
                                    self.color_list,
                                    G=nx.Graph(self.A))

            self.vacci_sum += self.vacci_num

        # Simulate the later part after the "A" matrix is changed.
        self.simulate(t[idx:])
        # Assign the new "A" sim result to the later part of the matrix.
        sim_ori[idx:] = self.sim
        # Let the new combined "sim_ori" as the self.sim (same as "simulate" func).
        self.sim = sim_ori
        # Split the results into a proper number of part.
        self.sims = np.array(np.split(self.sim, self.num_state, axis=1))
        # shape 3D: (num_state, num_period, num_node)
        return self.sims

    def count_neighbors(self, idx_nodes, target_state=1):
        # Create an array to record the number of targeted neighbors.
        nei_nums = np.zeros(len(idx_nodes))

        for i, nd in enumerate(idx_nodes):
            # Determine the neighbors of each selected node iteratively.
            neighbors = list(self.G.neighbors(nd))
            count = 0
            # Examine the state of those neighbor nodes.
            for nei in neighbors:
                # If the neighbor is an targeted neighbor ...
                if self.G.nodes[nei]['state'] == target_state:
                    # ... account the neighbor.
                    count += 1

            # Record the total targeted neighbors to the array.
            nei_nums[i] = count

        return nei_nums

    def state_transition(self, idx_nodes, probs, ori_s, offset):
        # Only when the probability of next state is larger than ...
        # ... a uniformly generated random number, the node state ...
        # ... will be changed.
        idx_trans = idx_nodes[probs > np.random.rand(len(probs))]
        self.state_list[idx_trans] = self.state_queue[int(ori_s + offset)]
        self.color_list[idx_trans] = self.state_colors[int(ori_s + offset)]

    def __propagate(self, idx, state_list, from_state, to_state,
                    neighbor=None, update_state=False):
        if update_state:
            # If the next state is not unique, state list should be updated.
            idx_nodes = np.where(self.state_list == from_state)[0]
        else:
            # Get the indices of those nodes at current state.
            idx_nodes = np.where(state_list == from_state)[0]

        if neighbor is not None:
            # Count how many infected neighbor nodes.
            nei_nums = self.count_neighbors(idx_nodes, target_state=neighbor)
            # Pick up only those nodes with infected neighbors.
            idx_nodes = idx_nodes[nei_nums > 0]

        # Get the probabilities of each node from the sim results.
        probs = self.sims[to_state, idx, idx_nodes]
        self.state_transition(idx_nodes, probs, from_state,
                              to_state - from_state)

    def set_propagate(self, from_state, to_state, neighbor=None, update_state=False):
        self.propagating.append([from_state, to_state,
                                 neighbor, update_state])

    def step(self, idx):
        assert self.sims is not None, '[!] simulation is not done yet.'
        # Things that changed in each "step" are only the attributes of nodes.

        # Record the original state so that the updated state will not be mixed.
        state_list = copy.deepcopy(self.state_list)

        for info in self.propagating:
            self.__propagate(idx, state_list, *info)

        self.G = self.set_graph(self.state_list,
                                self.color_list,
                                G=self.G)

    def propagation(self, reconstruct_index=[], reconstruct_param=None,
                    disconnect_index=[], disconnect_state=3):
        status = []
        probs = None

        for idx in range(len(self.t)):
            record = {}
            self.step(idx)

            if idx in reconstruct_index:
                probs = self.reconstruct(self.t, idx, param=reconstruct_param)

            if idx in disconnect_index:
                self.disconnect_links(state=disconnect_state)

            self.vacci_list.append(self.vacci_sum)

            record['iteration'] = idx
            record['state'] = copy.deepcopy(self.state_list)
            record['color'] = copy.deepcopy(self.color_list)
            record['adj'] = copy.deepcopy(self.A)
            record['graph'] = copy.deepcopy(self.G)
            record['vacci_sum'] = self.vacci_sum
            status.append(record)

        return status, probs

    def visualize(self, status, index, figsize, n_row, n_col, save=False):
        plt.figure(figsize=figsize)

        for n, idx in enumerate(index):
            plt.subplot(n_row, n_col, n + 1)
            plt.title('t = {}'.format(idx))

            nx.draw(status[idx]['graph'], self.pos,
                    with_labels=True, font_color='white',
                    font_size=6, node_size=60,
                    node_color=status[idx]['color'])
        if save:
            plt.savefig('nx.jpg')

########################################################################


if __name__ == '__main__':

    from equations import SI, SIS, SIR, SIRV

    # Seq: (The number of nodes,
    #       The number of random edges to add for each new node,
    #       Probability of adding a triangle after adding a random edge)
    G_ws = nx.watts_strogatz_graph(40, 5, 0.4)
    # G_pl = nx.powerlaw_cluster_graph(50, 5, 0.4)
    days = np.arange(0, 10, 0.1)

    # ------------------------------------------------------------------

    # SI   params: I0, beta
    # SIS  params: I0, beta, gamma
    # SIR  params: I0, R0, beta, gamma
    # SIRV params: I0, R0, V0, beta, gamma, eps

    # ------------------------------------------------------------------

    epi = EpiModel(G_ws, SIRV, num_state=4, params=[4, 2, 2, 0.4, 0.2, 0.3])
    s, i, r, v = epi.simulate(days)

    epi.set_propagate(0, 1, neighbor=1, update_state=False)
    epi.set_propagate(0, 3, neighbor=None, update_state=True)
    epi.set_propagate(1, 2, neighbor=None, update_state=False)

    status, _ = epi.propagation(disconnect_index=np.arange(16),
                                disconnect_state=3)
    # status, _ = epi.propagation(reconstruct_index=np.arange(16),
    #                             reconstruct_param=None)
    epi.visualize(status, np.arange(16), figsize=(15, 15), n_row=4, n_col=4)

    # Plot the data on three separate curves for S(t) and I(t).
    fig = plt.figure(facecolor='w')
    plt.plot(s, color='blue')
    plt.plot(i, color='red')
    plt.plot(r, color='green')
    plt.plot(v, color='orange')
    plt.xlabel('Time unit', fontsize=15)
    plt.ylabel('Probability (%)', fontsize=15)
    plt.xlim((0, 100))
    plt.ylim((0, 1))
    plt.grid(True)
    plt.legend()
    plt.show()
