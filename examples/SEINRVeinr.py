
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

from tqdm import tqdm
from epidemic import EpiModel
from equations import SEINRVeinr
from vaccination import VacciModel

from utils.plot import draw_probs
from utils.config import set_params

G = nx.watts_strogatz_graph(80, 5, 0.4)
days = np.arange(0, 50, 0.1)
colors = ['blue', 'yellow', 'red', 'gray', 'green',
          'orange', 'wheat', 'purple', 'lightgray', 'cyan']

# E0, I0, N0, R0, V0, e0, i0, n0, r0,
# beta1, theta1, eta1, gamma1, pi1, alpha1, eps,
# beta2, theta2, eta2, gamma2, pi2, alpha2
epi = VacciModel(G, SEINRVeinr, num_state=10, params=[
    4, 2, 2, 1, 1, 1, 1, 1, 1,
    0.3, 0.7, 0.8, 0.2, 0.4, 0.1, 0.1,
    0.6, 0.2, 0.1, 0.6, 0.7, 0.1],
    vacci_strategy=None, state_colors=colors)

# S, E, I, N, R, V, e, i, n, r
# 0, 1, 2, 3, 4, 5, 6, 7, 8, 9
probs = epi.simulate(days)
epi.set_propagate(0, 1, neighbor=2, update_state=False)
epi.set_propagate(0, 5, neighbor=None, update_state=True)
epi.set_propagate(1, 2, neighbor=None, update_state=False)
epi.set_propagate(1, 3, neighbor=None, update_state=True)
epi.set_propagate(2, 4, neighbor=None, update_state=False)
epi.set_propagate(3, 4, neighbor=None, update_state=False)
epi.set_propagate(4, 0, neighbor=None, update_state=False)

epi.set_propagate(5, 6, neighbor=2, update_state=False)
epi.set_propagate(6, 7, neighbor=None, update_state=False)
epi.set_propagate(6, 8, neighbor=None, update_state=True)
epi.set_propagate(7, 9, neighbor=None, update_state=False)
epi.set_propagate(8, 9, neighbor=None, update_state=False)
epi.set_propagate(9, 0, neighbor=None, update_state=False)

status, _ = epi.propagation()
epi.visualize(status, np.arange(0, 64, 4),
              figsize=(15, 15), n_row=4, n_col=4)
draw_probs(probs, colors)
