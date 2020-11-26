
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

from tqdm import tqdm
from epidemix.epidemic import EpiModel
from epidemix.equations import SIRVir1
from epidemix.vaccination import VacciModel

from epidemix.utils.plot import draw_probs
from epidemix.utils.config import set_params

G = nx.watts_strogatz_graph(40, 5, 0.4)
days = np.arange(0, 50, 0.1)
colors = ['blue', 'red', 'green', 'orange', 'purple', 'cyan']

# I0, R0, V0, i0, r0, beta1, gamma1, alpha1, beta2, gamma2, alpha2, eps
epi = VacciModel(G, SIRVir1, num_state=6, params=[
    4, 2, 0, 0, 0, 0.3, 0.1, 0.01, 0.1, 0.6, 0.01, 0.1],
    vacci_strategy=None, state_colors=colors)

probs = epi.simulate(days)
epi.set_propagate(0, 1, neighbor=1, update_state=False)
epi.set_propagate(0, 3, neighbor=None, update_state=True)
epi.set_propagate(1, 2, neighbor=None, update_state=False)
epi.set_propagate(2, 0, neighbor=None, update_state=False)

epi.set_propagate(3, 4, neighbor=1, update_state=False)
epi.set_propagate(4, 5, neighbor=None, update_state=False)
epi.set_propagate(5, 0, neighbor=None, update_state=False)

status, _ = epi.propagation()
epi.visualize(status, np.arange(16), figsize=(15, 15), n_row=4, n_col=4)
draw_probs(probs, colors)
