########################################################################

'''

Author: Kuo Chun-Lin
Email : guojl19@mails.tsinghua.edu.cn
'''

########################################################################

import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

from tqdm import tqdm
from epidemix.epidemic import EpiModel
from epidemix.equations import SIRV
from epidemix.vaccination import VacciModel

from epidemix.utils.plot import draw_probs, get_neighbor
from epidemix.utils.config import set_params

########################################################################

num_node = 100

# G = nx.watts_strogatz_graph(num_node, 5, 0.4)
# G = nx.powerlaw_cluster_graph(num_node, 5, 0.4)
G = nx.gnp_random_graph(num_node, 0.08)

days = np.arange(0, 50, 0.1)
colors = ['blue', 'red', 'green', 'orange']


beta_list = np.linspace(0, 1, 21)[12:13]
# for b in tqdm(beta_list, total=len(beta_list)):
for b in [0.1]:
    b = np.round(b, decimals=2)

    # [0.05, 0.10, 0.15, ..., 0.90, 0.95]
    # eps_list = np.linspace(0, 1, 21)[1:-1]
    # eps_list = np.linspace(0, 1, 6)[1:-1]
    success_rates = []
    vacc_avgs = []

    for e in [0.1]:    # -----> epsilon settings <----- #
        e = np.round(e, decimals=2)
        success_rate = []
        vacc_avg = []

        for i in range(1):   # -----> initial settings <----- #
            # for i in range(1, 3):
            I0 = int(num_node * 0.05)
            success_list = []
            vacc_list = []

            for _ in range(1):  # -----> repeat how many times <----- #

                # ======================================================
                # I0, R0, V0, beta, gamma, eps
                epi = VacciModel(G, SIRV, num_state=4, params=[
                    I0, 0, 0, b, 0.1, e],
                    vacci_strategy=None, state_colors=colors)

                probs = epi.simulate(days)
                epi.set_propagate(0, 1, neighbor=1, update_state=False)
                epi.set_propagate(0, 3, neighbor=None, update_state=True)
                epi.set_propagate(1, 2, neighbor=None, update_state=False)

                status, _ = epi.propagation()
                epi.visualize(status, np.arange(0, 64, 4),
                              figsize=(15, 15), n_row=4, n_col=4, save=True)
                draw_probs(probs, colors)
                # ======================================================

                nei_s = []
                s_trend = []

                for i in range(len(days)):
                    # Get the total number of I's neighbors which are S state.
                    idx = get_neighbor(
                        status[i]['adj'], status[i]['state'] == 1)
                    i_nei = np.sum(status[i]['state'][idx] == 0)
                    nei_s.append(i_nei)

                    # Get the total number of S under each moment.
                    s_tot = np.sum(status[i]['state'] == 0)
                    s_trend.append(s_tot)

                # Get the moment when there is no S nearby I.
                no_Snei = np.where(np.array(nei_s) == 0.0)[0][0]
                # Get the moment when there is no S in the graph.
                no_Sstat = np.where(np.array(s_trend) == 0.0)[0][0]
                # If nei_s get to 0 faster than s_trend, regard it as a success case.
                success = 1 if no_Snei < no_Sstat else 0
                success_list.append(success)

                # Count how many vacc is implemented.
                vacc_quantity = np.sum(status[no_Snei]['state'] == 3)
                vacc_list.append(vacc_quantity)

            # Get the success rate under diff init case.
            rate = np.sum(success_list) / len(success_list)
            success_rate.append(rate)

            # Get the minimum mean vaccination quantity under diff init cases.
            vacc_avg.append(np.mean(vacc_list))
            print(' >>> finished repeat <<<')

        success_rates.append(success_rate)
        vacc_avgs.append(vacc_avg)

    # df_succ = pd.DataFrame(np.array(success_rates).T, columns=eps_list)
    # df_succ.to_csv('succ_beta_{}.csv'.format(b))

    # df_vacc = pd.DataFrame(np.array(vacc_avgs).T, columns=eps_list)
    # df_vacc.to_csv('vacc_beta_{}.csv'.format(b))
