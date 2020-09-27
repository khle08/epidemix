########################################################################

'''
>>> attribute clustering:
"eva", "ilouvain"

>>> biparitte clustering:
"bimlpa"

>>> crisp partition
"louvain", "leiden", "rb_pots", "rber_pots", "cpm", "significance_communities", 
"surprise_communities", "greedy_modularity", "der", "label_propagation",
"async_fluid", "infomap", "walktrap", "girvan_newman", "em", "scan",
"gdmp2", "spinglass", "eigenvector", "agdl", "frc_fgsn", "sbm_dl",
"sbm_dl_nested", "markov_clustering", "edmot", "chinesewhispers"

>>> edge clustering
"hierarchical_link_community"

>>> overlapping partition
"ego_networks", "demon", "angel", "node_perception", "overlapping_seed_set_expansion",
"kclique", "lfm", "lais2", "congo", "conga", "lemon", "slpa", "multicom",
"big_clam", "danmf", "egonet_splitter", "nnsed", "nmnf", "aslpaw", "percomvc"

Author: Kuo Chun-Lin
Date  : 2020.06.06 - 2020.06.17
'''

########################################################################

import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

from cdlib import algorithms as algo
# import community
# import partition_networkx
# from networkx.algorithms import community as nxc

########################################################################


def count_neighbors(communities, A):
    out_com, in_com = [], []
    for nodes in communities:
        # Get the neighbor indices of each node.
        idx_nei = np.vstack(np.where(A[nodes] == 1)).T

        out_count, in_count = [], []
        for i in range(len(nodes)):
            # Get the actual neighbor nodes.
            neighbors = idx_nei[idx_nei[:, 0] == i][:, 1]
            # Find the number of links connecting to the other com.
            out_count.append(len(np.setdiff1d(neighbors, nodes)))
            # Find the number of links connecting in the com.
            in_count.append(len(np.intersect1d(neighbors, nodes)))

        out_com.append(out_count)
        in_com.append(in_count)

    return in_com, out_com


def formulate(lists):
    # Get how long is the list at each round.
    quantity = [len(i) for i in lists]
    # Create an empty "-1" matrix.
    table = -1 * np.ones((len(quantity), max(quantity)))

    for i, statistic in enumerate(lists):
        # Save the data to the created matrix.
        table[i, :len(statistic)] = statistic

    return table


def max_rank(communities, nei_statistic):
    # Formulate the list of lists into a 2D matrix.
    communities = formulate(communities)
    nei_statistic = formulate(nei_statistic)

    ranks = []
    # Iterate from max to min links.
    for val in np.unique(nei_statistic)[::-1][:-1]:
        # Get the actual node ID according to the num of link.
        rank = communities[np.where(nei_statistic == val)]
        rank = rank.astype(np.int).tolist()
        ranks.append(rank)

    return ranks


def get_amount(node_list, num=5):
    total = []
    for i in node_list:
        # Merge all the nodes into a vector.
        total.extend(i)

    while len(total) > 0:
        # Get the wanted number of nodes.
        yield list(total[:num])

        try:
            # Delete the obtained nodes.
            total = np.delete(total, np.arange(num))
        except:
            total = []


def target_nodes(communities, num, out_community=True):
    # Count the number of in / out neighbor number of each node.
    in_com, out_com = count_neighbors(communities, A)
    # Decide whether in or out statistic is used.
    statistic = out_com if out_community else in_com
    # Find the corresponding nodes according to the neighbor numbers.
    node_list = max_rank(communities, statistic)
    # Return the wanted number of node list.
    return list(get_amount(node_list, num=num))


########################################################################


if __name__ == '__main__':
    G = nx.watts_strogatz_graph(100, 5, 0.4)
    A = np.array(nx.adj_matrix(G).todense())
    pos = nx.spring_layout(G)

    # plt.figure()
    # nx.draw(G, pos, with_labels=True,
    #         font_color='white', node_color='blue',
    #         font_size=6, node_size=100)
    # plt.show()

    communities = algo.louvain(nx.Graph(A)).communities
    opt = target_nodes(communities, 5, out_community=True)

    print(opt)
