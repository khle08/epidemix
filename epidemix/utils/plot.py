########################################################################

'''

Author: Kuo Chun-Lin
Date  : 2020.06.06 - 2020.09.01
'''

########################################################################

import matplotlib.pyplot as plt

########################################################################


def draw_probs(probs, colors):
    # Plot the data on three separate curves for S(t) and I(t).
    fig = plt.figure(facecolor='w')

    for state, color in zip(probs, colors):
        plt.plot(state, color=color)

    plt.xlabel('Time unit', fontsize=15)
    plt.ylabel('Probability (%)', fontsize=15)
    plt.xlim((0, len(probs[0])))
    plt.ylim((0, 1))
    plt.grid(True)
    plt.legend()
    plt.show()


def get_neighbor(A, idx_node):
    return A[idx_node].sum(axis=0) > 0
