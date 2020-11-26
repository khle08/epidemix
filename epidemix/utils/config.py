########################################################################

'''

Author: Kuo Chun-Lin
Date  : 2020.06.06 - 2020.09.01
'''

########################################################################

import argparse

########################################################################


def set_params():
    parser = argparse.ArgumentParser()
    parser.add_argument('--model', type=str, default='SIR',
                        help='model.yaml path')
    parser.add_argument('--net', type=str, default='poisson',
                        help='model.yaml path')
    parser.add_argument('--vacci', type=str, default='acquaintance',
                        help='model.yaml path')
    parser.add_argument('--n', type=int, default=50,
                        help='num of nodes')
    parser.add_argument('--m', type=int, default=5,
                        help='num of random edges to add for each new node')
    parser.add_argument('--p', type=float, default=0.2,
                        help='Probability of adding a triangle after adding a random edge')
    parser.add_argument('--i0', type=int, default=8,
                        help='Initial infected number of people')
    parser.add_argument('--r0', type=int, default=0,
                        help='Initial recovered number of people')
    parser.add_argument('--beta', type=float, default=0.9,
                        help='transimmision rate on a contact')
    parser.add_argument('--gamma', type=float, default=0.2,
                        help='recovery rate')
    parser.add_argument('--max-day', type=int, default=100, help='')
    parser.add_argument('--vacci-start', type=int, default=5, help='')
    parser.add_argument('--vacci-end', type=int, default=95, help='')
    parser.add_argument('--vacci-freq', type=int, default=5, help='')
    parser.add_argument('--run-time', type=int, default=100, help='')
    return parser.parse_args()
