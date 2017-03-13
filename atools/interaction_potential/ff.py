from __future__ import division

import numpy as np
from numba import jit

class Forcefield(object):
    """ Mie-style forcefield.

    """

    def __init__(self, sigma, epsilon, n, m):
        super(Forcefield, self).__init__()

        self.sigma = sigma
        self.epsilon = epsilon
        self.n = n
        self.m = m

    @jit
    def calc_U(self, dists):
        n = self.n
        m = self.m
        epsilon = self.epsilon
        sigma = self.sigma

        C = (n/(n-m))*((n/m)**(m/(n-m)))
        '''
        sig_dist = [sigma/dist for dist in dists]
        U = sum([C*epsilon*((w**n)-(w**m)) for w in sig_dist])
        '''
        U = 0
        for dist in dists:
            w = sigma/dist
            U += C*epsilon*((w**n)-(w**m))
        return U

    def show_plot(self):
        return 0
