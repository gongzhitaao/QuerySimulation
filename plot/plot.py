#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

if __name__ == '__main__':
    data = np.genfromtxt('hitrate.txt')

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(data[:, 0], data[:, 1])

    ax.set_xlabel('query window size')
    ax.set_ylabel('hit rate')

    plt.show()
