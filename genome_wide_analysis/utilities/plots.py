'''
Created on 2013-05-19

@author: jyeung
Various plotting methods
'''

import matplotlib.pyplot as plt
# from scipy.stats import gaussian_kde
import numpy as np


def plot_binned_bar_graph(xlist, bins, xlabel, ylabel, title, full_save_path):
    hist, bins = np.histogram(xlist, bins=bins)
    width = 0.7*(bins[1]-bins[0])
    center = (bins[:-1]+bins[1:])/2
    plt.bar(center, hist, align = 'center', width = width)
    
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    # plt.xlim([-50000, 50000])
    plt.savefig(full_save_path)