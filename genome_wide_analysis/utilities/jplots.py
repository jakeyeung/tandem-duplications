'''
Created on 2013-05-19

@author: jyeung
Various plotting methods
'''

import matplotlib.pyplot as plt
# from scipy.stats import gaussian_kde
import numpy as np


def plot_binned_bar_graph(xlist, bins, xmin, xmax, xlabel, 
                          ylabel, title, full_save_path):
    '''
    Plot almost like a histogram but bar graph, user specify bins
    '''
    hist, bins = np.histogram(xlist, bins=bins)
    width = 0.7*(bins[1]-bins[0])
    center = (bins[:-1]+bins[1:])/2
    plt.bar(center, hist, align = 'center', width = width)
    
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.xlim([xmin, xmax])
    plt.savefig(full_save_path)
    
def plot_tandem_tss(tandem_list_pairs, tss_list, chromosome):
    '''
    For a chromosome, plot TSS_list as vertical lines and 
    tandem_list_pairs as vertical bars (start and end pos)
    '''
    chr_tss_list = [int(i) for i in tss_list[chromosome]]
    chr_tandem_list_pairs = [(int(i[0]), int(i[1])) for i in tandem_list_pairs[chromosome]]
    
    plt.figure()
    plt.axis((min(chr_tss_list), max(chr_tss_list), 0, 1))
    plt.xlabel('Chromosome location')
    plt.title('Chromosome %s TSS locations in dots, ' \
              'tandem duplications in vertical lines' %chromosome)
    
    plt.plot(chr_tss_list, [0.5]*(len(chr_tss_list)), 'ro')
    
    for pairs in chr_tandem_list_pairs:
        plt.axvspan(pairs[0], pairs[1], facecolor='g', alpha=0.5)
    
    print('Data crunched, plotting...')
    plt.show()
    
def plot_tandem_tss_nopairs(tandem_list, tss_list, chromosome):
    '''
    For a chromosome, plot TSS_list as vertical lines and 
    tandem_list as vertical bars (start and end pos)
    '''
    chr_tss_list = [int(i) for i in tss_list[chromosome]]
    chr_tandem_list = [i for i in tandem_list[chromosome]]
    
    plt.figure()
    plt.axis((min(chr_tss_list), max(chr_tss_list), 0, 1))
    plt.xlabel('Chromosome location')
    plt.title('Chromosome %s TSS locations in dots, ' \
              'tandem duplications in vertical lines' %chromosome)
    
    plt.plot(chr_tss_list, [0.5]*(len(chr_tss_list)), 'ro')
    
    for pos in chr_tandem_list:
        plt.axvline(pos)
    
    print('Data crunched, plotting...')
    plt.show()
        
    
    