'''
Created on 2013-05-19

@author: jyeung
Various plotting methods
'''

import matplotlib.pyplot as plt
# from scipy.stats import gaussian_kde
import numpy as np
import sys


def plot_binned_bar_graph(xlist, bins, xmin, xmax, xlabel, 
                          ylabel, title, full_save_path, autoxlim=False):
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
    if autoxlim==False:
        plt.xlim([xmin, xmax])
    elif autoxlim==True:
        pass
    else:
        print('autoxlim is neither true nor false')
        sys.exit()
    plt.savefig(full_save_path)
    
def plot_tandem_tss(tandem_list_pairs, tss_list, chromosome, full_save_path):
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
    
    print('Data crunched, plotting for %s' %chromosome)
    plt.savefig(full_save_path)
    
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
    
def plot_vertical_lines_genome_wide(tandem_pairs_dict, chr_list, 
                                    chr_lengths_dic,
                                    save_output_path):
    '''
    Visualize tandem duplicates on the genome scale. 
    '''
    plt.figure()
    # Create x-axis spanning whole genome.
    # Initialize list containing base pair at which new chr arises
    chr_breaks = []
    current_bp = 0
    chr_breaks.append(current_bp)
    for c in chr_list:
        current_bp += chr_lengths_dic[c]
        chr_breaks.append(current_bp)
    
    # Initialize plot
    plt.axis((0, chr_breaks[-1], 0, 1))
    # Add dashed lines at chromosome separations
    for chrstart in chr_breaks:
        plt.axvline(x=chrstart, color='black', ls='dashed')
    # Add texts to label chromosomes and convert tandem-dic to list
    i = 0
    tandem_pairs_list = []
    for c in chr_list:
        # Add text to label chromosomes on figure
        xpos = float(chr_breaks[i] + chr_breaks[i+1])/2 - 30000000    # Shifted to left
        plt.text(xpos, 1, c)
        # Convert tandem_dic to list of tuple pairs
        for tup in tandem_pairs_dict[c]:
            # Add chr_breaks[i+1] to each element in tuple so when you plot
            # the tandem duplication on the genome, it would be shifted correctly
            # in the chromosome. Because tandem duplication location is 
            # recorded as a position relative to chromosome. If we plot the whole
            # genome, we must shift all the positions by where the chromosome
            # is located in the genome. 
                tandem_pairs_list.append(tuple([(tandempos + chr_breaks[i]) for tandempos in tup]))
        i += 1
        
        
    # Add tandem locations
    # Convert tandem_dic to list
    for pairs in tandem_pairs_list:
        plt.axvspan(pairs[0], pairs[1], alpha=0.25, color='blue')
    plt.show()
        
    
    
    
    