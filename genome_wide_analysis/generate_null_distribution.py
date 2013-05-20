'''
Created on 2013-05-19

@author: jyeung
'''

import sys
import os
from utilities import set_directories, \
chr_tools, genome_info, random_genome_locs, jplots, integrate_data

# Set directories
_cur_dir, \
_proj_dir, \
_input_dir, \
_output_dir, \
_plot_dir = set_directories.set_directories('inputs', 
                                            'outputs', 
                                            'plots')

# Set constants
# chromosome = raw_input('Insert chromosome in chr format (chr1, chrX, chrY): ')
plot_output_fname = os.path.join(_plot_dir, 'tandem_dupe_chr_distribution_null_10kb.pdf')
numb_rands_per_chr = 1000
bins_in_plot = 20000
xmin = -10000
xmax = 10000


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('TSS data must be given on the command line.')
        sys.exit()
    genome_fname = sys.argv[1]
    
    genome_path = os.path.join(_input_dir, genome_fname)
    
    chromosome_list = chr_tools.get_chr_list()
    
    all_tss = genome_info.get_all_tss_locations(genome_path, 
                                                'chrom', 
                                                'txStart', 
                                                chromosome_list)
    
    random_locs = random_genome_locs.generate_random_chr_pos(numb_rands_per_chr)
    
    tss_distances_dict, _ = integrate_data.calc_tss_dist(random_locs, all_tss, 
                                                         chromosome_list)
    
    tss_distances_list = []
    for dist_list in tss_distances_dict.values():    # Each value is a list
        tss_distances_list.extend(dist_list)
        
    jplots.plot_binned_bar_graph(tss_distances_list, bins_in_plot, 
                                 xmin,
                                 xmax,
                                 'distance from TSS', 
                                 'frequency', 
                                 'Random distribution of TSS distances '\
                                 'across genome\n%s points per chromosome' \
                                 %numb_rands_per_chr, 
                                 plot_output_fname)
    
    
    
    
    
    
    