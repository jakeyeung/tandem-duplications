'''
Created on 2013-05-17

@author: jyeung

Investigating tandem duplications in prostate cancer. 

Reads filtered tandem duplications in 445RT and place it in the context of the genome using 
a text file downloaded from UCSC. See readme.md for input file information. 
'''


import sys
import time
import os
from utilities import jplots, list_tools, genome_info, tandem_data


# Set directories
_cur_dir = os.path.dirname(os.path.realpath(__file__))
_upthree_dir = os.path.dirname(os.path.dirname(os.path.dirname(_cur_dir)))
_input_dir = os.path.join(_upthree_dir, 'inputs')
_output_dir = os.path.join(_upthree_dir, 'outputs')
_plot_dir = os.path.join(_output_dir, 'plots')


# Set constants
chromosome_list = [str('chr%s' %i) for i in range(1, 23)]
chromosome_list.append('chrX')
chromosome_list.append('chrY')
plot_output_fname = 'TSS_wide_distribution.pdf'


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print('Tandem dupes file and TSS data must be given on the command line.')
        sys.exit()
    starttime = time.time()
    
    tandem_fname = sys.argv[1]
    genome_fname = sys.argv[2]
    
    tandem_path = os.path.join(_input_dir, tandem_fname)
    genome_path = os.path.join(_input_dir, genome_fname)
    
    merged_tss_list = []
    
    for chromosome in chromosome_list:
        genome_tss = genome_info.get_tss_locations(genome_path, 'chrom', 'txStart', 
                                                   chromosome, chromosome_list)
        dist_from_tss_list = tandem_data.get_tandem_tss_distribution(tandem_path, 
                                                                     'start_1', 
                                                                     'chromosome_1', 
                                                                     genome_tss, 
                                                                     chromosome, 
                                                                     chromosome_list)
        merged_tss_list.extend(dist_from_tss_list)
    
    jplots.plot_binned_bar_graph(merged_tss_list, 100, 
                                 'Distance from a known TSS',
                                 'Frequency', 
                                 'Distribution of filtered tandem duplications by TSS', 
                                 os.path.join(_plot_dir, plot_output_fname))
    
    total = len(merged_tss_list)
    
    abs_merged_tss_list = [abs(i) for i in merged_tss_list]
    print abs_merged_tss_list
    print len(abs_merged_tss_list)
    
    chunked_list = list_tools.chunks(sorted(abs_merged_tss_list), 48)
    print chunked_list
    print len(chunked_list)
    print [[min(i), max(i)] for i in chunked_list]
    
    
    
    