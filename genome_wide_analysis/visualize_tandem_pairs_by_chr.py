'''
Created on 2013-05-19

@author: jyeung
'''

import sys
import time
import os
from utilities import genome_info, tandem_data, jplots


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
# plot_output_fname = 'tandem_dupe_chr_distribution.pdf'


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print('Tandem dupes file and TSS data must be given on the command line.')
        sys.exit()
    starttime = time.time()
    
    tandem_fname = sys.argv[1]
    genome_fname = sys.argv[2]
    
    tandem_path = os.path.join(_input_dir, tandem_fname)
    genome_path = os.path.join(_input_dir, genome_fname)
    
    all_tss = genome_info.get_all_tss_locations(genome_path, 'chrom', 'txStart', chromosome_list)
    all_tandem_dupes_pairs = tandem_data.tandem_dupes_all_pairs(tandem_path,
                                                                 'start_1', 
                                                                 'end_1', 
                                                                 'chromosome_1', 
                                                                 chromosome_list)
    
    for chromosome in chromosome_list:
        full_save_path = os.path.join(_plot_dir, 
                                      'tandem_dupe_chr_distribution_%s.pdf' %chromosome)
        jplots.plot_tandem_tss(all_tandem_dupes_pairs, all_tss, chromosome, full_save_path)
    