'''
Created on 2013-05-26

@author: jyeung

Take tandem pairs and plot them across the genome. 
'''

import sys
import os
from utilities import chr_tools, jplots, set_directories, tandem_data


# Initialize directories
mydirs = set_directories.my_directories()

# Set constants
inputfolder = 'inputs'
outputfolder = 'outputs'
plotfolder = 'plots'

chr_length_fname = 'chromosome_lengths.txt'
chrname_colname = 'chromosome'    # colname in the fname
chrlength_colname = 'total_length'    # colname in the fname

save_output_fname = 'genome_wide_distribution_tandem_dupes.pdf'
save_output_path = os.path.join(mydirs.root, outputfolder, plotfolder, 
                                save_output_fname)


# Initialize chromosome information
chromosome_list = chr_tools.get_chr_list()

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('Tandem data must be given on the command line.')
        sys.exit()
    tandem_fname = sys.argv[1]
    tandem_path = os.path.join(mydirs.root, inputfolder, tandem_fname)
    
    # Calculate chromosome lengths
    chr_length_dic = chr_tools.get_chr_lengths(chr_length_fname, 
                                               chrname_colname, chrlength_colname)
    
    # Get tandem_pairs
    all_tandem_dupes_pairs_dic = tandem_data.tandem_dupes_all_pairs(tandem_path,
                                                                    'start_1', 
                                                                    'end_1', 
                                                                    'chromosome_1',
                                                                    chromosome_list)
        
    jplots.plot_vertical_lines_genome_wide(all_tandem_dupes_pairs_dic, 
                                           chromosome_list, 
                                           chr_length_dic, 
                                           save_output_path)