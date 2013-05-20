'''
Created on 2013-05-19

@author: jyeung
'''

import sys
import os
from utilities import tandem_data, set_directories, jplots

# Set directories
_cur_dir, \
_proj_dir, \
_input_dir, \
_output_dir, \
_plot_dir = set_directories.set_directories('inputs', 
                                            'outputs', 
                                            'plots')

# Set constants
save_path = os.path.join(_plot_dir, 'gc_content_of_dupes.pdf')

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('Tandem data must be given on the command line.')
        sys.exit()
    tandem_fname = sys.argv[1]
    tandem_path = os.path.join(_input_dir, tandem_fname)
    
    gc_chr_pairlist = tandem_data.calculate_tandem_gc_content(tandem_path, 'sequence', 'chromosome_1')
    
    just_gc = [i[0] for i in gc_chr_pairlist]
    
    jplots.plot_binned_bar_graph(just_gc, 30, 0, 1, 'gc_content', 'frequency', 'gc_content_across_tandem_dupes', save_path)
    
    
    
    