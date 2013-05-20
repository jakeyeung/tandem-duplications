'''
Created on 2013-05-19

@author: jyeung
'''

import random
import chr_tools


def generate_random_chr_pos(samples_per_chr):
    random.seed(1)
    chr_lengths = chr_tools.get_chr_lengths('chromosome_lengths.txt', 
                                            'chromosome', 'total_length')
    chromosome = chr_tools.get_chr_list()
    
    random_locs = {}
    for c in chromosome:
        length_of_chr = chr_lengths[c]
        random_locs[c] = [random.randint(0, length_of_chr) \
                          for _ in xrange(samples_per_chr)]
    
    return random_locs
    
    
    
        