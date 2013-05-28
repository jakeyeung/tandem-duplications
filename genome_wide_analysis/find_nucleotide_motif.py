'''
Created on 2013-05-28

@author: jyeung
'''

import os
import sys
import csv
from utilities import set_directories


# Set directories
mydirs = set_directories.my_directories()
inputdir = mydirs.joinpath(mydirs.root, 'inputs')
outputdir = mydirs.joinpath(mydirs.root, 'outputs')

# Set constants
seq_rowname = 'sequence'

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('Reference tandem info must be given on the command line.')
        sys.exit()
    tandem_fname = sys.argv[1]
    
    break_seq_start_list = []
    break_seq_end_list = []
    with open(os.path.join(inputdir, tandem_fname), 'rb') as tandemfile:
        tandemreader = csv.reader(tandemfile, delimiter='\t')
        tandemcolnames = tandemreader.next()
        for row in tandemreader:
            # Initialize empty strings
            start_seq = ''
            end_seq = ''
            sequence = row[tandemcolnames.index(seq_rowname)]    # String
            # Find | separator indicating breakpoint. 
            try:
                breakpt_index = sequence.index('|')
            except ValueError:
                print('No breakpoint found, skipping to next sequence')
            # Find first 120 base pairs after breakpoint
            for i in range((breakpt_index+1), (breakpt_index+1) + 120):
                # Find base pairs after breakpoint
                start_seq += str(sequence[i])
                # Find base pairs counting backwards from breakpoint
                end_seq += str(sequence[(2*breakpt_index-i)])
            break_seq_start_list.append(start_seq)
            break_seq_end_list.append(end_seq)
            print sequence
            print start_seq
            print end_seq
            break
        