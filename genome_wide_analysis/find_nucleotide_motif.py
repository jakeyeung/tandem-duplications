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
seq_length = 60    # how long from beginning (or from end) to extract sequence
seq_rowname = 'sequence'
chr_rowname = 'chromosome_1'
start1_rowname = 'start_1'
end1_rowname = 'end_1'
start2_rowname = 'start_2'
end2_rowname = 'end_2'
start_seq_rowname = 'seq_start'
end_seq_rowname = 'seq_end'
output_path = mydirs.joinpath(outputdir, 'parsed_sequence.txt')

def read_and_write_sequences(tandem_path, output_path, 
                             colnames_list_output, colnames_list_input, 
                             seq_rowname):
    '''
    Read tandem duplication file, get sequence of duplication. Get the first 
    '''
    # Initialize output dict
    start_end_dict = {}
    start_end_dict['starts'] = []
    start_end_dict['ends'] = []
    
    # Open a write and read file and parse sequences...
    with open(tandem_path, 'rb') as tandemfile, open(output_path, 'wb') as writefile:
        tandemreader = csv.reader(tandemfile, delimiter='\t')
        filewriter = csv.writer(writefile, delimiter='\t')
        # Read tandem colnames
        tandemcolnames = tandemreader.next()
        # Write colnames to file
        cname_list = []
        for cname in colnames_list_output:
            cname_list.append(cname) 
        filewriter.writerow(cname_list)
        # Iterate rows, getting sequence and parsing it
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
            # Find first seq_length number of base pairs after breakpoint
            for i in range((breakpt_index+1), (breakpt_index+1) + seq_length):
                # Find base pairs after breakpoint
                start_seq += str(sequence[i])
                # Find base pairs counting backwards from breakpoint
                end_seq += str(sequence[(2*breakpt_index-i)])
            # Append data to lists
            start_end_dict['starts'].append(start_seq)
            start_end_dict['ends'].append(end_seq)
            # Write data to row
            coldat_list = []
            for colname in colnames_list_input:
                coldat_list.append(row[tandemcolnames.index(colname)])
            coldat_list.append(start_seq)
            coldat_list.append(end_seq)
            filewriter.writerow(coldat_list)
    return start_end_dict

def calc_base_counts(seq_list):
    '''
    Given a list of sequences, find for each position, the abundance of
    either A, C, G or T
    '''
    # Initialize dict
    bp_count_dict = {}
    bp_list = ['A', 'C', 'G', 'T']
    for base in bp_list:
        # seq_list index has all same length, choose 0 arbitrarily.
        bp_count_dict[base] = [0] * len(seq_list[0])
    
    for seq in seq_list:
        bp_index = 0
        for bp in seq:
            # Check bp is a valid nucleotide.
            if bp not in bp_list:
                print('%s not a known nucleotide sequence.' %bp)
                sys.exit()
            # Update dict
            for base in bp_list:
                if base == bp:
                    bp_count_dict[base][bp_index] += 1
                else:
                    pass
            bp_index += 1
    return bp_count_dict

def check_calc_base_counts(bp_count_dict, seq_length):
    colsum_list = []
    for i in range(0, seq_length):
        colsum = 0
        for l in bp_count_dict.values():
            colsum += l[i]
        colsum_list.append(colsum)
    return len(set(colsum_list)) <= 1
    
    
if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('Reference tandem info must be given on the command line.')
        sys.exit()
    tandem_fname = sys.argv[1]
    tandem_path = os.path.join(inputdir, tandem_fname)
    # Make list containing rownames
    colnames_list_output = [seq_rowname, chr_rowname, start1_rowname, end1_rowname, 
                            start2_rowname, end2_rowname, start_seq_rowname, end_seq_rowname]
    colnames_list_input = [seq_rowname, chr_rowname, start1_rowname, end1_rowname, 
                            start2_rowname, end2_rowname]
    start_end_dict = read_and_write_sequences(tandem_path, output_path, 
                                              colnames_list_output, colnames_list_input, seq_rowname)


    