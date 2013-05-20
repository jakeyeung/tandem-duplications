'''
Created on 2013-05-19

@author: jyeung
'''

import set_directories
import csv
import sys
import os


def get_chr_list():
    '''
    Generate a chromosome list chr1, chr2... chrX, chrY
    '''
    chromosome_list = [str('chr%s' %i) for i in range(1, 23)]
    chromosome_list.append('chrX')
    chromosome_list.append('chrY')
    return chromosome_list

def get_chr_lengths(chr_length_fname, chrname_colname, chrlength_colname):
    '''
    Get chromosome lengths from a text file. 
    '''
    _input_dir = set_directories.set_input_dir('inputs')
    chr_lengths_full_path = os.path.join(_input_dir, chr_length_fname)
    
    with open(chr_lengths_full_path, 'rb') as chrlen_file:
        chrlen_reader = csv.reader(chrlen_file, delimiter='\t')
        chrlen_colnames = chrlen_reader.next()
        try: 
            chrname_index = chrlen_colnames.index(chrname_colname)
        except ValueError:
            print('ValueError: couldnt find match %s to ' \
             'chromosome column name.' %chrname_colname)
            sys.exit()
            
        try:
            chrlength_index = chrlen_colnames.index(chrlength_colname)
        except ValueError:
            print('ValueError: couldnt find match %s to ' \
             'TSS column name.' %chrlength_colname)
            sys.exit()
        
        chr_len_dic = {}
        for row in chrlen_reader:
            chromosome = row[chrname_index]
            chr_len_dic[chromosome] = int(row[chrlength_index])
        
        return chr_len_dic