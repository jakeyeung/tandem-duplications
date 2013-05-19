'''
Created on 2013-05-19

@author: jyeung
'''

import csv
import sys


def get_tss_locations(genome_path, chr_col_name, tss_col_name, 
                      chromosome, chromosome_list):
    
    with open(genome_path, 'rb') as genome_file:
        genome_reader = csv.reader(genome_file, delimiter='\t')
        genome_colnames = genome_reader.next()
        
        try: 
            genome_chr_index = genome_colnames.index(chr_col_name)
        except ValueError:
            print('ValueError: couldnt find match %s to ' \
             'chromosome column name.' %chr_col_name)
            sys.exit()
            
        try:
            genome_txstart_index = genome_colnames.index(tss_col_name)
        except ValueError:
            print('ValueError: couldnt find match %s to ' \
             'TSS column name.' %tss_col_name)
            sys.exit()
        
        genome_tss = []
        target_chr_index = chromosome_list.index(chromosome)
        
        for genome_row in genome_reader:
            genome_chr = genome_row[genome_chr_index]    # chrX format
            
            # Get its index relative to chromosome_list. 
            file_chr_index = chromosome_list.index(genome_chr)
            
            if file_chr_index == target_chr_index:
                genome_tss.append(int(genome_row[genome_txstart_index]))
            
            elif file_chr_index < target_chr_index:
                pass
            
            elif file_chr_index > target_chr_index:
                break
                
            else:
                sys.exit('Could not find desired chromosome, exiting...')
                
        print('%s: %s known transcript start sites ' \
        'loaded, %s uniques' %(chromosome, 
                               len(genome_tss), 
                               len(set(genome_tss))))
        genome_tss = sorted(list(set(genome_tss)))
    return genome_tss