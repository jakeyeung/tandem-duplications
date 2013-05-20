'''
Created on 2013-05-19

@author: jyeung
'''

import csv
import sys


def get_tss_locations(genome_path, chr_col_name, tss_col_name, 
                      chromosome, chromosome_list):
    '''
    Return TSS locations from ONE chromosome, in list form. 
    '''
    
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


def get_all_tss_locations(genome_path, chr_col_name, tss_col_name, chromosome_list):
    '''
    Return a dictionary with chr:tsslist key:value pairs. 
    '''
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
        
        # Create empty list with key as chromosome as dictionary
        all_tss = {}
        for c in chromosome_list:
            all_tss[c] = []
        
        for genome_row in genome_reader:
            genome_chr = genome_row[genome_chr_index]    # chrX format
            genome_tss_start = genome_row[genome_txstart_index]
            if genome_chr in all_tss:
                all_tss[genome_chr].append(int(genome_tss_start))
            else:
                print('Warning, %s not found in chromosome list' %genome_chr)
                sys.exit('Exiting...')
                
        count = 0
        for l in all_tss.values():
            count += len(l)
        print('%s TSS locations found for all chromosomes.' %count)
        return all_tss
            
            
        
    
