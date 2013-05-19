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
import csv
from utilities import plots, list_tools


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
        
def get_tandem_tss_distribution(file_path, file_start_colname, 
                                file_chr_colname, 
                                genome_tss, chromosome, chromosome_list):
    with open(file_path, 'rb') as read_file:
        file_reader = csv.reader(read_file, delimiter='\t')
        file_colnames = file_reader.next()
        
        try: 
            file_start_index = file_colnames.index(file_start_colname)
        except ValueError:
            print('ValueError: couldnt find match %s to ' \
             'chromosome column name.' %file_start_colname)
            sys.exit()
        
        try: 
            file_chr_index = file_colnames.index(file_chr_colname)
        except ValueError:
            print('ValueError: couldnt find match %s to ' \
             'chromosome column name.' %file_chr_colname)
            sys.exit()
            
        dist_from_tss_list = []
        tss_list = []
        target_chr_index = chromosome_list.index(chromosome)
        
        for file_row in file_reader:
            
            row_chr = file_row[file_chr_index]
            row_chr_index = chromosome_list.index(row_chr)
            
            if row_chr_index < target_chr_index:
                pass    # Next row, havent reached desired chr
            
            if row_chr_index > target_chr_index:
                break    # Went past desired chr, break out of loop
            
            if row_chr_index == target_chr_index:
                dist_prev = 0
                dist_curr = 0
                tss_curr = 0
                tss_prev = 0
                tandem_curr = int(file_row[file_start_index])
                index_start = 0
                    
                count = 0
                for i in xrange(index_start, len(genome_tss)):
                    # print int(tandem_row[tandem_start_index]), genome_tss[i]
                    tss_curr = genome_tss[i]
                    dist_curr = tandem_curr - tss_curr
                    
                    if dist_curr <= 0 and count > 1:
                        # Create summary list of current and previous data
                        dist_list = [dist_curr, dist_prev]
                        abs_dist_list = [abs(i) for i in dist_list]
                        tss_list_currprev = [tss_curr, tss_prev]
                        
                        # Find minimum absolute distance and its corresponding index
                        min_abs_dist = min(abs_dist_list)
                        dist_index = abs_dist_list.index(min_abs_dist)
                        
                        # Append minimum absolute distance and at which TSS. 
                        dist_from_tss_list.append(dist_list[dist_index])
                        tss_list.append(tss_list_currprev[dist_index])
                        
                        # Print summary, debugging purposes
                        '''
                        print('%s closest TSS is %s far away at %s. ' \
                              ' i is %s' %(tandem_curr, 
                                          dist_list[dist_index], 
                                          tss_list_currprev[dist_index], 
                                          i))
                        print('tandem: %s, dist_curr: %s, dist_prev: ' \
                              '%s, tss_curr: %s, tss_prev: %s' %(tandem_curr, 
                                                                 dist_curr, 
                                                                 dist_prev, 
                                                                 tss_curr,
                                                                  tss_prev))
                        '''
                        break
                    elif i == len(genome_tss):
                        print('Could not find solution for tandem location %s, ' \
                         'farthest TSS was %s' %(tandem_curr, tss_curr))
                        sys.exit()
                    tss_prev = tss_curr
                    dist_prev = dist_curr
                    count += 1
        
    return dist_from_tss_list

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
        genome_tss = get_tss_locations(genome_path, 'chrom', 'txStart', 
                                       chromosome, chromosome_list)
        dist_from_tss_list = get_tandem_tss_distribution(tandem_path, 'start_1', 
                                                         'chromosome_1', genome_tss, chromosome, chromosome_list)
        merged_tss_list.extend(dist_from_tss_list)
    
    plots.plot_binned_bar_graph(merged_tss_list, 100, 
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
    
    
    
    