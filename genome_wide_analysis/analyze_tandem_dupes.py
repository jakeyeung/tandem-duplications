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
import matplotlib.pyplot as plt
import scipy.stats


_cur_dir = os.path.dirname(os.path.realpath(__file__))
_upthree_dir = os.path.dirname(os.path.dirname(os.path.dirname(_cur_dir)))
_plot_dir = os.path.join(_upthree_dir, 'inputs')


# Set constants
# chromosome = raw_input('Enter chromosome (1, 2, 3... X, Y): ')
chromosome = 1


def place_by_tss(tandem_path, genome_path, chromosome):
    
    full_chr = 'chr%s' %chromosome    # converts 1 to chr1, X to chrX ...     
    
    with open(tandem_path, 'rb') as tandem_file, open(genome_path, 'rb') as genome_file:
        
        tandem_reader = csv.reader(tandem_file, delimiter='\t')
        genome_reader = csv.reader(genome_file, delimiter='\t')
        
        tandem_colnames = tandem_reader.next()
        genome_colnames = genome_reader.next()
        
        # print tandem_colnames
        # print genome_colnames
        
        tandem_chr_index = tandem_colnames.index('chromosome_1')
        tandem_start_index = tandem_colnames.index('start_1')
        tandem_end_index = tandem_colnames.index('end_2')
        genome_chr_index = genome_colnames.index('chrom')
        genome_txstart_index = genome_colnames.index('txStart')
        genome_txend_index = genome_colnames.index('txEnd')
        genome_cdstart_index = genome_colnames.index('cdsStart')
        genome_cdend_index = genome_colnames.index('cdsEnd')
        genome_exoncount_index = genome_colnames.index('exonCount')
        
        genome_tss = []
        
        for genome_row in genome_reader:
            if genome_row[genome_chr_index] == full_chr:
                genome_tss.append(int(genome_row[genome_txstart_index]))
                # print int(genome_row[genome_txstart_index])
            else:
                break
        
        # print max(genome_tss)
        genome_tss = sorted(list(set(genome_tss)))
        # print max(genome_tss)
        
        dist_from_tss_list = []
        tss_list = []
        # diff_list = []
        
        for tandem_row in tandem_reader:
            dist_prev = 0
            dist_curr = 0
            tss_curr = 0
            tss_prev = 0
            for i in xrange(len(genome_tss)):
                # print int(tandem_row[tandem_start_index]), genome_tss[i]
                tandem_curr = int(tandem_row[tandem_start_index])
                tss_curr = genome_tss[i]
                dist_curr = tandem_curr - tss_curr
                if dist_curr <= 0 and i > 0:
                    # Create summary list of current and previous data
                    dist_list = [dist_curr, dist_prev]
                    abs_dist_list = [abs(i) for i in dist_list]
                    tss_list = [tss_curr, tss_prev]
                    
                    # Find minimum absolute distance and its corresponding index
                    min_abs_dist = min(abs_dist_list)
                    dist_index = abs_dist_list.index(min_abs_dist)
                    
                    # Append minimum absolute distance and at which TSS. 
                    # min_tss = tss_list[dist_index]
                    # min_dist = dist_list[dist_index]
                    dist_from_tss_list.append(dist_list[dist_index])
                    tss_list.append(tss_list[dist_index])
                    
                    # Print summary
                    print('%s closest TSS is %s far away at %s.' %(tandem_curr, dist_list[dist_index], tss_list[dist_index]))
                    print('tandem: %s, dist_curr: %s, dist_prev: %s, tss_curr: %s, tss_prev: %s' %(tandem_curr, dist_curr, dist_prev, tss_curr, tss_prev))
                    break
                # print tandem_curr, tss_curr, dist_curr, tss_prev
                tss_prev = tss_curr
                dist_prev = dist_curr
                
            raw_input('Press enter to continue')

            # plt.plot(range(len(diff_list)-1), diff_list[1:len(diff_list)])
            # plt.show()
            # break
        print dist_from_tss_list
                
            
                
        
    


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print('Tandem dupes file and TSS data must be given on the command line.')
        sys.exit()
    starttime = time.time()
    
    tandem_fname = sys.argv[1]
    genome_fname = sys.argv[2]
    
    tandem_path = os.path.join(_plot_dir, tandem_fname)
    genome_path = os.path.join(_plot_dir, genome_fname)
    
    # print _plot_dir
    # print tandem_path
    # print genome_path
    
    place_by_tss(tandem_path, genome_path, chromosome)
    
    