'''
Created on 2013-05-19

@author: jyeung
'''


import sys
import csv

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