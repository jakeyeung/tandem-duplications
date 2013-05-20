'''
Created on 2013-05-19

@author: jyeung
'''


def calc_tss_dist(all_tandems, all_tss, chromosome_list):
    '''
    Return distance of tandem duplication from TSS, grouped by chromosomes
    '''
    # Initialize dictionary
    tss_distances = {}
    closest_tss_locs = {}
    
    # Loop through chromosomes
    for c in chromosome_list:
        tss_distances[c] = []
        closest_tss_locs[c] = []
        chr_tss = sorted(all_tss[c])
        
        for curr_tandem in sorted(all_tandems[c]):
            # Initialize variables
            tss_count = 0
            prev_dist = 0
            prev_tandem = 0
            for tss in chr_tss:
                curr_dist = curr_tandem - tss
                if curr_dist < 0 and tss_count > 0:
                    # Find whether current dist or prev dist is closer
                    curr_and_prev = [curr_dist, prev_dist]
                    abs_curr_and_prev = [abs(i) for i in curr_and_prev]
                    min_abs_dist = min(abs(i) for i in curr_and_prev)
                    
                    # Find whether prev or current tandem was closer
                    min_index = abs_curr_and_prev.index(min_abs_dist)
                    closest_tss = [curr_tandem, prev_tandem][min_index]
                    min_dist = curr_and_prev[min_index]
                    
                    # Append dictionaries/list
                    tss_distances[c].append(min_dist)
                    closest_tss_locs[c].append(closest_tss)
                    
                    # Min dist found, move to next tandem duplication
                    break
                
                tss_count += 1
                prev_dist = curr_dist
                prev_tandem = curr_tandem
    return tss_distances, closest_tss_locs


    