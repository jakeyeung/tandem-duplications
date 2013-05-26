'''
Created on 2013-05-19

@author: jyeung
'''

import sys
import time
import os
from utilities import genome_info, tandem_data, jplots, integrate_data


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
plot_output_fname = os.path.join(_plot_dir, 
                                 'tandem_dupe_chr_distribution_cdsStart_10kbxrange.pdf')
text_output_path = os.path.join(_output_dir, 'distance_to_tss.txt')
tandem_position_colname = 'start_1'
tandem_chrcolname = 'chromosome_1'
genomedat_chrcolname = 'chrom'
genomedat_startcolname = 'cdsEnd'
bins_in_plot = 1800   # 2000 for xlim[-10000, 10000], 100 if no xlims
xmin = -10000
xmax = 10000
autoxlim = False    # True or False

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print('Tandem dupes file and TSS data must be given on the command line.')
        sys.exit()
    starttime = time.time()
    
    tandem_fname = sys.argv[1]
    genome_fname = sys.argv[2]
    
    tandem_path = os.path.join(_input_dir, tandem_fname)
    genome_path = os.path.join(_input_dir, genome_fname)
    
    all_tss = genome_info.get_all_tss_locations(genome_path, 
                                                genomedat_chrcolname, 
                                                genomedat_startcolname, 
                                                chromosome_list)
    
    all_tandem_dupes = tandem_data.tandem_dupes_all(tandem_path, 
                                                    tandem_position_colname, 
                                                    tandem_chrcolname, 
                                                    all_tss, 
                                                    chromosome_list)
    
    
    min_dists, closest_tss = integrate_data.calc_tss_dist(all_tandem_dupes, 
                                                          all_tss, 
                                                          chromosome_list, 
                                                          text_output_path)
    
    min_dists_list = []
    for min_dists_chr_list in min_dists.values():
        min_dists_list.extend(min_dists_chr_list)
    
    lessthan2kb, btwn2kb10kb, grtrthan10kb = \
    integrate_data.bin_distances(min_dists_list)
    
    # print len(lessthan2kb), len(btwn2kb10kb), len(grtrthan10kb)

    jplots.plot_binned_bar_graph(min_dists_list, bins_in_plot,
                                 xmin,
                                 xmax,
                                 'Min Distance to %s' %genomedat_startcolname, 
                                 'Frequency', 
                                 'Tandem Duplications by Min Dist to %s' %genomedat_startcolname, 
                                 plot_output_fname, 
                                 autoxlim=autoxlim)
    
    largest_tss_dic = {}
    for k, l in all_tandem_dupes.iteritems():
        largest_tss_dic[k] = max(l)
    
    print 'Largest distances by chromosome: '
    print largest_tss_dic
    
    
    
    