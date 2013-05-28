'''
Created on 2013-05-27

@author: jyeung
'''


import sys
import os
import csv
from find_exon_locations import find_exon_coordinates, calc_distance_from_exon
from utilities import chr_tools, ref_and_tandem, set_directories, random_genome_locs

# Set directories
mydirs = set_directories.my_directories()
inputdir = mydirs.joinpath(mydirs.root, 'inputs')
outputdir = mydirs.joinpath(mydirs.root, 'outputs')

# Set constants
nrandopms_per_chr = 1000
# Set randomized genome loc colnames (same as actual tandem file)
tandem_fname = 'random_tandems_genome_locations.pdf'
tandem_chr_colname = 'chromosome_1'
tandem_start_colname = 'start_1'
# Set output constants
output_fname = 'exon_distance_null.txt'
output_tandem_pos_colname = 'tandem_dupe_position'
output_dist_colname = 'distance_to_exon'
output_exonstart_colname = 'closest_exon_start'
output_exonend_colname = 'closest_exon_end'
output_event_colname = 'exon_or_nonexon'


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('Reference genome info must be given on the command line.')
        sys.exit()
    ref_fname = sys.argv[1]
    
    chromosome_list = chr_tools.get_chr_list()
    
    # Create randomized tandem distribution file. 
    random_genome_locs = random_genome_locs.generate_random_chr_pos(nrandopms_per_chr)
    # Create a mock textfile with appropriate rownames
    with open(os.path.join(outputdir, tandem_fname), 'wb') as randomfile:
        randomwriter = csv.writer(randomfile, delimiter='\t')
        # Writer colnames
        header = [tandem_chr_colname, tandem_start_colname]
        randomwriter.writerow(header)
        for c in chromosome_list:
            for pos in random_genome_locs[c]:
                randomwriter.writerow([c, pos])
    
    # Initialize rf_data class, using tandem_path we just created above.
    rf_data = ref_and_tandem.data(ref_path=mydirs.joinpath(inputdir, 
                                                           ref_fname), 
                                  tandem_path=os.path.join(outputdir, tandem_fname),
                                  tandem_output_path=mydirs.joinpath(outputdir, 
                                                                     output_fname))
    
    with rf_data:
        # First initialize write header colnames.
        rf_data.writecolnames(output_tandem_pos_colname,
                              output_dist_colname, output_exonstart_colname, 
                              output_exonend_colname, output_event_colname)
        # Supply its first row
        # print rf_data.refheaders
        # print rf_data.tandemheaders
        
        reffirstrow = rf_data.refnext()
        tandemfirstrow = rf_data.tandemnext()
        
        distances_dict = {}
        coordinates_dict = {}
        locations_dict = {}
        for chromosome in chromosome_list:
            # Initialize empty list
            # Input first row of chr1, get firstrow of chr2 back and exon_dict
            # containing exonStart and End positions of chr1. Repeat for all
            # chromosomes. 
            exon_coordinates, \
            reffirstrow = find_exon_coordinates(rf_data, 
                                                chromosome, 
                                                reffirstrow)
            # From tandem dupe, calculate distance from nearest exon. 
            distances_list, \
            coordinates_list, \
            locations_list, \
            tandemfirstrow = calc_distance_from_exon(rf_data,
                                                     exon_coordinates, 
                                                     chromosome, 
                                                     tandemfirstrow)
            distances_dict[chromosome] = distances_list
            coordinates_dict[chromosome] = coordinates_list
            locations_dict[chromosome] = locations_list
            print('Done for %s' %chromosome)