'''
Created on 2013-05-25

@author: jyeung

FIXME and POTENTIAL BUGS:
If you have no tandem duplication event OR if you have no exon coordinates in
your reference data, weird things may happen. For example, the StopIteration
error handler would be unable to assign the variable lastrow to anything! 
'''


import sys
from utilities import set_directories, ref_and_tandem, chr_tools, jplots, list_tools


# Set directories
mydirs = set_directories.my_directories()
inputdir = mydirs.joinpath(mydirs.root, 'inputs')
outputdir = mydirs.joinpath(mydirs.root, 'outputs')
plotsavedir = mydirs.joinpath(outputdir, 'plots', 'exon_distance')
save_path_name = mydirs.joinpath(plotsavedir, 'exon_distances.pdf')

# Define constants
ref_chr_colname = 'chrom'
ref_exonstarts_colname = 'exonStarts'
ref_exonends_colname = 'exonEnds'
tandem_chr_colname = 'chromosome_1'    # Equal to 'chromosome_2'
tandem_start_colname = 'start_1'
tandem_end_colname = 'end_2'
output_tandem_pos_colname = 'tandem_dupe_position'
output_dist_colname = 'distance_to_exon'
output_exonstart_colname = 'closest_exon_start'
output_exonend_colname = 'closest_exon_end'
output_event_colname = 'exon_or_nonexon'

def find_exon_coordinates(rf_data, chromosome, firstrow):
    '''
    Find exon coordinates, start and stop locations in a tuple for 
    a specific chromosome.  
    '''
    # Initialze list
    exon_list = []
    
    # Process first row
    cur_chr = firstrow[rf_data.refheaders.index(ref_chr_colname)]
    if cur_chr == chromosome:
        exonstart_str = firstrow[rf_data.refheaders.index(ref_exonstarts_colname)]
        exonend_str = firstrow[rf_data.refheaders.index(ref_exonends_colname)]
        # If last character is a comma, remove it. 
        if exonstart_str[-1] == ',':
            exonstart_str = exonstart_str[:-1]
        if exonend_str[-1] == ',':
            exonend_str = exonend_str[:-1]
        # Append to list
        exonstart_list = [int(i) for i in exonstart_str.split(',')]
        exonend_list = [int(i) for i in exonend_str.split(',')]
        
        # Check their lengths are equal, if so, form tuples
        # between start and end. 
        if len(exonstart_list) == len(exonend_list):
            for i in range(len(exonstart_list)):
                exon_list.append((exonstart_list[i], exonend_list[i]))
        else:
            print('Error: Exon start and end list lengths not equal.')
            sys.exit()
    else:
        print('First row chr, %s, does not match %s' %(cur_chr, chromosome))
        sys.exit()
    
    # Process all other rows until the row chromosome does not match
    # the chromosome of interest. 
    while True:
        try:
            row = rf_data.refnext()
            cur_chr = row[rf_data.refheaders.index(ref_chr_colname)]
            if cur_chr == chromosome:
                exonstart_str = row[rf_data.refheaders.index(ref_exonstarts_colname)]
                exonend_str = row[rf_data.refheaders.index(ref_exonends_colname)]
                # If last character is a comma, remove it. 
                if exonstart_str[-1] == ',':
                    exonstart_str = exonstart_str[:-1]
                if exonend_str[-1] == ',':
                    exonend_str = exonend_str[:-1]
                # Append to list
                exonstart_list = [int(i) for i in exonstart_str.split(',')]
                exonend_list = [int(i) for i in exonend_str.split(',')]
                # Check their lengths are equal, if so, form tuples
                # between start and end. 
                if len(exonstart_list) == len(exonend_list):
                    for i in range(len(exonstart_list)):
                        exon_list.append((exonstart_list[i], exonend_list[i]))
                else:
                    print('Error: Exon start and end list lengths not equal.')
                    sys.exit()
            else:
                print('Current refrow is now in %s, breaking...' %cur_chr)
                lastrow = row
                break
        except StopIteration:
            print('Ref row count: %s. No more rows to iterate.' %rf_data.refrowcount)
            lastrow = row
            break
    
    # Sort and set exon_list
    return sorted(set(exon_list)), lastrow

def calc_distance_from_exon(exon_coordinates, chromosome, firstrow):
    '''
    From exon_coordinates, find shortest distance to an exon. Use
    tandem start only...
    '''
    # Initialize distances
    distances_list = []
    coordinates_list = []
    locations_list = []
    
    # Calculate for first row...
    cur_chr = firstrow[rf_data.tandemheaders.index(tandem_chr_colname)]
    # If not equal to chromosome, error exit.
    if cur_chr == chromosome:
        tandem_pos = int(firstrow[rf_data.tandemheaders.index(tandem_start_colname)])
    else:
        sys.exit('First row chr, %s, not equal to chromosome.')
    
    # First first row, find distance, coordinate, location.
    distance, coordinate, location = check_exon_or_not(tandem_pos, exon_coordinates)
    
    # Append to list
    distances_list.append(distance)
    coordinates_list.append(coordinate)
    locations_list.append(location)
    
    # Loop for all other tandem dupes, until either you run out of rows or 
    # current chromosome does not equal to inputed chromosome.
    
    while True:
        try:
            row = rf_data.tandemnext()
            cur_chr = row[rf_data.tandemheaders.index(tandem_chr_colname)]
            tandem_pos = int(row[rf_data.tandemheaders.index(tandem_start_colname)])
            if cur_chr == chromosome:
                # Get distance, coordinate and whether its exon or nonexon
                distance, coordinate, location = check_exon_or_not(tandem_pos, exon_coordinates)
                distances_list.append(distance)
                coordinates_list.append(coordinate)
                locations_list.append(location)
                # Append results to current row
                # If you append more things in the row, ensure that 
                # you've written enough colnames in your headers to match
                # the output. 
                row.append(tandem_pos)
                row.append(distance)
                row.append(coordinate[0])    # Exon start
                row.append(coordinate[1])    # Exon end
                row.append(location)    # Exon or nonexon
                # Write appended row to file
                # Check if row length matches writefile header
                if len(row) == len(rf_data.outheaders):
                    rf_data.writerow(row)
                else: 
                    sys.exit('Error: length of row to write '\
                             'does not match length of colnames')
            else:
                print('Current row chr is %s, breaking.' %cur_chr)
                lastrow = row
                break
        except StopIteration:
            print('Tandem row count %s. No more rows to iterate.' %rf_data.tandemrowcount)
            try:
                lastrow = row
            except UnboundLocalError:
                print('Only one tandem event in %s, returning first row' %chromosome)
                lastrow = firstrow
            break
    return distances_list, coordinates_list, locations_list, lastrow

def check_exon_or_not(tandem_pos, exon_coordinates):
    # Calculate distance for first exon, then loop. Dist = tandem - exonstart
    prev_exonstart =  exon_coordinates[0][0]
    prev_exonend = exon_coordinates[0][1]
    
    # Check if it is between exon
    if tandem_pos > prev_exonstart and tandem_pos < prev_exonend:
        location = 'exon'
        distance = 0
        coordinate = exon_coordinates[0]
        raw_input('%s is between %s' %(tandem_pos, exon_coordinates[0]))
    
    # Loop for all other distances
    # List index, starts at 1 because we did first exon already
    for i in range(1, len(exon_coordinates)):
        cur_exonstart = exon_coordinates[i][0]
        cur_exonend = exon_coordinates[i][1]
        
        if tandem_pos > cur_exonstart and tandem_pos < cur_exonend:
            location = 'exon'
            distance = 0
            coordinate = exon_coordinates[0]
            break
        
        elif tandem_pos > prev_exonend and tandem_pos < cur_exonstart:
            location = 'nonexon'
            # Calculate shortest distance, is it prev or cur exon?
            dist_prev = exon_coordinates[i-1][1] - tandem_pos
            dist_cur = exon_coordinates[i][0] - tandem_pos
            if abs(dist_cur) < abs(dist_prev):    # Choose cur exon
                coordinate = exon_coordinates[i]
                distance = dist_cur
            elif abs(dist_cur) > abs(dist_prev):    # Choose prev
                coordinate = exon_coordinates[i-1]
                distance = dist_prev
            else:
                raw_input('Distances are equal, enter to put cur exon as shortest.')
                coordinate = exon_coordinates[i]
                distance = dist_cur
            break
        else:    # Haven't found closest exon yet, keep looping.
            prev_exonstart = cur_exonstart
            prev_exonend = cur_exonend
    else:
        # Loop fell without finding an exon or a closest exon.
        sys.exit('Loop fell without finding an exon \n %s \n %s' %(tandem_pos, exon_coordinates))
    return distance, coordinate, location
    
    
if __name__ == '__main__':
    if len(sys.argv) < 3:
        print('Reference genome info and tandem dupes file must be given on the command line.')
        sys.exit()

    output_fname = sys.argv[3]
    ref_fname = sys.argv[2]    
    tandem_fname = sys.argv[1]
    
    chromosome_list = chr_tools.get_chr_list()
    
    rf_data = ref_and_tandem.data(ref_path=mydirs.joinpath(inputdir, 
                                                           ref_fname), 
                                  tandem_path=mydirs.joinpath(inputdir, 
                                                              tandem_fname),
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
            tandemfirstrow = calc_distance_from_exon(exon_coordinates, 
                                                     chromosome, 
                                                     tandemfirstrow)
            distances_dict[chromosome] = distances_list
            coordinates_dict[chromosome] = coordinates_list
            locations_dict[chromosome] = locations_list
            print('Done for %s' %chromosome)
            
        print distances_dict
        print coordinates_dict
        print locations_dict
    
    # Convert dictionary to list
    distances = []
    for dist_list in distances_dict.values():
        distances.append(dist_list)
    
    # Flatten list
    distances = list_tools.flatten_list(distances)
    # Plot results
    jplots.plot_binned_bar_graph(distances, 100, -100000, 100000, 
                                 'Distance from Exon', 'Frequency', 
                                 'Title', save_path_name,
                                 autoxlim=True)
    print('Plot saved in %s' %save_path_name)
            
        
    
            

        
        
        
    
    
    
    
