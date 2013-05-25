'''
Created on 2013-05-25

@author: jyeung
'''


import sys
from utilities import set_directories, ref_and_tandem, chr_tools, list_tools


# Set directories
mydirs = set_directories.my_directories()
inputdir = mydirs.joinpath(mydirs.root, 'inputs')
outputdir = mydirs.joinpath(mydirs.root, 'outputs')

# Define constants
ref_chr_colname = 'chrom'
ref_exonstarts_colname = 'exonStarts'
ref_exonends_colname = 'exonEnds'
tandem_chr_colname = 'chromosome_1'    # Equal to 'chromosome_2'
tandem_start_colname = 'start_1'
tandem_end_colname = 'end_2'


def find_exon_coordinates(rf_data, chromosome, firstrow):
    '''
    Find exon coordinates, start and stop locations in a tuple for 
    a specific chromosome.  
    '''
    # Initialze some dicts
    exon_dict = {}
    
    # Process first row
    cur_chr = firstrow[rf_data.refheaders.index(ref_chr_colname)]
    if cur_chr == chromosome:
        for colname in [ref_exonstarts_colname, ref_exonends_colname]:
            # Make empty list for corresponding key
            exon_dict[colname] = []
            exon_str = firstrow[rf_data.refheaders.index(colname)]
            # If last character is a comma, remove it. 
            if exon_str[-1] == ',':
                exon_str = exon_str[:-1]
            # Add to dict a list of integers representing coordinates
            exon_dict[colname].append([int(i) for i in exon_str.split(',')])
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
                for colname in [ref_exonstarts_colname, ref_exonends_colname]:
                    exon_str = row[rf_data.refheaders.index(colname)]
                    # If last character is a comma, remove it. 
                    if exon_str[-1] == ',':
                        exon_str = exon_str[:-1]
                    # Add to dict a list of integers representing coordinates
                    exon_dict[colname].append([int(i) for i in exon_str.split(',')])
            else:
                print('Current row is now in %s, breaking...' %cur_chr)
                lastrow = row
                break
        except StopIteration:
            print('Row count: %s. No more rows to iterate.' %rf_data.refrowcount)
            break
    
    # Flatten list of list before returning it to user.
    for colname, nested_list in exon_dict.iteritems():
        exon_dict[colname] = list_tools.flatten_list(nested_list)
    return exon_dict, lastrow
    
    
    '''
    while True:
        try:
            row = rf_data.refnext()
            curr_chr = rf_data.refheaders.index(ref_chr_colname)
            if curr_chr != prev_chr:
                pass
            exon_dict = {}
            # Convert a string separated by commas 
            # into separate lists, store into dict
            for colname in [ref_exonstarts_colname, ref_exonends_colname]:
                exon_str = row[rf_data.refheaders.index(colname)]
                # If last character is a comma, remove it. 
                if exon_str[-1] == ',':
                    exon_str = exon_str[:-1]
                # Add to dict a list of integers representing coordinates
                exon_dict[colname] = [int(i) for i in exon_str.split(',')]
            
            prev_chr = curr_chr
            break
        except StopIteration:
            print('Row count: %s. No more rows to iterate.' %rowcount)
            break
    print exon_dict
    '''
    
    
if __name__ == '__main__':
    if len(sys.argv) < 3:
        print('Reference genome info and tandem dupes file must be given on the command line.')
        sys.exit()

    ref_fname = sys.argv[2]    
    tandem_fname = sys.argv[1]
    
    chromosome_list = chr_tools.get_chr_list()
    
    rf_data = ref_and_tandem.data(ref_path=mydirs.joinpath(inputdir, 
                                                           ref_fname), 
                                  tandem_path=mydirs.joinpath(inputdir, 
                                                              tandem_fname))
    
    with rf_data:
        # Supply its first row
        print rf_data.refheaders
        print rf_data.tandemheaders
        # print rf_data.refnext()
        # print rf_data.tandemnext()
        firstrow = rf_data.refnext()
        for chromosome in chromosome_list:
            # Input first row of chr1, get firstrow of chr2 back and exon_dict
            # containing exonStart and End positions of chr1. Repeat for all
            # chromosomes. 
            exon_dict, firstrow = find_exon_coordinates(rf_data, chromosome, firstrow)
            
            break
            
            
        
    
            

        
        
        
    
    
    
    
