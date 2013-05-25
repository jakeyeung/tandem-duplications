'''
Created on 2013-05-25

@author: jyeung
'''


import sys
from utilities import set_directories, ref_and_tandem


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


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print('Reference genome info and tandem dupes file must be given on the command line.')
        sys.exit()

    ref_fname = sys.argv[2]    
    tandem_fname = sys.argv[1]
    
    rf_data = ref_and_tandem.data(ref_path=mydirs.joinpath(inputdir, 
                                                           ref_fname), 
                                  tandem_path=mydirs.joinpath(inputdir, 
                                                              tandem_fname))
    with rf_data:
        # Define column names
        ref_headers = rf_data.refnext()
        tandem_headers = rf_data.tandemnext()
        
        # Define column indices
        ref_chr_index = ref_headers.index(ref_chr_colname)
        ref_exonstarts_index = ref_headers.index(ref_exonstarts_colname)
        ref_exonends_index = ref_headers.index(ref_exonends_colname)
        
        
        # Look at exon starts and ends.
        rowcount = 0
        while True:
            try:
                row = rf_data.refnext()
                
                # Convert exon starts into a list
                exonstarts_str = row[ref_exonstarts_index]
                # If last character is a comma, remove it. 
                if exonstarts_str[-1] == ',':
                    # print('Remove comma from %s' %exonstarts_str)
                    exonstarts_str = exonstarts_str[:-1]
                    # print('After removal, %s' %exonstarts_str)
                exonstarts_list = [int(i) for i in exonstarts_str.split(',')]
                print exonstarts_list
                
                # Convert exon ends into a list
                exonends_str = row[ref_exonends_index]
                # If last character is a comma, remove it.
                if exonends_str[-1] == ',':
                    # print('Remove comma from %s' %exonends_str)
                    exonends_str = exonends_str[:-1]
                    # print('After removal, %s' %exonends_str)
                exonends_list = [int(i) for i in exonends_str.split(',')]
                print exonends_list
                
                rowcount += 1
                break
            except StopIteration:
                print('Row count: %s. No more rows to iterate.' %rowcount)
                break
            

        
        
        
    
    
    
    
