'''
Created on 2013-05-19

@author: jyeung
'''


import sys
import csv


def get_tandem_tss_distribution(file_path, file_start_colname, 
                                file_chr_colname, 
                                genome_tss, chromosome, chromosome_list):
    '''
    An oldish piece of code replaced by 
    tandem_dupes_all and tandem_dupes_all_pairs.
    '''
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

    
def tandem_dupes_all(file_path, file_start_colname,
                        file_chr_colname, 
                        all_genome_tss, chromosome_list):
    
    with open(file_path, 'rb') as read_file:
        file_reader = csv.reader(read_file, delimiter='\t')
        file_colnames = file_reader.next()
        
        try: 
            file_start_index = file_colnames.index(file_start_colname)
        except ValueError:
            print('ValueError: couldnt find match %s to ' \
             'start column name. Possible values are: ' %file_start_colname)
            print file_colnames
            sys.exit()
        
        try: 
            file_chr_index = file_colnames.index(file_chr_colname)
        except ValueError:
            print('ValueError: couldnt find match %s to ' \
             'chromosome column name. Possible values are: ' %file_chr_colname)
            print file_colnames
            sys.exit()
            
        # Initialize with empty list
        tandem_dupe_pos = {}
        for c in chromosome_list:
            tandem_dupe_pos[c] = []
        
        for file_row in file_reader:
            current_chr = file_row[file_chr_index]
            current_pos = file_row[file_start_index]
            
            if current_chr in tandem_dupe_pos:
                tandem_dupe_pos[current_chr].append(int(current_pos))
            else:
                print('Warning: cannot find %s in chromosome list' %current_chr)
                sys.exit('Exiting...')
        
        count = 0
        for l in tandem_dupe_pos.values():
            count += len(l)
        print('%s tandem dupes found' %count)
        return tandem_dupe_pos
    
    
def tandem_dupes_all_pairs(file_path, file_start_colname,
                           file_end_colname,
                           file_chr_colname, chromosome_list):
    
    with open(file_path, 'rb') as read_file:
        file_reader = csv.reader(read_file, delimiter='\t')
        file_colnames = file_reader.next()
        
        try: 
            file_start_index = file_colnames.index(file_start_colname)
        except ValueError:
            print('ValueError: couldnt find match %s to ' \
             'start column name. Possible values are: ' %file_start_colname)
            print file_colnames
            sys.exit()
        try: 
            file_end_index = file_colnames.index(file_end_colname)
        except ValueError:
            print('ValueError: couldnt find match %s to ' \
             'end column name. Possible values are: ' %file_end_colname)
            print file_colnames
            sys.exit()
        try: 
            file_chr_index = file_colnames.index(file_chr_colname)
        except ValueError:
            print('ValueError: couldnt find match %s to ' \
             'chromosome column name. Possible values are: ' %file_chr_colname)
            print file_colnames
            sys.exit()
            
        # Initialize with empty list
        tandem_dupe_pos_pairs = {}
        for c in chromosome_list:
            tandem_dupe_pos_pairs[c] = []
            
        for file_row in file_reader:
            current_chr = file_row[file_chr_index]
            current_pos = file_row[file_start_index]
            current_pos_end = file_row[file_end_index]
        
            if current_chr in tandem_dupe_pos_pairs:
                tandem_dupe_pos_pairs[current_chr].append((int(current_pos), 
                                                           int(current_pos_end)))
            else:
                print('Warning: cannot find %s in chromosome list' %current_chr)
                sys.exit('Exiting...')
                
        count = 0
        for l in tandem_dupe_pos_pairs.values():
            count += len(l)
        print('%s tandem dupes pairs found' %count)
        return tandem_dupe_pos_pairs
    
    
def calculate_tandem_gc_content(file_path, file_seq_colname, file_chr_colname, 
                                output_path, *output_additional_colnames):
    '''
    Reads tandem duplication file, calculates gc content and sequence length.
    Outputs gc content and length of sequence to writefile. 
    '''
    with open(file_path, 'rb') as read_file:
        with open(output_path, 'wb') as write_file:
            # Initialize and get colnames
            file_reader = csv.reader(read_file, delimiter='\t')
            file_colnames = file_reader.next()
            output_writer = csv.writer(write_file, delimiter='\t')
            # Create output colnames by adding additional colnames
            output_colnames = file_colnames
            for cname in output_additional_colnames:
                # Should have three additional colnames:
                # sequence, gc_content and length_of_sequence
                output_colnames.append(cname)
            # Write colnames to output file.
            output_writer.writerow(output_colnames)
            
            # Get index numbers for reader           
            try: 
                file_gc_index = file_colnames.index(file_seq_colname)
            except ValueError:
                print('ValueError: couldnt find match %s to ' \
                 'sequence column name.' %file_seq_colname)
                sys.exit()
            try: 
                file_chr_index = file_colnames.index(file_chr_colname)
            except ValueError:
                print('ValueError: couldnt find match %s to ' \
                 'chromosome column name.' %file_chr_colname)
                sys.exit()
            
            gc_content = []    # List of tuples, content and chromosome.
            for file_row in file_reader:
                sequence = file_row[file_gc_index]
                chromosome = file_row[file_chr_index]
                gc_count = 0
                for base in sequence:
                    if base in ['G', 'C']:
                        gc_count += 1
                gc_frac = float(gc_count)/len(sequence)
                gc_content.append((gc_frac, chromosome, len(sequence)))
                # Write to file
                file_row.append(gc_frac)
                file_row.append(len(sequence))
                if len(file_row) == len(output_colnames):
                    output_writer.writerow(file_row)
                else:
                    print('File row length %s not equal to length of colnames %s' %(len(file_row), len(output_colnames)))
                    sys.exit()
    return gc_content
                    
            
            
        
    