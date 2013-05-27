# Jake Yeung 
# May 26 2013
# bind_results.R
# Bind some outputs from python into a single dataframe.

# SetDirectories ----------------------------------------------------------

cur_dir <- getwd()
root_dir <- dirname(dirname(dirname(dirname(cur_dir))))
input_dir <- file.path(root_dir, 'inputs')
output_dir <- file.path(root_dir, 'outputs')


# SetFileNames ------------------------------------------------------------

seq_fname <- 'sequence_summary.txt'
exon_dist_fname <- 'exon_distances_from_start.txt'
tss_dist_fname <- 'distance_to_tss.txt'
output_filename <- 'tandem_dupes_data_appended.txt'


# SetConstants ------------------------------------------------------------
exon_dist_cnames <- c('tandem_dupe_position', 
                      'distance_to_exon', 
                      'closest_exon_start', 
                      'closest_exon_end', 
                      'exon_or_nonexon')

tss_dist_cnames <- c('distance_to_tss', 'tss_location')



# DefineFunctions ---------------------------------------------------------

jReadTable <- function(filename){
    df <- read.table(file.path(output_dir, filename), header=TRUE, sep='\t')
    return(df)
}


# LoadFiles ---------------------------------------------------------------

seq_dat <- jReadTable(seq_fname)
exon_dist_dat <- jReadTable(exon_dist_fname)
tss_dist_dat <- jReadTable(tss_dist_fname)


# BindColumns -------------------------------------------------------------

binded_df <- cbind(seq_dat, 
                   exon_dist_dat[, exon_dist_cnames], 
                   tss_dist_dat[, tss_dist_cnames])


# WriteDFtoInputFolder ----------------------------------------------------

write.table(binded_df, file=file.path(input_dir, output_filename), 
            quote=FALSE, sep='\t', row.names=FALSE)






