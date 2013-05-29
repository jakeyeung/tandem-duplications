# Jake Yeung
# 28 May 2013
# analyze_startend_motifs.R



# ImportLibraries ---------------------------------------------------------


# SetDirectories ----------------------------------------------------------

cur_dir <- getwd()
root_dir <- dirname(dirname(dirname(dirname(cur_dir))))
source_dir <- file.path(cur_dir, 'utilities')
input_dir <- file.path(root_dir, 'inputs')
output_dir <- file.path(root_dir, 'outputs')


# SetConstants ------------------------------------------------------------

seq_file_path <- file.path(output_dir, 'parsed_sequence.txt')
motif_length = 30


# Functions ---------------------------------------------------------------

source(file.path(source_dir, 'motif_weights.R'))
source(file.path(source_dir, 'plot_seq_logo.R'))

# LoadData ----------------------------------------------------------------

seq_dat <- read.table(seq_file_path, header=TRUE, sep='\t')


# PlotSeqLogo -------------------------------------------------------------

# Get sequences
seq_start <- seq_dat[, 'seq_start']
seq_end <- seq_dat[, 'seq_end']

base_frac_matrix_start <- motif_weights(seq_start, motif_length)
base_frac_matrix_end <- motif_weights(seq_end, motif_length)

pStart <- plot_seq_logo(base_frac_matrix_start, 'startSeq')
pEnd <- plot_seq_logo(base_frac_matrix_end, 'endSeq')
