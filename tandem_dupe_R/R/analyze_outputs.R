# Jake Yeung
# May 26 2013
# analyze_outputs.R
# use ggplot2 to look at some python results.


# ImportingLibraries ------------------------------------------------------

library(ggplot2)



# SetDirectories ----------------------------------------------------------

cur_dir <- getwd()
root_dir <- dirname(dirname(dirname(dirname(cur_dir))))
input_dir <- file.path(root_dir, 'inputs')
output_dir <- file.path(root_dir, 'outputs')



# SetConstants ------------------------------------------------------------

filename <- 'tandem_dupes_data_appended.txt'


# ReadData ----------------------------------------------------------------

tandem_dat <- read.table(file.path(input_dir, filename), header=TRUE, sep='\t')



# PlotDat -----------------------------------------------------------------

qplot(distance_to_tss, data=tandem_dat, binwidth=50000) + facet_wrap(~chromosome_1)
qplot(seq_length, data=tandem_dat, binwidth=20, fill=exon_or_nonexon)
qplot(distance_to_exon, data=tandem_dat, binwidth=50000) + facet_wrap(~chromosome_1)
qplot(gc_content, data=tandem_dat, binwidth=0.02) + facet_wrap(~chromosome_1)



