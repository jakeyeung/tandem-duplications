# Jake Yeung
# May 26 2013
# analyze_outputs.R
# use ggplot2 to look at some python results.


# ImportingLibraries ------------------------------------------------------

library(ggplot2)
# library(seqLogo)
library(scales)
options(scipen=1000)    # Replaces scientific notation to numbers in plots



# SetDirectories ----------------------------------------------------------

cur_dir <- getwd()
root_dir <- dirname(dirname(dirname(dirname(cur_dir))))
input_dir <- file.path(root_dir, 'inputs')
output_dir <- file.path(root_dir, 'outputs')



# SetConstants ------------------------------------------------------------

filename <- 'tandem_dupes_data_appended.txt'
filename_nulldist_tss <- 'tss_distance_null_distribution.txt'
filename_nulldist_exon <- 'exon_distance_null.txt'
plotfolder <- 'plots'
histoplots_fname <- 'ggplot2_tandem_duplication_summary.pdf'
seqplots_fname <- 'seqlogo_beginning_and_end2.pdf' 
tandem_and_null_fname <- 'ggplot2_tandem_and_null_bychr_withexon2.pdf'
motif_length <- 30


# DefineFunctions ---------------------------------------------------------


# ReadData ----------------------------------------------------------------

tandem_dat <- read.table(file.path(input_dir, filename), header=TRUE, sep='\t')
null_tts_dist_dat <- read.table(file.path(output_dir, filename_nulldist_tss), header=TRUE, sep='\t')
null_exon_dist_dat <- read.table(file.path(output_dir, filename_nulldist_exon), header=TRUE, sep='\t')


# AppendtoOneDF -----------------------------------------------------------

# There are a few more rows in exon than tss, remove some to make them the same.
row_diff <- nrow(null_exon_dist_dat) - nrow(null_tts_dist_dat)
if(row_diff > 0){
    random_rows <- sample(1:nrow(null_exon_dist_dat), row_diff)
    null_exon_dist_dat <- null_exon_dist_dat[-random_rows, ]
} else if(row_diff < 0){
    random_rows <- sample(1:nrow(null_tts_dist_dat), row_diff)
    null_tts_dist_dat <- null_tts_dist_dat[-random_rows, ]
}
# Do nothing if row_diff == 0

# Cbind the rows together **FIX ME** CHROMOSOMES MAY NOT MATCH
null_dat <- data.frame(null_exon_dist_dat[, c('chromosome_1', 'distance_to_exon')])
null_dat <- cbind(null_dat, null_tts_dist_dat$distance_to_tss)


# cbind whether it is tandem or it is null to each respective df
tandem_dat.df <- data.frame(tandem_dat, 'type'='tandem')
null_dat.df <- data.frame(null_dat, 'type'='null')


extract_cnames <- c('chromosome_1', 'distance_to_exon', 'distance_to_tss', 'type')
colnames(null_dat.df) <- extract_cnames
# Rbind the two DFs
full_dat <- rbind(tandem_dat.df[, extract_cnames], null_dat.df)


# KolmogorovSmirnov -------------------------------------------------------

# Transcription start site test
ks_results_list_tss <- vector('list', length(levels(tandem_dat$chromosome_1)))
i <- 0
for(chromosome in levels(tandem_dat$chromosome_1)){
    i <- i + 1
    dat_sub_tandem <- tandem_dat[which(tandem_dat$chromosome_1 == chromosome), ]
    dat_sub_nulltss <- null_tts_dist_dat[which(null_tts_dist_dat$chromosome == chromosome), ]
    ks_results_list_tss[i] <- ks.test(dat_sub_tandem$distance_to_tss, dat_sub_nulltss$distance_to_tss)$p.value
}

# Exon location test
ks_results_list_exon <- vector('list', length(levels(tandem_dat$chromosome_1)))
i <- 0
for(chromosome in levels(tandem_dat$chromosome_1)){
    i <- i + 1
    dat_sub_tandem <- tandem_dat[which(tandem_dat$chromosome_1 == chromosome), ]
    dat_sub_nullexon <- null_exon_dist_dat[which(null_exon_dist_dat$chromosome_1 == chromosome), ]
    ks_results_list_exon[i] <- ks.test(dat_sub_tandem$distance_to_exon, dat_sub_nullexon$distance_to_exon)$p.value
}


# PlotDat -----------------------------------------------------------------

pdf(file.path(output_dir, plotfolder, histoplots_fname), width=20, height=15)
# Distance to tss
qplot(distance_to_tss, data=tandem_dat, binwidth=500, 
      fill=exon_or_nonexon) + xlim(-10000, 10000)
qplot(distance_to_tss, data=tandem_dat, binwidth=1000, 
      fill=exon_or_nonexon) + facet_wrap(~chromosome_1) + xlim(-10000, 10000)
# Sequence length
qplot(seq_length, data=tandem_dat, binwidth=10)
qplot(seq_length, data=tandem_dat, binwidth=20) + facet_wrap(~chromosome_1)
# Distance to exon
qplot(distance_to_exon, data=tandem_dat, binwidth=500, 
      fill=exon_or_nonexon) + xlim(-10000, 10000)
qplot(distance_to_exon, data=tandem_dat, binwidth=1000, 
      fill=exon_or_nonexon) + facet_wrap(~chromosome_1) + xlim(-10000, 10000)
# GC content
qplot(gc_content, data=tandem_dat, binwidth=0.02, 
      fill=exon_or_nonexon) + xlim(0, 1)
qplot(gc_content, data=tandem_dat, binwidth=0.04, 
      fill=exon_or_nonexon) + facet_wrap(~chromosome_1) + xlim(0, 1)
dev.off()


# CompareWithNull ---------------------------------------------------------

pdf(file.path(output_dir, plotfolder, tandem_and_null_fname), width=20, height=15)
# Distance to tss
# Genome wide view
ggplot(full_dat, aes(x=distance_to_tss, fill=type)) + 
    geom_density(alpha=0.5) + xlim(-1000000, 1000000) + 
    geom_rug(data=full_dat[full_dat$type=='tandem', ], aes(x=distance_to_tss)) + xlab('Distance to TSS') + 
    ggtitle('Density plot: distance to TSS. Rug indicates tandem duplications')
# By chromosome
ggplot(full_dat, aes(x=distance_to_tss, fill=type)) + 
    geom_density(alpha=0.5) + xlim(-100000, 100000) + ylim(0, 4e-5) + 
    geom_rug(data=full_dat[full_dat$type=='tandem', ], aes(x=distance_to_tss)) + 
    facet_wrap(~chromosome_1) + xlab('Distance to TSS') + 
    ggtitle('Density plot: distance to TSS. Rug indicates tandem duplications')
# Distance to some exon
# Genome wide view
ggplot(full_dat, aes(x=distance_to_exon, fill=type)) + 
    geom_density(alpha=0.5) + xlim(-100000, 100000) + 
    geom_rug(data=full_dat[full_dat$type=='tandem', ], aes(x=distance_to_tss)) + xlab('Distance to Exon') + 
    ggtitle('Density plot: distance to Exon. Rug indicates tandem duplications')
# By chromosome
ggplot(full_dat, aes(x=distance_to_exon, fill=type)) + 
    geom_density(alpha=0.5) + xlim(-10000, 10000) + ylim(0, 4e-4) + 
    geom_rug(data=full_dat[full_dat$type=='tandem', ], aes(x=distance_to_tss)) + 
    facet_wrap(~chromosome_1) + xlab('Distance to an Exon') + 
    ggtitle('Density plot: distance to an Exon. Rug indicates tandem duplications')
dev.off()
