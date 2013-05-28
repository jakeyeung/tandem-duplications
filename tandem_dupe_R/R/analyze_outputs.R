# Jake Yeung
# May 26 2013
# analyze_outputs.R
# use ggplot2 to look at some python results.


# ImportingLibraries ------------------------------------------------------

library(ggplot2)
library(seqLogo)



# SetDirectories ----------------------------------------------------------

cur_dir <- getwd()
root_dir <- dirname(dirname(dirname(dirname(cur_dir))))
input_dir <- file.path(root_dir, 'inputs')
output_dir <- file.path(root_dir, 'outputs')



# SetConstants ------------------------------------------------------------

filename <- 'tandem_dupes_data_appended.txt'
filename_nulldist <- 'tss_distance_null_distribution.txt'
plotfolder <- 'plots'
histoplots_fname <- 'ggplot2_tandem_duplication_summary.pdf'
seqplots_fname <- 'seqlogo_beginning_and_end2.pdf' 
tandem_and_null_fname <- 'ggplot2_tandem_and_null_bychr2.pdf'
motif_length <- 30


# DefineFunctions ---------------------------------------------------------

motif_weights <- function(sequences, motif_length){
    # Calculate Weights of Each Sequence by Position
    base_frac_matrix <- matrix(0, nrow=4, ncol=motif_length)
    # rownames(base_frac_matrix) <- c('A', 'T', 'G', 'C')
    rownames(base_frac_matrix) <- c('A', 'C', 'G', 'T')
    for(i in 1:motif_length){
        Acount <- 0
        Tcount <- 0
        Gcount <- 0
        Ccount <- 0
        for(j in 1:length(sequences)){
            base <- substring(sequences[j], i, i)
            if (base == 'A'){
                Acount <- Acount + 1
            } else if (base == 'T'){
                Tcount <- Tcount + 1
            } else if (base == 'G'){
                Gcount <- Gcount + 1
            } else if (base == 'C'){
                Ccount <- Ccount + 1
            } else{
                stop(c('Unknown base of: ', base))
            }
        }
        total_count <- Acount + Tcount + Gcount + Ccount
        base_frac_matrix['A', i] <- Acount / total_count
        base_frac_matrix['T', i] <- Tcount / total_count
        base_frac_matrix['G', i] <- Gcount / total_count
        base_frac_matrix['C', i] <- Ccount / total_count
    }
    return(base_frac_matrix)
}
# ReadData ----------------------------------------------------------------

tandem_dat <- read.table(file.path(input_dir, filename), header=TRUE, sep='\t')
null_dat <- read.table(file.path(output_dir, filename_nulldist), header=TRUE, sep='\t')


# AppendtoOneDF -----------------------------------------------------------

# cbind whether it is tandem or it is null to each respective df
tandem_dat.df <- data.frame(tandem_dat, 'type'='tandem')
null_dat.df <- data.frame(null_dat, 'type'='null')


extract_cnames <- c('chromosome_1', 'tandem_dupe_position', 'distance_to_tss', 'tss_location', 'type')
colnames(null_dat.df) <- extract_cnames
# Rbind the two DFs
full_dat <- rbind(tandem_dat.df[, extract_cnames], null_dat.df)


# KolmogorovSmirnov -------------------------------------------------------

ks_results_list <- vector('list', length(levels(tandem_dat$chromosome_1)))
i <- 0
for(chromosome in levels(tandem_dat$chromosome_1)){
    i <- i + 1
    dat_sub_tandem <- tandem_dat[which(tandem_dat$chromosome_1 == chromosome), ]
    dat_sub_null <- null_dat[which(null_dat$chromosome_1 == chromosome), ]
    ks_results_list[i] <- ks.test(dat_sub_tandem$distance_to_tss, null_dat$distance_to_tss[rand_rows])$p.value
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
qplot(distance_to_tss, data=full_dat, fill=type, geom='density', 
      alpha=I(.5)) + xlim(-100000, 100000) + facet_wrap(~chromosome_1) + geom_rug()
dev.off()

# SeqLogo -----------------------------------------------------------------

# Look at first positions of the sequence.
sequences <- list(tandem_dat$sequence)
startSeq <- sapply(sequences, function(x){
    return(substring(as.character(x), 1, motif_length))
})
# Look at last positions of the sequence.
endSeq <- sapply(sequences, function(x){
    x <- as.character(x)
    return(substring(x, 
                     nchar(x)-(motif_length-1), 
                     nchar(x)))
})

base_weight_frac_start <- motif_weights(startSeq, motif_length)
base_weight_frac_start
base_weight_frac_end <- motif_weights(endSeq, motif_length)
base_weight_frac_end

# Use SeqLogo library to make plot
pdf(file.path(output_dir, plotfolder, seqplots_fname))
pStart <- makePWM(base_weight_frac_start)
seqLogo(pStart, ic.scale=FALSE)
grid.text(sprintf('Sequence from Beginning to Base 30'), x=0.5, y=1, hjust=0.5, vjust=1)
# Plot
pEnd <- makePWM(base_weight_frac_end)
seqLogo(pEnd, ic.scale=FALSE)
grid.text(sprintf('Sequence from End-30 to End'), x=0.5, y=1, hjust=0.5, vjust=1)
dev.off()