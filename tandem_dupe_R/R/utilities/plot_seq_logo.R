# Jake Yeung
# 28 May 2013
# plot_seq_logo.R

# After using motif_weights.R function, we take our output matrix
# containing a weighted fraction of A, C, G, T (respectively) at 
# each position and put it into plot_seq_logo.R function to 
# plot a seqLogo.


# loadLibs ----------------------------------------------------------------

library(seqLogo)

plot_seq_logo <- function(base_frac_matrix, title, 
                          x=0.5, y=1, hjust=0.5, vjust=1){
    p <- makePWM(base_frac_matrix)
    seqLogo(p, ic.scale=FALSE)
    grid.text(title, x=0.5, y=1, hjust=0.5, vjust=1)
    return(p)
}