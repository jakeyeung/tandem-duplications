# Jake Yeung
# 28 May 2013
# motif_weights.R


motif_weights <- function(sequences, motif_length){
    # Calculate Weights of Each Sequence by Position
    base_frac_matrix <- matrix(0, nrow=4, ncol=motif_length)
    # rownames(base_frac_matrix) <- c('A', 'T', 'G', 'C')
    rownames(base_frac_matrix) <- c('A', 'C', 'G', 'T')    # Order is important!
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