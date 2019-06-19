### This script performs block jackknifing to obtain significance of genome-wide estimates of Patterson's D and is adapted from the tutorial here http://evomics.org/learning/population-and-speciation-genomics/2018-population-and-speciation-genomics/abba-baba-statistics/

### Read in the allele frequency data calculated using this script 'https://github.com/simonhmartin/genomics_general/freq.py'. The 'ancestral' allele is obtained based on sites that are fixed within M. clevelandii.

freq_table = read.table('all_9_taxa_G1_vars_geno.freq.derived.LG_only.txt', header =T, as.is = T)
freq_table <- na.omit(freq_table)

### Read the coordinates of genomic windows for jackknifing
genome_wins <- read.table('1Mb_wins.txt')
genome_wins <- genome_wins[order(genome_wins$V1, genome_wins$V2),]


### Calculate the observed D statistic

abba = function(p1, p2, p3) (1 - p1) * p2 * p3

baba = function(p1, p2, p3) p1 * (1 - p2) * p3

D.stat = function(dataframe) (sum(dataframe$ABBA) - sum(dataframe$BABA)) / (sum(dataframe$ABBA) + sum(dataframe$BABA))

P1 = "LON"
P2 = "Y"
P3 = "ARI"

ABBA = abba(freq_table[,P1], freq_table[,P2], freq_table[,P3])
BABA = baba(freq_table[,P1], freq_table[,P2], freq_table[,P3])

ABBA_BABA_df = as.data.frame(cbind(ABBA,BABA))

D = D.stat(ABBA_BABA_df)




D_pseudo <- c()

### Re-calculate D while iteratively removing one 1Mb window at a time

for(i in 1: length(genome_wins[,1])) {
  curr_lg <- genome_wins[[1]][i]
  start <- genome_wins[[2]][i]
  stop <- genome_wins[[3]][i]
  
  window_frq <- freq_table[freq_table$scaffold == curr_lg & freq_table$position >= start & freq_table$position <= stop,]
  
  first_pos <- which(freq_table$scaffold == curr_lg & freq_table$position==window_frq[[2]][1])
  last_pos <- which(freq_table$scaffold == curr_lg & freq_table$position==window_frq[[2]][length(window_frq[,1])])
  
  freq_table_minus_win <- freq_table[-c(first_pos:last_pos),]
  
  ABBA = abba(freq_table_minus_win[,P1], freq_table_minus_win[,P2], freq_table_minus_win[,P3])
  BABA = baba(freq_table_minus_win[,P1], freq_table_minus_win[,P2], freq_table_minus_win[,P3])
  
  ABBA_BABA_df = as.data.frame(cbind(ABBA,BABA))
  
  D_jack <- D.stat(ABBA_BABA_df)
  
  D_diff <- D - D_jack
  
  D_pseudo <- append(D_pseudo, D_diff)
}

### Calculate p value
D_sd <- sd(D_pseudo)
D_err <- D_sd/sqrt(length(genome_wins[,1]))
D_Z <- D/D_err
D_p <- 2*pnorm(-abs(D_Z))






