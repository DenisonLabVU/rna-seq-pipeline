library(dplyr)
library(Biostrings)
##Load in data, save quantification of rows as variable, and slice start and stop sequences
dat <- read.table("/Users/jennifergribble/Dropbox/P250_recombination/Passage_Populations/NT_frequency/WTP3A_DVGs_sequences.txt", header = FALSE)
rows = nrow(dat)
dat_start <- select(dat, V9)
dat_stop <- select(dat, V10)
##generate matrix of sequences
new <- matrix(nrow = 41, ncol = rows)
for(i in 1:41){
  for(j in 1:rows){
    new [i,j] <- substring(dat_start[j,], i, i)
  }
}
##Count matrix
countTable <- matrix(nrow = 41, ncol = 4)
for(i in 1:41){
  columnSeq <- DNAStringSet(paste0(new[i,], collapse = ""))
  columnCounts <- letterFrequency(columnSeq, letters = "ACGT", OR = 0)
  countTable[i,] <- columnCounts
}
##Rename columns, calculate frequency, and save
colnames(countTable) <- c("A", "C", "G", "T")
freqTable <- countTable/rows
df<- round(t(freqTable), digit = 4)
df <- df * 100
write.table(df, file = "/Users/jennifergribble/Dropbox/P250_recombination/Passage_Populations/NT_frequency/WTP3A_start_%ACGU.txt", sep = "\t", quote = FALSE)
