library(dplyr)
library(Biostrings)
##Load in data, save quantification of rows as variable, and slice start and stop sequences
dat <- read.table("/Users/jennifergribble/Dropbox/P250_recombination/Passage_Populations/NT_frequency/WTP3B_DVGs_sequences.txt", header = FALSE)
n = nrow(dat)
dat_start <- select(dat, V9)
dat_stop <- select(dat, V10)
##generate matrix of sequences
new <- matrix(nrow = 41, ncol = n)
for(i in 1:41){
  for(j in 1:n){
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
colnames(countTable) <- c("A", "C", "G", "U")
freqTable <- countTable/n
df<- round(t(freqTable), digit = 4)
df <- df * 100
df <- as.data.frame(t(df))
##Uncomment below if you are doing start or stop sequences. +1 always indicates the junction participating nucleotide. Position numbers indicate upstream sequence, negative numbers indicate downstream sequence.
#vec_start <- c("+21", "+20", "+19", "+18", "+17", "+16", "+15", "+14", "+13", "+12", "+11", "+10", "+9", "+8", "+7", "+6", "+5", "+4", "+3", "+2", "+1", "-1", "-2", "-3", "-4", "-5", "-6", "-7", "-8", "-9", "-10", "-11", "-12", "-13", "-14", "-15", "-16", "-17", "-18", "-19", "-20")
#df$Position <- vec_start
#vec_stop <- c("-20", "-19", "-18", "-17", "-16", "-15", "-14", "-13", "-12", "-11", "-10", "-9", "-8", "-7", "-6", "-5", "-4", "-3", "-2", "-1", "+1", "+2", "+3", "+4", "+5", "+6", "+7", "+8", "+9", "+10", "+11", "+12", "+13", "+14", "+15", "+16", "+17", "+18", "+19", "+20")
#df$Position <- vec_stop
df <- df[c(5,1,2,3,4)]
write.table(df, file = "/Users/jennifergribble/Dropbox/P250_recombination/Passage_Populations/NT_frequency/WTP3B_start_%ACGU.txt", sep = "\t", quote = FALSE, row.names = FALSE)
