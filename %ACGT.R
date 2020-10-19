dat <- read.table("sample_DVG_start_sequences.txt", header = TRUE)

rows = nrow(dat)

dat <- dat[4]

new <- matrix(nrow = 61, ncol = rows)
for(i in 1:61){
  for(j in 1:rows){
    new [i,j] <- substring(dat[j,], i, i)
  }
}

countTable <- matrix(nrow = 61, ncol = 4)
for(i in 1:61){
  columnSeq <- DNAStringSet(paste0(new[i,], collapse = ""))
  columnCounts <- letterFrequency(columnSeq, letters = "ACGT", OR = 0)
  countTable[i,] <- columnCounts
}

colnames(countTable) <- c("A", "C", "G", "T")

freqTable <- countTable/rows
df<- round(t(freqTable), digit = 4)
df <- df * 100
write.table(df, file = "sample_DVG_start_%ACGT.txt", sep = "\t")
