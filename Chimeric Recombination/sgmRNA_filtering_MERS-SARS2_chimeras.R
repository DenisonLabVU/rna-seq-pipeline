#Authored by Jennifer Gribble. Last edited November 18, 2020.
#This script filters junctions forming putative chimeric sgmRNAs between MERS-CoV (JX869059.2) and SARS-CoV-2 (MT020881.1).
library(dplyr)
#Load data with 3 columns: "Start", "Stop", and "Depth". Hard code name changes for each sample below.
df_MERS_to_SARS2 <- read.table("/Users/gribblj/Dropbox/Perlman_MERS_SARS2/RNAseq/Junction_Files/0-01A_MERS_to_SARS2_junctions.txt", header = TRUE)
df_SARS2_to_MERS <- read.table("/Users/gribblj/Dropbox/Perlman_MERS_SARS2/RNAseq/Junction_Files/0-01A_SARS2_to_MERS_junctions.txt", header = TRUE)
#MERS-CoV TRS is located at nt 32-97
df_MERS_TRSL <- filter(df_MERS_to_SARS2, between(df_MERS_to_SARS2$Start, 32, 97))
df_MERS_to_SARS2_sgmRNA2 <- filter(df_MERS_TRSL, between(df_MERS_TRSL$Stop, 21525, 21592)) #SARS2 TRS2
df_MERS_to_SARS2_sgmRNA3 <- filter(df_MERS_TRSL, between(df_MERS_TRSL$Stop, 25354, 25421)) #SARS2 TRS3
df_MERS_to_SARS2_sgmRNA4 <- filter(df_MERS_TRSL, between(df_MERS_TRSL$Stop, 26206, 26273)) #SARS2 TRS4
df_MERS_to_SARS2_sgmRNA5 <- filter(df_MERS_TRSL, between(df_MERS_TRSL$Stop, 26442, 26509)) #SARS2 TRS5
df_MERS_to_SARS2_sgmRNA6 <- filter(df_MERS_TRSL, between(df_MERS_TRSL$Stop, 27010, 27077)) #SARS2 TRS6
df_MERS_to_SARS2_sgmRNA7 <- filter(df_MERS_TRSL, between(df_MERS_TRSL$Stop, 27357, 27424)) #SARS2 TRS7
df_MERS_to_SARS2_sgmRNA8 <- filter(df_MERS_TRSL, between(df_MERS_TRSL$Stop, 27857, 27924)) #SARS2 TRS8
df_MERS_to_SARS2_sgmRNA9 <- filter(df_MERS_TRSL, between(df_MERS_TRSL$Stop, 28229, 28296)) #SARS2 TRS9
#SARS-CoV-2 TRS-L and sgmRNA filtering. SARS-CoV-2 TRS-L between nt 40 - 105.
df_SARS2_TRSL <- filter(df_SARS2_to_MERS, between(df_SARS2_to_MERS$Start, 40, 105))
df_SARS2_to_MERS_sgmRNA2 <- filter(df_SARS2_TRSL, between(df_SARS2_TRSL$Stop, 21373, 21440)) #MERS TRS2
df_SARS2_to_MERS_sgmRNA3 <- filter(df_SARS2_TRSL, between(df_SARS2_TRSL$Stop, 25489, 25556)) #MERS TRS3
df_SARS2_to_MERS_sgmRNA4 <- filter(df_SARS2_TRSL, between(df_SARS2_TRSL$Stop, 25811, 25878)) #MERS TRS4
df_SARS2_to_MERS_sgmRNA5 <- filter(df_SARS2_TRSL, between(df_SARS2_TRSL$Stop, 26801, 26868)) #MERS TRS5
df_SARS2_to_MERS_sgmRNA6 <- filter(df_SARS2_TRSL, between(df_SARS2_TRSL$Stop, 27551, 27618)) #MERS TRS6
df_SARS2_to_MERS_sgmRNA7 <- filter(df_SARS2_TRSL, between(df_SARS2_TRSL$Stop, 27806, 27873)) #MERS TRS7
df_SARS2_to_MERS_sgmRNA8 <- filter(df_SARS2_TRSL, between(df_SARS2_TRSL$Stop, 28513, 28580)) #MERS TRS8
#Add column identifying sgmRNA species
df_MERS_to_SARS2_sgmRNA2$Type <- "sgmRNA2"
df_MERS_to_SARS2_sgmRNA3$Type <- "sgmRNA3"
df_MERS_to_SARS2_sgmRNA4$Type <- "sgmRNA4"
df_MERS_to_SARS2_sgmRNA5$Type <- "sgmRNA5"
df_MERS_to_SARS2_sgmRNA6$Type <- "sgmRNA6"
df_MERS_to_SARS2_sgmRNA7$Type <- "sgmRNA7"
df_MERS_to_SARS2_sgmRNA8$Type <- "sgmRNA8"
df_MERS_to_SARS2_sgmRNA9$Type <- "sgmRNA9"
#Add column identifying sgmRNA species
df_SARS2_to_MERS_sgmRNA2$Type <- "sgmRNA2"
df_SARS2_to_MERS_sgmRNA3$Type <- "sgmRNA3"
df_SARS2_to_MERS_sgmRNA4$Type <- "sgmRNA4"
df_SARS2_to_MERS_sgmRNA5$Type <- "sgmRNA5"
df_SARS2_to_MERS_sgmRNA6$Type <- "sgmRNA6"
df_SARS2_to_MERS_sgmRNA7$Type <- "sgmRNA7"
df_SARS2_to_MERS_sgmRNA8$Type <- "sgmRNA8"
#Make total sgmRNA table
df_MERS_to_SARS2_sgmRNA <- rbind(df_MERS_to_SARS2_sgmRNA2, df_MERS_to_SARS2_sgmRNA3, df_MERS_to_SARS2_sgmRNA4, df_MERS_to_SARS2_sgmRNA5, df_MERS_to_SARS2_sgmRNA6, df_MERS_to_SARS2_sgmRNA7, df_MERS_to_SARS2_sgmRNA8, df_MERS_to_SARS2_sgmRNA9)
df_MERS_to_SARS2_sgmRNA$Total <- sum(df_MERS_to_SARS2_sgmRNA$Depth)
df_SARS2_to_MERS_sgmRNA <- rbind(df_SARS2_to_MERS_sgmRNA2, df_SARS2_to_MERS_sgmRNA3, df_SARS2_to_MERS_sgmRNA4, df_SARS2_to_MERS_sgmRNA5, df_SARS2_to_MERS_sgmRNA6, df_SARS2_to_MERS_sgmRNA7, df_SARS2_to_MERS_sgmRNA8)
df_SARS2_to_MERS_sgmRNA$Total <- sum(df_SARS2_to_MERS_sgmRNA$Depth)
#Slice DVGs
df_MERS_to_SARS2_DVG <- anti_join(df_MERS_to_SARS2, df_MERS_to_SARS2_sgmRNA, by = c("Start", "Stop"))
df_MERS_to_SARS2_DVG$Total <- sum(df_MERS_to_SARS2_DVG$Depth)
df_SARS2_to_MERS_DVG <- anti_join(df_SARS2_to_MERS, df_SARS2_to_MERS_sgmRNA, by = c("Start", "Stop"))
df_SARS2_to_MERS_DVG$Total <- sum(df_SARS2_to_MERS_DVG$Depth)
#Write sliced dataframes. Hard code name changes below.
write.table(df_MERS_to_SARS2_sgmRNA, file = "/Users/gribblj/Dropbox/Perlman_MERS_SARS2/RNAseq/0-01A_total_MERS_to_SARS2_sgmRNAs.txt", sep = "\t", row.names = FALSE)
write.table(df_MERS_to_SARS2_DVG, file = "/Users/gribblj/Dropbox/Perlman_MERS_SARS2/RNAseq/0-01A_MERS_to_SARS2_DVGs.txt", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(df_SARS2_to_MERS_sgmRNA, file = "/Users/gribblj/Dropbox/Perlman_MERS_SARS2/RNAseq/0-01A_total_SARS2_to_MERS_sgmRNAs.txt", sep = "\t", row.names = FALSE)
write.table(df_SARS2_to_MERS_DVG, file = "/Users/gribblj/Dropbox/Perlman_MERS_SARS2/RNAseq/0-01A_SARS2_to_MERS_DVGs.txt", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)