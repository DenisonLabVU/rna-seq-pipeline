#Authored by Jennifer Gribble. November 18, 2020.
#Generates summary data for MERS-CoV and SARS-CoV-2 co-infection data. Must run separately  for each sample. Aggregate data can be output at final sample.
library(dplyr)
#Hardcode different file names.
df_MERS_to_SARS2 <- read.table("/Users/gribblj/Dropbox/Perlman_MERS_SARS2/RNAseq/Junction_Files/0-01C_MERS_to_SARS2_junctions.txt", header = TRUE)
df_SARS2_to_MERS <- read.table("/Users/gribblj/Dropbox/Perlman_MERS_SARS2/RNAseq/Junction_Files/0-01C_SARS2_to_MERS_junctions.txt", header = TRUE)
df_total <- read.table("/Users/gribblj/Dropbox/Perlman_MERS_SARS2/RNAseq/BED_Files/0-01C_virema__Virus_Recombination_Results.bed", skip = 1, header = FALSE)
df_total <- df_total[c(1,2,3,5)]
df_total <- df_total %>% rename("Genome" = V1, "Start" = V2, "Stop" = V3, "Depth" = V5)
df_MERS <- df_total[which(df_total$Genome=="JX869059.2"), ]
df_SARS2 <- df_total[which(df_total$Genome=="MT020881.1"), ]
#make new dataframe and add calculated data
df_summary <- data.frame(row.names=1)
df_summary$Sample <- "0.01C" #Change sample name for each sample
df_summary$Total_Junctions <- nrow(df_total)
df_summary$Total_Depth <- sum(df_total$Depth)
df_summary$MERS_Junctions <- nrow(df_MERS)
df_summary$MERS_Depth <- sum(df_MERS$Depth)
df_summary$SARS2_Junctions <- nrow(df_SARS2)
df_summary$SARS2_Depth <- sum(df_SARS2$Depth)
df_summary$MERS_to_SARS2_Junctions <- nrow(df_MERS_to_SARS2)
df_summary$MERS_to_SARS2_Depth <- sum(df_MERS_to_SARS2$Depth)
df_summary$SARS2_to_MERS_Junctions <- nrow(df_SARS2_to_MERS)
df_summary$SARS2_to_MERS_Depth <- sum(df_SARS2_to_MERS$Depth)
#Make new dataframe for aggregate data. Uncomment below ONLY for first sample. Leave commented for subsequent samples.
#df_aggregate <- data.frame(row.names=1)
df_aggregate <- rbind(df_aggregate, df_summary)
#Output file. Uncomment below when all samples added to dataframe.
write.table(df_aggregate, file = "/Users/gribblj/Dropbox/Perlman_MERS_SARS2/RNAseq/Junction_Files/ViReMa_summary.txt", row.names = FALSE, sep = "\t")
