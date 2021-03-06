library(dplyr)
df <- read.table("/Users/jennifergribble/Dropbox/Specific_Infectivity_paper/SARSCoV2_NHC_RNAseq/Junction_Files/SARS2-05-3_forward_junctions.txt", header = TRUE)
df_TRSL <- filter(df, between(df$Start, 39, 106))
df_sgmRNA2 <- filter(df_TRSL, between(df_TRSL$Stop, 21525, 21592))
df_sgmRNA3 <- filter(df_TRSL, between(df_TRSL$Stop, 25354, 25421))
df_sgmRNA4 <- filter(df_TRSL, between(df_TRSL$Stop, 26206, 26273))
df_sgmRNA5 <- filter(df_TRSL, between(df_TRSL$Stop, 26442, 26509))
df_sgmRNA6 <- filter(df_TRSL, between(df_TRSL$Stop, 27010, 27077))
df_sgmRNA7 <- filter(df_TRSL, between(df_TRSL$Stop, 27357, 27424))
df_sgmRNA8 <- filter(df_TRSL, between(df_TRSL$Stop, 27857, 27924))
df_sgmRNA9 <- filter(df_TRSL, between(df_TRSL$Stop, 28229, 28296))
#Add column identifying sgmRNA species
df_sgmRNA2$Type <- "sgmRNA2"
df_sgmRNA3$Type <- "sgmRNA3"
df_sgmRNA4$Type <- "sgmRNA4"
df_sgmRNA5$Type <- "sgmRNA5"
df_sgmRNA6$Type <- "sgmRNA6"
df_sgmRNA7$Type <- "sgmRNA7"
df_sgmRNA8$Type <- "sgmRNA8"
df_sgmRNA9$Type <- "sgmRNA9"
#Slice canonical sgmRNA species
sgmRNA2_canonical <- df_sgmRNA2 %>% arrange(desc(df_sgmRNA2$Depth)) %>% slice(1)
sgmRNA3_canonical <- df_sgmRNA3 %>% arrange(desc(df_sgmRNA3$Depth)) %>% slice(1)
sgmRNA4_canonical <- df_sgmRNA4 %>% arrange(desc(df_sgmRNA4$Depth)) %>% slice(1)
sgmRNA5_canonical <- df_sgmRNA5 %>% arrange(desc(df_sgmRNA5$Depth)) %>% slice(1)
sgmRNA6_canonical <- df_sgmRNA6 %>% arrange(desc(df_sgmRNA6$Depth)) %>% slice(1)
sgmRNA7_canonical <- df_sgmRNA7 %>% arrange(desc(df_sgmRNA7$Depth)) %>% slice(1)
sgmRNA8_canonical <- df_sgmRNA8 %>% arrange(desc(df_sgmRNA8$Depth)) %>% slice(1)
sgmRNA9_canonical <- df_sgmRNA9 %>% arrange(desc(df_sgmRNA9$Depth)) %>% slice(1)
#Create concatenated dataframe of canonical sgmRNAs
df_canonical <- rbind(sgmRNA2_canonical, sgmRNA3_canonical, sgmRNA4_canonical, sgmRNA5_canonical, sgmRNA6_canonical, sgmRNA7_canonical, sgmRNA8_canonical, sgmRNA9_canonical)
df_canonical$Total <- sum(df_canonical$Depth)
df_sgmRNA <- rbind(df_sgmRNA2, df_sgmRNA3, df_sgmRNA4, df_sgmRNA5, df_sgmRNA6, df_sgmRNA7, df_sgmRNA8, df_sgmRNA9)
df_sgmRNA$Total <- sum(df_sgmRNA$Depth)
#Print list of alternative sgmRNAs
df_alternative <- anti_join(df_sgmRNA, df_canonical, by = "Depth")
df_alternative$Total <- sum(df_alternative$Depth)
df_alt_summary <- df_alternative %>% group_by(Type) %>% summarise(Sum = sum(Depth))
#Slice DVGs and turn format into BED
df_DVG <- anti_join(df, df_sgmRNA, by = c("Start", "Stop"))
df_DVG$Duplication <- "Duplication"
df_DVG$Strand <- "+"
df_DVG$Start2 <- df_DVG$Start
df_DVG$Stop2 <- df_DVG$Stop
df_DVG <- df_DVG[c(1,2,3,8,4,9,10,11)]
#Write tables
write.table(df_canonical, file = "/Users/jennifergribble/Dropbox/Specific_Infectivity_paper/SARSCoV2_NHC_RNAseq/sgmRNAS_DVGs/05-3_canonical_sgmRNAs.txt", sep = "\t", row.names = FALSE)
write.table(df_alternative, file = "/Users/jennifergribble/Dropbox/Specific_Infectivity_paper/SARSCoV2_NHC_RNAseq/sgmRNAS_DVGs/05-3_alternative_sgmRNAs.txt", sep = "\t", row.names = FALSE)
write.table(df_sgmRNA, file = "/Users/jennifergribble/Dropbox/Specific_Infectivity_paper/SARSCoV2_NHC_RNAseq/sgmRNAS_DVGs/05-3_total_sgmRNAs.txt", sep = "\t", row.names = FALSE)
write.table(df_alt_summary, file = "/Users/jennifergribble/Dropbox/Specific_Infectivity_paper/SARSCoV2_NHC_RNAseq/sgmRNAS_DVGs/05-3_alt_sgmRNA_summary.txt", sep = "\t", row.names = FALSE)
write.table(df_DVG, file = "//Users/jennifergribble/Dropbox/Specific_Infectivity_paper/SARSCoV2_NHC_RNAseq/sgmRNAS_DVGs/05-3_DVGs.bed.txt", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)