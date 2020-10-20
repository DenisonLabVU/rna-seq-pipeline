library(dplyr)
df <- read.table("/Users/jennifergribble/Dropbox/P250_recombination/Passage_Populations/XN_passages_virion/Junction_Files/XNP3A_forward_junctions.txt", header = TRUE)
df_TRSL <- filter(df, between(df$Start, 31, 103))
df_sgmRNA2 <- filter(df_TRSL, between(df_TRSL$Stop, 21713, 21785))
df_sgmRNA3 <- filter(df_TRSL, between(df_TRSL$Stop, 23888, 23960))
df_sgmRNA4 <- filter(df_TRSL, between(df_TRSL$Stop, 27901, 27973))
df_sgmRNA5 <- filter(df_TRSL, between(df_TRSL$Stop, 28284, 28356))
df_sgmRNA6 <- filter(df_TRSL, between(df_TRSL$Stop, 28924, 28996))
df_sgmRNA7 <- filter(df_TRSL, between(df_TRSL$Stop, 29621, 29693))
#Add column identifying sgmRNA species
df_sgmRNA2$Type <- "sgmRNA2"
df_sgmRNA3$Type <- "sgmRNA3"
df_sgmRNA4$Type <- "sgmRNA4"
df_sgmRNA5$Type <- "sgmRNA5"
df_sgmRNA6$Type <- "sgmRNA6"
df_sgmRNA7$Type <- "sgmRNA7"
#Slice canonical sgmRNA species
sgmRNA2_canonical <- df_sgmRNA2 %>% arrange(desc(df_sgmRNA2$Depth)) %>% slice(1)
sgmRNA3_canonical <- df_sgmRNA3 %>% arrange(desc(df_sgmRNA3$Depth)) %>% slice(1)
sgmRNA4_canonical <- df_sgmRNA4 %>% arrange(desc(df_sgmRNA4$Depth)) %>% slice(1)
sgmRNA5_canonical <- df_sgmRNA5 %>% arrange(desc(df_sgmRNA5$Depth)) %>% slice(1)
sgmRNA6_canonical <- df_sgmRNA6 %>% arrange(desc(df_sgmRNA6$Depth)) %>% slice(1)
sgmRNA7_canonical <- df_sgmRNA7 %>% arrange(desc(df_sgmRNA7$Depth)) %>% slice(1)
#Create concatenated dataframe of canonical sgmRNAs
df_canonical <- rbind(sgmRNA2_canonical, sgmRNA3_canonical, sgmRNA4_canonical, sgmRNA5_canonical, sgmRNA6_canonical, sgmRNA7_canonical)
df_canonical$Total <- sum(df_canonical$Depth)
df_sgmRNA <- rbind(df_sgmRNA2, df_sgmRNA3, df_sgmRNA4, df_sgmRNA5, df_sgmRNA6, df_sgmRNA7)
df_sgmRNA$Total <- sum(df_sgmRNA$Depth)
#Print list of alternative sgmRNAs
df_alternative <- anti_join(df_sgmRNA, df_canonical, by = "Depth")
df_alternative$Total <- sum(df_alternative$Depth)
df_alt_summary <- df_alternative %>% group_by(Type) %>% summarise(Sum = sum(Depth))
#Write tables
write.table(df_canonical, file = "/Users/jennifergribble/Dropbox/P250_recombination/Passage_Populations/sgmRNAs_DVGs/XNP3A_canonical_sgmRNAs.txt", sep = "\t", row.names = FALSE)
write.table(df_alternative, file = "/Users/jennifergribble/Dropbox/P250_recombination/Passage_Populations/sgmRNAs_DVGs/XNP3A_alternative_sgmRNAs.txt", sep = "\t", row.names = FALSE)
write.table(df_sgmRNA, file = "/Users/jennifergribble/Dropbox/P250_recombination/Passage_Populations/sgmRNAs_DVGs/XNP3A_total_sgmRNAs.txt", sep = "\t", row.names = FALSE)
