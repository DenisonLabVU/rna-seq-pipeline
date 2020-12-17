#!/usr/bin/env Rscript
print("Usage: Lofreq.R input_file input_coverage_file output_variant_freq output_summary output_transitions output_transversions")
args <- commandArgs(trailingOnly = TRUE)
library(tidyverse)
library(dplyr)
df <- read.table(args[1], header = FALSE)
#Separate column fields
df <- df %>% separate(V8, c("DP", "AF", "SB", "DP4"), sep=";")
#Remove identifies from string to leave just the numbers
df <- df %>% mutate(DP = gsub("DP=", "", DP))
df <- df %>% mutate(AF = gsub("AF=", "", AF))
df <- df %>% mutate(SB = gsub("SB=", "", SB))
df <- df %>% mutate(DP4 = gsub("DP4=", "", DP4))
#Further slice depth fields and remove original column
df <- df %>% separate(DP4, c("Ref_F", "Ref_R", "Variant_F", "Variant_R"), sep=",", remove = TRUE)
#Isolate only columns necessary
df <- df[c(1,2,4,5,11,12,13,14)]
df <- df %>% rename(c("Genome" = V1, "Position" = V2, "Reference" = V4, "Variant" = V5))
#Calculate total depth and variant frequency
df$Total = as.numeric(df$Ref_F) + as.numeric(df$Ref_R) + as.numeric(df$Variant_F) + as.numeric(df$Variant_R)
df$Variant_Total = as.numeric(df$Variant_F) + as.numeric(df$Variant_R)
df$Frequency = df$Variant_Total / df$Total
#Slice transition data
df_transitions <- data.frame()
df_transitions <- df %>% filter(Reference == "A" & Variant == "G") %>% rbind(., df_transitions)
df_transitions <- df %>% filter(Reference == "C" & Variant == "T") %>% rbind(., df_transitions)
df_transitions <- df %>% filter(Reference == "G" & Variant == "A") %>% rbind(., df_transitions)
df_transitions <- df %>% filter(Reference == "T" & Variant == "C") %>% rbind(., df_transitions)
#Slice transversion data
df_transversions <- data.frame()
df_transversions <- df %>% filter(Reference =="A" & Variant == "C") %>% rbind(., df_transversions)
df_transversions <- df %>% filter(Reference =="C" & Variant == "A") %>% rbind(., df_transversions)
df_transversions <- df %>% filter(Reference =="A" & Variant == "T") %>% rbind(., df_transversions)
df_transversions <- df %>% filter(Reference =="T" & Variant == "A") %>% rbind(., df_transversions)
df_transversions <- df %>% filter(Reference =="C" & Variant == "G") %>% rbind(., df_transversions)
df_transversions <- df %>% filter(Reference =="G" & Variant == "C") %>% rbind(., df_transversions)
df_transversions <- df %>% filter(Reference =="G" & Variant == "T") %>% rbind(., df_transversions)
df_transversions <- df %>% filter(Reference =="T" & Variant == "G") %>% rbind(., df_transversions)
#Calculate summary transition data
df_summary <- data.frame()
df_summary <- df_transitions %>% summarise(Transitions_Total = sum(Variant_Total))
df_summary <- df_transitions %>% filter(Reference == "A" & Variant == "G") %>% summarise(A_to_G = sum(Variant_Total)) %>% cbind(., df_summary)
df_summary <- df_transitions %>% filter(Reference == "G" & Variant == "A") %>% summarise(G_to_A = sum(Variant_Total)) %>% cbind(., df_summary)
df_summary <- df_transitions %>% filter(Reference == "C" & Variant == "T") %>% summarise(C_to_T = sum(Variant_Total)) %>% cbind(., df_summary)
df_summary <- df_transitions %>% filter(Reference == "T" & Variant == "C") %>% summarise(T_to_C = sum(Variant_Total)) %>% cbind(., df_summary)
#Calculate summary transversion data
df_summary <- df_transversions %>% summarise(Transversions_Total = sum(Variant_Total)) %>% cbind(., df_summary)
df_summary <- df_transversions %>% filter(Reference == "A" & Variant == "C") %>% summarise(A_to_C = sum(Variant_Total)) %>% cbind(., df_summary)
df_summary <- df_transitions %>% filter(Reference == "C" & Variant == "A") %>% summarise(C_to_A = sum(Variant_Total)) %>% cbind(., df_summary)
df_summary <- df_transitions %>% filter(Reference == "A" & Variant == "T") %>% summarise(A_to_T = sum(Variant_Total)) %>% cbind(., df_summary)
df_summary <- df_transitions %>% filter(Reference == "T" & Variant == "A") %>% summarise(T_to_A = sum(Variant_Total)) %>% cbind(., df_summary)
df_summary <- df_transversions %>% filter(Reference == "G" & Variant == "C") %>% summarise(G_to_C = sum(Variant_Total)) %>% cbind(., df_summary)
df_summary <- df_transitions %>% filter(Reference == "C" & Variant == "G") %>% summarise(C_to_G = sum(Variant_Total)) %>% cbind(., df_summary)
df_summary <- df_transitions %>% filter(Reference == "G" & Variant == "T") %>% summarise(G_to_T = sum(Variant_Total)) %>% cbind(., df_summary)
df_summary <- df_transitions %>% filter(Reference == "T" & Variant == "G") %>% summarise(T_to_G = sum(Variant_Total)) %>% cbind(., df_summary)
#Calculate overall statistics
df_coverage <- read.table(args[2], header = FALSE)
df_coverage <- df_coverage %>% rename(c("Genome" = V1, "Position" = V2, "Depth" = V3))
df_summary <- df_coverage %>% summarise(Total_Mapped = sum(Depth)) %>% cbind(., df_summary)
df_summary <- df %>% summarise(Total_Variant = sum(Variant_Total)) %>% cbind(., df_summary)
df_summary$Mutation_Freq <- (df_summary$Total_Variant / df_summary$Total_Mapped) * 10000
df_summary$Transitions_Freq <- df_summary$Transitions_Total / df_summary$Total_Mapped
df_summary$Transitions_Percent <- (df_summary$Transitions_Total / df_summary$Total_Variant) * 100
df_summary$Transversions_Freq <- df_summary$Transversions_Total / df_summary$Total_Mapped
df_summary$Transversions_Percent <- (df_summary$Transversions_Total / df_summary$Total_Variant) * 100
df_summary$A_to_C_freq <- df_summary$A_to_C / df_summary$Total_Mapped
df_summary$A_to_C_percent <- (df_summary$A_to_C / df_summary$Total_Variant) * 100
df_summary$C_to_A_freq <- df_summary$C_to_A / df_summary$Total_Mapped
df_summary$C_to_A_percent <- (df_summary$C_to_A / df_summary$Total_Variant) * 100
df_summary$A_to_G_freq <- df_summary$A_to_G / df_summary$Total_Mapped
df_summary$A_to_G_percent <- (df_summary$A_to_G / df_summary$Total_Variant) * 100
df_summary$G_to_A_freq <- df_summary$G_to_A / df_summary$Total_Mapped
df_summary$G_to_A_percent <- (df_summary$G_to_A / df_summary$Total_Variant) * 100
df_summary$A_to_T_freq <- df_summary$A_to_T / df_summary$Total_Mapped
df_summary$A_to_T_percent <- (df_summary$A_to_T / df_summary$Total_Variant) * 100
df_summary$T_to_A_freq <- df_summary$T_to_A / df_summary$Total_Mapped
df_summary$T_to_A_percent <- (df_summary$T_to_A / df_summary$Total_Variant) * 100
df_summary$C_to_G_freq <- df_summary$C_to_G / df_summary$Total_Mapped
df_summary$C_to_G_percent <- (df_summary$C_to_G / df_summary$Total_Variant) * 100
df_summary$G_to_C_freq <- df_summary$G_to_C / df_summary$Total_Mapped
df_summary$G_to_C_percent <- (df_summary$G_to_C / df_summary$Total_Variant) * 100
df_summary$C_to_T_freq <- df_summary$C_to_T / df_summary$Total_Mapped
df_summary$C_to_T_percent <- (df_summary$C_to_T / df_summary$Total_Variant) * 100
df_summary$T_to_C_freq <- df_summary$T_to_C / df_summary$Total_Mapped
df_summary$T_to_C_percent <- (df_summary$T_to_C / df_summary$Total_Variant) * 100
df_summary$G_to_T_freq <- df_summary$G_to_T / df_summary$Total_Mapped
df_summary$G_to_T_percent <- (df_summary$G_to_T / df_summary$Total_Variant) * 100
df_summary$T_to_G_freq <- df_summary$T_to_G / df_summary$Total_Mapped
df_summary$T_to_G_percent <- (df_summary$T_to_G / df_summary$Total_Variant) * 100
#Write table with frequencies
write.table(df, args[3], sep = "\t", row.names = FALSE)
write.table(df_summary, args[4], sep = "\t", row.names = FALSE)
write.table(df_transitions, args[5], sep = "\t", row.names = FALSE)
write.table(df_transversions, args[6], sep = "\t", row.names = FALSE)
