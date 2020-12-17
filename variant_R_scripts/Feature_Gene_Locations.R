library(tidyverse)
library(dplyr)
df <- read.table("/Users/jennifergribble/Dropbox/SARS2_RDV_passage/21-9133_variant_freq.txt", header = TRUE)
df <- df %>% select(Genome, Position, Reference, Variant, Frequency)
df$Gene <- if_else(between(df$Position, 0, 266), "5UTR",
                  if_else(between(df$Position, 265, 806), "nsp1",
                          if_else(between(df$Position, 805, 2720), "nsp2",
                                  if_else(between(df$Position, 2719, 8555), "nsp3",
                                          if_else(between(df$Position, 8554, 10055), "nsp4",
                                                  if_else(between(df$Position, 10054, 10973), "nsp5",
                                                          if_else(between(df$Position, 10972, 11843), "nsp6",
                                                                  if_else(between(df$Position, 11842, 12092), "nsp7",
                                                                          if_else(between(df$Position, 12091, 12686), "nsp8",
                                                                                  if_else(between(df$Position, 12685, 13025), "nsp9",
                                                                                          if_else(between(df$Position, 13024, 13442), "nsp10",
                                                                                                  if_else(between(df$Position, 13441, 16237), "nsp12",
                                                                                                          if_else(between(df$Position, 16236, 18040), "nsp13",
                                                                                                                  if_else(between(df$Position, 18039, 19621), "nsp14",
                                                                                                                          if_else(between(df$Position, 19620, 20659), "nsp15",
                                                                                                                                  if_else(between(df$Position, 20658, 21553), "nsp16",
                                                                                                                                          if_else(between(df$Position, 21562, 25385), "S protein",
                                                                                                                                                  if_else(between(df$Position, 25392, 26221), "ORF3a",
                                                                                                                                                          if_else(between(df$Position, 26244, 26473), "E protein",
                                                                                                                                                                  if_else(between(df$Position, 26522, 27192), "M protein",
                                                                                                                                                                          if_else(between(df$Position, 27201, 27388), "ORF6",
                                                                                                                                                                                  if_else(between(df$Position, 27393, 27888), "ORF7ab",
                                                                                                                                                                                          if_else(between(df$Position, 27893, 28260), "ORF8",
                                                                                                                                                                                                  if_else(between(df$Position, 28273, 29534), "N protein",
                                                                                                                                                                                                          if_else(between(df$Position, 29557, 29675), "ORF10",
                                                                                                                                                                                                                  ifelse(between(df$Position, 29674, 29870), "3UTR", "NA"))))))))))))))))))))))))))
write.table(df, "/Users/jennifergribble/Dropbox/SARS2_RDV_passage/21-9133_variants.txt", sep = "\t", row.names = FALSE)
