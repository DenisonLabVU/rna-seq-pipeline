library(ggplot2)
data <- read.table("/Users/jennifergribble/Dropbox/Specific_Infectivity_paper/SARSCoV2_NHC_RNAseq/Junction_Files/SARS2-0125-2_junctions.txt", header = TRUE)
data <- data[order(data$Depth), ]
data_forward <- data[which(data$Start < data$Stop), ]
data_forward$Total = sum(data_forward$Depth)
data_forward$Frequency = data_forward$Depth / data_forward$Total
data_forward$logFreq = log10(data_forward$Frequency)
data_graph <- ggplot(data_forward, aes(Stop, Start, colour = logFreq, alpha = logFreq)) + geom_point(size = 2) + theme_linedraw(base_size = 20) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.title = element_text(size = 16), legend.text = element_text(size = 14), legend.justification = "center") + labs(x = "Stop Position (nt)", y = "Start Position (nt)") + guides(alpha = FALSE)
print(data_graph + scale_color_gradientn(colours = rainbow(6)))
ggsave(filename = "SARS2-0125-2_junctions.tiff", plot = last_plot(), device = "tiff", path = "/Users/jennifergribble/Dropbox/Specific_Infectivity_paper/SARSCoV2_NHC_RNAseq/Junction_Plots/", scale = 1, width = 6, height = 4, units = "in", dpi = 1200, limitsize = TRUE)
write.table(data_forward, file="/Users/jennifergribble/Dropbox/Specific_Infectivity_paper/SARSCoV2_NHC_RNAseq/Junction_Files/SARS2-0125-2_forward_junctions.txt", sep = "\t", row.names = FALSE)