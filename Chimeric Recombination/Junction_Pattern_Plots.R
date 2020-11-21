library(ggplot2)
data <- read.table("/Users/gribblj/Dropbox/Perlman_MERS_SARS2/RNAseq/0-01B_MERS_to_SARS2_junctions.txt", header = TRUE)
data <- data[order(data$Depth), ]
data$Total = sum(data$Depth)
data$Frequency = data$Depth / data$Total
data$logFreq = log10(data$Frequency)
data_graph <- ggplot(data, aes(Stop, Start, colour = logFreq, alpha = logFreq)) + geom_point(size = 2) + xlim(0, 30500) + ylim(0, 30500) + theme_linedraw(base_size = 20) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.title = element_text(size = 16), legend.text = element_text(size = 14), legend.justification = "center") + labs(x = "Stop Position (nt)", y = "Start Position (nt)") + guides(alpha = FALSE)
print(data_graph + scale_color_gradientn(colours = rainbow(6)))
ggsave(filename = "0-01B_MERS_to_SARS2_junctions.tiff", plot = last_plot(), device = "tiff", path = "/Users/gribblj/Dropbox/Perlman_MERS_SARS2/RNAseq/Junction_Plots/", scale = 1, width = 6, height = 4, units = "in", dpi = 1200, limitsize = TRUE)
