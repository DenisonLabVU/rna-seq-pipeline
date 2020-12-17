#!/usr/bin/env Rscript
print("Usage: Junction_Pattern_Plots.R input_file output_tiff_name output_tiff_path output_file")
args <- commandArgs(trailingOnly = TRUE)
library(ggplot2)
data <- read.table(args[1], header = TRUE)
data <- data[order(data$Depth), ]
data_forward <- data[which(data$Start < data$Stop), ]
data_forward$Total = sum(data_forward$Depth)
data_forward$Frequency = data_forward$Depth / data_forward$Total
data_forward$logFreq = log10(data_forward$Frequency)
data_graph <- ggplot(data_forward, aes(Stop, Start, colour = logFreq, alpha = logFreq)) + geom_point(size = 2) + theme_linedraw(base_size = 20) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.title = element_text(size = 16), legend.text = element_text(size = 14), legend.justification = "center") + labs(x = "Stop Position (nt)", y = "Start Position (nt)") + guides(alpha = FALSE)
print(data_graph + scale_color_gradientn(colours = rainbow(6)))
ggsave(filename = args[2], plot = last_plot(), device = "tiff", path = args[3], scale = 1, width = 6, height = 4, units = "in", dpi = 1200, limitsize = TRUE)
write.table(data_forward, file=args[4], sep = "\t", row.names = FALSE)