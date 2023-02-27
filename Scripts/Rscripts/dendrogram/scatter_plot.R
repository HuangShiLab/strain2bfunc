p <- c("ggplot2", "ggpubr", "stringr", "RColorBrewer")
usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep=TRUE, repos="https://cloud.r-project.org/")
  suppressWarnings(suppressMessages(invisible(require(p, character.only=TRUE))))
}
invisible(lapply(p, usePackage))

dist <- read.table("all.dist.txt", sep = "\t", header = T, row.names = 1)
meta <- read.table("Staphylococcus.txt", sep = "\t", header = T, row.names = 1)

out <- gsub("_([^_]+)_([^_]+)$", " \\1_\\2", rownames(dist))
out <- str_split_fixed(out, " ", 2)
dist$V1 <- out[, 1]
dist$V2 <- out[, 2]

dist$group <- paste(meta[dist$V1, "species"], meta[dist$V2, "species"], sep = "_v.s._")
dist[which(dist$group == "SE_v.s._SA"), "group"] = "SA_v.s._SE"

max_value <- max(dist[c("WMS", "BcgI", "rRNA")])

plot <- ggplot(dist, aes(x=WMS, y=BcgI)) +
  geom_point(mapping = aes(color = group), shape = 21) +
  scale_colour_manual(values = brewer.pal(8, "Pastel2")) +
  geom_smooth(method=lm) +
  stat_cor(data=dist, method = "spearman") +
  xlim(0, max_value) +
  ylim(0, max_value) +
  xlab("Whole genome") +
  ylab("2bRAD (BcgI)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

ggsave(filename=paste("genome_vs_bcgi.pdf", sep=""), plot=plot, width=3.5, height=2.3)

plot <- ggplot(dist, aes(x=WMS, y=rRNA)) + 
  geom_point(mapping = aes(color = group), shape = 21) + 
  scale_colour_manual(values = brewer.pal(8, "Pastel2")) +
  geom_smooth(method=lm) +
  stat_cor(data=dist, method = "spearman") +
  xlim(0, max_value) +
  ylim(0, max_value) +
  xlab("Whole genome") +
  ylab("16S rRNA") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

ggsave(filename=paste("genome_vs_16S.pdf", sep=""), plot=plot, width=3.5, height=2.3)
