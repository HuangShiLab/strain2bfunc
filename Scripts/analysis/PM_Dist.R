#!/usr/bin/env Rscript
#################################################################
# Function: Distance calculation
# Call: Rscript PM_Dist.R -i abund_file [-d dist_type -t taxonomy_file -o outfile]
# Updated at Apr. 12, 2024
# Updated by ZHANG Yufeng <yfz96@connect.hku.hk>, HUANG Shi <shihuang@hku.hk>
# Faculty of Dentistry, HKU
#################################################################

## install necessary libraries
p <- c("optparse","vegan","ade4","ggplot2","grid","RColorBrewer", "readr", "reshape2", "dplyr")
usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep=TRUE, repos="http://cran.us.r-project.org/")
  suppressWarnings(suppressMessages(invisible(require(p, character.only=TRUE))))
}
invisible(lapply(p, usePackage))

## clean R environment
rm(list = ls())
setwd('./')

## parsing arguments
args <- commandArgs(trailingOnly=TRUE)

# make option list and parse command line
option_list <- list(
  make_option(c("-i", "--abund_file"), type="character", help="Input feature table with relative abundance (*.Abd) [Required]"),
  make_option(c("-d", "--dist_type"), type="character", default="bray", help="Input the distance algorithm, e.g., \"bray\" for Bray-Curtis distance, \"taxUMAP\" for TaxUMAP distance, \"euclidean\" for Euclidean distance, \"jaccard\" for Jaccard distance, etc.\" [default %default]"),
  make_option(c("-t", "--taxonomy_file"), type="character", help="Input the taxonomy annotation file [Required when dist_type is equal to \"taxUMAP\"]"),
  make_option(c("-o", "--outfile"), type="character", default='dist.txt', help="Output distance matrix [default %default]")
)
opts <- parse_args(OptionParser(option_list=option_list), args=args)

# paramenter checking
if(is.null(opts$abund_file)) stop('Please input a feature table (*.Abd)')

# import abundance file
abund_orig <- read.table(file=opts$abund_file, sep="\t", header=TRUE, row.names=1)
dist_type <- opts$dist_type

if(dist_type == "taxUMAP") {
  if(is.null(opts$taxonomy_file)) stop('Please input the taxonomy annotation file to calculate taxUMAP distance')
  
  strain_dist <- as.matrix(vegdist(t(abund_orig), method = "bray"))
  
  taxonomy <- read_delim(opts$taxonomy_file, delim = "\t", col_names = TRUE, comment = "", show_col_types = FALSE)
  taxonomy <- as.data.frame(taxonomy)
  
  match_indices <- match(rownames(abund_orig), taxonomy$strain)
  abund_taxonomy <- taxonomy[match_indices, c("kingdom", "phylum", "class", "order", "family", "genus", "specie", "strain")]
  abund_orig <- cbind(abund_taxonomy, abund_orig)
  
  melt_abund <- melt(abund_orig)
  
  family_abd <- melt_abund %>%
                group_by(variable, family) %>%
                summarise(family_abd = sum(value))
  family_abd <- dcast(family_abd, variable ~ family)
  rownames(family_abd) <- family_abd$variable
  family_abd <- family_abd[, -1, drop = FALSE]
  family_dist <- as.matrix(vegdist(family_abd, method = "bray"))
  
  phylum_abd <- melt_abund %>%
    group_by(variable, phylum) %>%
    summarise(phylum_abd = sum(value))
  phylum_abd <- dcast(phylum_abd, variable ~ phylum)
  rownames(phylum_abd) <- phylum_abd$variable
  phylum_abd <- phylum_abd[, -1, drop = FALSE]
  phylum_dist <- as.matrix(vegdist(phylum_abd, method = "bray"))

  dist <- strain_dist + family_dist + phylum_dist
} else {
  abund_orig <- t(abund_orig)
  dist <- vegdist(abund_orig, method = dist_type)
  dist <- as.matrix(dist)
}

write.table(dist, opts$outfile, sep = "\t", quote = F, row.names = T, col.names = NA)

