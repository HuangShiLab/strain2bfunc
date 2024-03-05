#!/usr/bin/env Rscript
#################################################################
# Function: PCA
# Call: Rscript PM_pca.R -m meta_data -i abund_file [-l T/F -o outfile -a axesfile]
# R packages used: optparse vegan ggplot2 grid
# update 2: 2020-10-18, 
# Authors: Yuzhu Chen, Yufeng Zhang, Zheng Sun, Yanhai Gong,Wang Honglei, Xiaoquan Su
# Updated at Aug. 20, 2021
# Updated by Yuzhu Chen
# Bioinformatics Group, College of Computer Science & Technology, Qingdao University
#################################################################

## install necessary libraries
p <- c("optparse","vegan","ade4","ggplot2","grid","RColorBrewer")
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
  make_option(c("-d", "--dist_type"), type="character", default="bray", help="Input the distance algorithm, e.g., \"bray\" for Bray-Curtis distance, \"euclidean\" for Euclidean distance, \"jaccard\" for Jaccard distance, etc.\" [default %default]"),
  make_option(c("-o", "--outfile"), type="character", default='dist.txt', help="Output distance matrix [default %default]")
)
opts <- parse_args(OptionParser(option_list=option_list), args=args)

# paramenter checking
if(is.null(opts$abund_file)) stop('Please input a feature table (*.Abd)')

# import abundance file
abund_orig <- read.table(file=opts$abund_file, header=TRUE, row.names=1)
abund_orig <- t(abund_orig)
dist_type <- opts$dist_type

dist <- vegdist(abund_orig, method = dist_type)
dist <- as.matrix(dist)

write.table(dist, opts$outfile, sep = "\t", quote = F, row.names = T, col.names = NA)

