#!/usr/bin/env Rscript
#################################################################
# Function: Merge the strain-level abundance table of different species
# Call: Rscript PM_pca.R -m meta_data -i abund_file [-l T/F -o outfile -a axesfile]
# R packages used: optparse vegan ggplot2 grid
# update 2: 2020-10-18, 
# Authors: Yuzhu Chen, Yufeng Zhang, Zheng Sun, Yanhai Gong,Wang Honglei, Xiaoquan Su
# Updated at Aug. 20, 2021
# Updated by Yuzhu Chen
# Bioinformatics Group, College of Computer Science & Technology, Qingdao University
#################################################################

## install necessary libraries
p <- c("optparse","vegan","ade4","ggplot2","grid","RColorBrewer", "readr")
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
  make_option(c("-s", "--species_list_file"), type="character", help="Input the species list file [Required]"),
  make_option(c("-p", "--prefix"), type = "character", help="Prefix of the path of the strain-level abundance file [Required]"),
  make_option(c("-o", "--outfile"), type="character", default='merged_strain_level_abd.txt', help="Output merged strain-level abundance table [default %default]")
)
opts <- parse_args(OptionParser(option_list=option_list), args=args)

# paramenter checking
if(is.null(opts$species_list_file)) stop('Please input a species list file')
if(is.null(opts$prefix)) stop('Please input the prefix of the path of the strain-level abundance file')

species_list <- read.table(opts$species_list_file, sep = "\t", header = F)

prefix <- opts$prefix

strain_file_path <- paste(prefix, species_list[, 1], "strain_level_abd.txt", sep = "/")
strain_file_path <- strain_file_path[file.exists(strain_file_path)]

if(length(strain_file_path) == 0) stop("This is no strain-level abundance table")

results = data.frame()
for(i in 1:length(strain_file_path)) {
  table <- read.table(strain_file_path[i], sep = "\t", header = T, row.names = 1)
  table <- table[, order(colnames(table)), drop = FALSE]
  if(nrow(results) == 0) {
    results = table
  } else {
    if(! identical(colnames(results), colnames(table))) {
      next
    }
    else {
      results <- rbind(results, table)
    }
  }
}

results <- as.data.frame(apply(results, 2, function(x) x/sum(x)))

write.table(results, opts$outfile, sep = "\t", quote = F, row.names = T, col.names = NA)


