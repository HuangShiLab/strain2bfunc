options(warn=-1)
## install necessary libraries
p <- c("optparse", "Matrix")

usePackage <- function(p){
        if (!is.element(p, installed.packages()[,1]))
                install.packages(p, dep=TRUE, repos="http://cran.us.r-project.org/")
        suppressWarnings(suppressMessages(invisible(require(p, character.only=TRUE))))
}
invisible(lapply(p, usePackage))

## parsing arguments
args <- commandArgs(trailingOnly=TRUE)

# make option list and parse command line
option_list <- list(
        make_option(c("-l", "--sample_list_file"), type="character", 
                    help="Input the samples list file [Required]"),
        make_option(c("-a", "--species_abundance_file"), type="character", 
                    help="Input the species abundance file [Required when mode is equal to 0]"),
        make_option(c("-t", "--threshold"), type="numeric", default = 0.001,
                    help="The threshold of species abundance for strain-level identification, [Required when mode is equal to 0, default %default]"),
        make_option(c("-s", "--species_list_file"), type="character", 
                    help="Input the identified species list [Required when mode is equal to 1]"),
        make_option(c("-m", "--mode"), type="integer", default = 0, 
                    help="choices=c(0, 1), 0 for global strain level profiling, 1 for strain level profiling for identified species [default %default]"),
        make_option(c("-o", "--output_path"), type="character", default='.', 
                    help="Output file [default %default]")
)

opts <- parse_args(OptionParser(option_list=option_list), args=args)

source("utility.R")

if(is.null(opts$sample_list_file)) stop('Please input a sample list file')
if (opt$mode == 0) {
  if (is.null(opt$species_abundance_file)) {
    stop("Please specify a species abundance table file!")
  } else {
    Sample_List_Pipeline(opts$sample_list_file, opts$species_abundance_file, opts$output_path, mode = 0, threshold = threshold)
  }
} else if (opt$mode == 1) {
  if (is.null(opt$species_list_file)) {
    stop("Please specify a species list file!")
  } else {
    Sample_List_Pipeline(opts$sample_list_file, opts$species_list_file, opts$output_path, mode = 1)
  }
} else {
  stop("Invalid mode option. Please specify 0 or 1.")
}



