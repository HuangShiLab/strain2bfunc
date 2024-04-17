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
	      make_option(c("-m", "--mode"), type="integer", default = 0,
                    help="choices=c(0, 1), 0 for profiling based on the 2bGDB, 1 for profiling based on a customized copy number matrix [default %default]"),
	      make_option(c("-d", "--cnm_matrix_dir"), type="character",
                          help="The directory of copy number matrix [Required when mode equals to 0]"),
	      make_option(c("-s", "--species_list_file"), type="character", 
	                  help="Input the identified species list [Required when mode is equal to 0]"),
	      make_option(c("-f", "--cnm_file"), type="character", 
	                  help="The path of copy number matrix [Required when mode equals to 1]"),
        make_option(c("-o", "--output_path"), type="character", default='.', 
                    help="Output file [default %default]")
)

opts <- parse_args(OptionParser(option_list=option_list), args=args)

Env <-Sys.getenv("Strain2bFunc")
if(nchar(Env)<1){
  cat('Please set the environment variable \"Strain2bFunc\" to the directory\n')
 }

source(paste0(Env, "/Scripts/strain2b/utility.R"))

if (is.null(opts$sample_list_file)) stop('Please input a sample list file')
if (opts$mode == 0) {
  if (is.null(opts$cnm_matrix_dir)) {
    stop("Please specify a directory of copy number matrix!")
  } else {
    Sample_List_Pipeline(opts$sample_list_file, opts$species_list_file, opts$output_path, cnm_matrix_dir = opts$cnm_matrix_dir, mode = opts$mode)
  }
} else if (opts$mode == 1) {
  if (is.null(opts$cnm_file)) {
    stop("Please specify a copy number matrix!")
  } else {
    Sample_List_Pipeline(opts$sample_list_file, opts$species_list_file, opts$output_path, cnm_file = opts$cnm_file, mode = opts$mode)
  }
} else {
  stop("Invalid mode option. Please specify 0 or 1.")
}

