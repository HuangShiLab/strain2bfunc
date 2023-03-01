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
        make_option(c("-l", "--sample_list_file"), type="character", help="Input the samples list file [Required]"),
        make_option(c("-s", "--species_list_file"), type="character", help="Input the identified species list [Required]"),
        make_option(c("-o", "--output_path"), type="character", default='.', help="Output file [default %default]")
)

opts <- parse_args(OptionParser(option_list=option_list), args=args)

# Error checking
if(is.null(opts$sample_list_file)) stop('Please input a sample list file')
if(is.null(opts$species_abd_file)) stop('Please input an identified species list')

source("/lustre1/g/aos_shihuang/Strain2b/src/2/utility.R")

Sample_List_Pipeline(opts$sample_list_file, opts$species_list_file, opts$output_path)
