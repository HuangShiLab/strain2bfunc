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
	      make_option(c("-c", "--cnm_file"), type="character", 
                    help="The path of copy number matrix [Required]"),
	      make_option(c("-r", "--rcm_file"), type="integer", default = 0,
                    help="The path of read count matrix [Required]"),
        make_option(c("-o", "--output_path"), type="character", default='.', 
                    help="Output file [default %default]")
)

opts <- parse_args(OptionParser(option_list=option_list), args=args)

Env <-Sys.getenv("Strain2bFunc")
if(nchar(Env)<1){
  cat('Please set the environment variable \"Strain2bFunc\" to the directory\n')
}
source(paste0(Env, "/Scripts/strain2b/utility.R"))

if (is.null(opts$cnm_file)) {
  stop("Please specify a copy number matrix!")
} 
if (is.null(opts$rcm_file)) {
  stop("Please specify a read count matrix!")
}

cnm <- read.table(opts$cnm_file, sep = "\t", header = T, row.names = 1)
tags_count <- read.table(opts$rcm_file, sep = "\t", header = T, row.names = 1, check.names = F, comment.char = "")

# cnm_file <- "../test_data/2bRAD_tag_matrix.xls"
# rcm_file <- "../test_data/read_count_matrix.txt"
# cnm <- read.table(cnm_file, sep = "\t", header = T, row.names = 1)
# tags_count <- read.table(rcm_file, sep = "\t", header = T, row.names = 1, check.names = F, comment.char = "")

matrix <- Filter_CNM(cnm, tags_count)

tags_count <- matrix$tags_count
cnm <- matrix$cnm

strain_abd_matrix <- Strain_Level_Profiling(tags_count, cnm)
write.table(strain_abd_matrix, opts$output_path, sep = "\t", row.names = T, col.names = NA, quote = F)

