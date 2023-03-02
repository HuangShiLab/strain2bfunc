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
        make_option(c("-g", "--ground_truth_file"), type="character", help="Input the ground truth file [Required]"),
        make_option(c("-s", "--predicted_abd_file"), type="character", help="Input the predicted abundance file [Required]"),
        make_option(c("-o", "--output_file"), type="character", default='abd.txt', help="Output file [default %default]")
)

opts <- parse_args(OptionParser(option_list=option_list), args=args)

# Error checking
if(is.null(opts$ground_truth_file)) stop('Please input a sample list file')
if(is.null(opts$predicted_abd_file)) stop('Please input an identified species list')
source("/lustre1/g/aos_shihuang/Strain2b/src/2/composition.R")

Merge_Two_Profiling_Matrix <- function(abd1, abd2) {
        result <- merge(abd1, abd2, by = 1, all = TRUE)
        result[is.na(result)] <- 0
        return (result)
}

ground_truth_abd <- read.table(opts$ground_truth_file, sep = "\t", header = T, row.names = 1)
predicted_abd <- read.table(opts$predicted_abd_file, sep = "\t", header = T, row.names = 1)

result <- Merge_Two_Profiling_Matrix(ground_truth_abd, predicted_abd)
write.table(result, opts$output_file, sep = "\t", quote = F, row.names = 1, col.names = NA)
