Env <-Sys.getenv("Strain2bFunc")
if(nchar(Env)<1){
  cat('Please set the environment variable \"Strain2bFunc\" to the directory\n')
}
#print(Env)

source(paste0(Env, "/Scripts/strain2b/composition.R"))


Read_Copynumber_Matrix <- function(species, cnm_matrix_dir) {
  
  if(length(species) == 0) {
    result = NULL
  }
  else {
    cnm_file <- paste0(cnm_matrix_dir, "/", species, ".CNM.xls")
    
    if (!file.exists(cnm_file)) {
      stop(paste0("The copy number matrix of ", species, " does not exist. Please make sure the file exists and try again."))
    }
    
    result = read.table(cnm_file, sep = "\t", header = T)
  }
  
  return (result)
}


Merge_Copynumber_Matrix <- function(all_species, cnm_matrix_dir) { # merge the copynumber matrix of all species
  
  if(length(all_species) == 0) {
    result <- NULL
  }
  else if(length(all_species) == 1) {
    result <- Read_Copynumber_Matrix(all_species[1], cnm_matrix_dir)
  }
  else {
    result <- data.frame()
    
    for (i in 1:length(all_species)) {
      cnm <- Read_Copynumber_Matrix(all_species[i], cnm_matrix_dir)
      
      if (nrow(result) == 0) {
        result <- cnm
      } else {
        result <- merge(result, cnm, by = 1, all = TRUE)
        result[is.na(result)] <- 0
      }
      
    }
  }
  
  rownames(result) <- result[, 1]
  result <- result[ , -1, drop = F]
  return (result)
}

Rename_sample <- function(sample_name) {
# Since sequence names for vsearch need to include file names, 
# but the file names cannot contain ".", 
# the "." in file names need to be replaced with "_"
  sample_name <- gsub("\\.", "_", sample_name)
  sample_name <- gsub("-", "_", sample_name)

  return (sample_name)
}


Rename_Fasta <- function(sample_name, sample_fa, new_sample_fa) {
  fasta <- read.fasta(sample_fa, forceDNAtolower = F)	
  
  sample_name <- Rename_sample(sample_name)
  
  names(fasta) <- paste(sample_name, 1:length(fasta), sep = ".")
  
  write.fasta(sequences = fasta, names = names(fasta), file.out = new_sample_fa)
}


Vsearch <- function(cnm, new_sample_fa, similarity, output_tag_path, tags_count_file) {
  sink(output_tag_path)
  
  for (tag in rownames(cnm)) {
    cat(paste0(">", tag, "\n"))
    cat(paste0(tag, "\n"))
  }
  
  sink()
  
  cmd <- paste0("vsearch --usearch_global ", new_sample_fa, " --db ", output_tag_path, 
                " --id ", similarity, " --iddef 4 --strand both -otutabout " , tags_count_file, 
                " --threads 2") # --threads: number of threads to use, zero for all cores (0)
  
  system(cmd, intern = TRUE)
  
  #file.remove(new_sample_fa) #Delete intermediate result files
  #file.remove(output_tag_path) #Delete intermediate result files
}

corrDist <- function(X){
  csd <- apply(X, 2, sd)
  if(sum(csd == 0) > 0){
    stop("Columns", which(csd == 0), "have no variance! Cannot compute correlation distances")
  }
  CSD <- csd %*% t(csd)
  cm <- colMeans(X)
  CM <- cm %*% t(cm)
  n <- nrow(X)
  XX <- crossprod(X)
  D <- 1 - (as.matrix(XX) - n * CM) / ((n - 1) * CSD)
  return(D)
}

#' @name deduplicate_cnm
#' 
#' @description Check the auto-correlation between columns, and deduplicate the columns correlated with any others from the copy number matrix.
#'
#' @param copy_number_matrix A copy-number matrix (Matrix) with no columns has a variance of 0.
#'
#' @return An updated copy-number matrix with columns that are all uncorrelated.
#'
#' @author 
#'
#' @importFrom Matrix crossprod colMeans
#' @example
#' cnm <- matrix(c(1, 1, 2, 2, 1,
#'                 1, 1, 2, 2, 1,
#'                 1, 1, 5, 5, 1,
#'                 2, 2, 4, 4, 2),
#'               nrow = 4, byrow = T)
#' colnames(cnm) <- paste0("g_", 1:5)
#' new_cnm <- deduplicate_cnm(cnm, threshold = 0.2)

deduplicate_cnm <- function(copy_number_matrix, threshold = 0.2) {
  csd <- apply(copy_number_matrix, 2, sd)
  if(sum(csd == 0) > 0){
    stop("Columns", which(csd == 0), "have no variance! Cannot compute correlation distances")
  }
  cor_dist <- corrDist(copy_number_matrix)
  index <- which(cor_dist < threshold, arr.ind = TRUE)
  index <- index[index[, 1] != index[, 2], ]
  new_index <- index
  new_index[, 1] <- apply(index, 1, function(x) {min(x)})
  new_index[, 2] <- apply(index, 1, function(x) {max(x)})
  delete_index <- unique(new_index[, 1])
  copy_number_matrix <- copy_number_matrix[, -delete_index]
  return(copy_number_matrix)  
}

#' @name check_corr_in_cnm
#' 
#' @description Check the auto-correlation between columns, and deduplicate the columns correlated with any others from the copy number matrix.
#'
#' @param copy_number_matrix A copy-number matrix (Matrix).
#'
#' @return An updated copy-number matrix with columns that are all uncorrelated.
#'
#' @author 
#'
#' @importFrom Matrix crossprod colMeans
#' @example
#' cnm <- matrix(c(0, 1, 2, 1, 1, 2, 2, 1,
#'                 0, 1, 2, 1, 1, 2, 2, 1,
#'                 0, 1, 2, 1, 1, 5, 5, 1,
#'                 0, 1, 2, 2, 2, 4, 4, 2),
#'               nrow = 4, byrow = T)
#' colnames(cnm) <- paste0("g_", 1:8)
#' new_cnm <- check_corr_in_cnm(cnm, threshold = 0.2)
#' 

check_corr_in_cnm <- function(copy_number_matrix, threshold = 0.2) {
  copy_number_matrix <- as.matrix(copy_number_matrix)
  csum <- colSums(copy_number_matrix)
  
  new_copy_number_matrix <- copy_number_matrix[, which(csum != 0)]
  
  csd <- apply(new_copy_number_matrix, 2, sd)
  
  #If there are no columns with zero variance, 
  #then directly deduplicate the columns correlated with any others from the copy number matrix
  if(sum(csd == 0) == 0) {
    new_copy_number_matrix <- deduplicate_cnm(new_copy_number_matrix, threshold)
  } 
  # #If there is one column with zero variance, then retain the column,
  # #and then deduplicate the columns correlated with any others from the retained copy number matrix
  # else if(sum(csd == 0) == 1) {
  #   csd_0 <- which(csd == 0)
  #   col_csd_0 <- new_copy_number_matrix[, csd_0, drop = F]
  #   new_copy_number_matrix <- new_copy_number_matrix[, -csd_0, drop = F]
  #   new_copy_number_matrix <- deduplicate_cnm(new_copy_number_matrix, threshold)
  #   new_copy_number_matrix <- data.frame(col_csd_0, new_copy_number_matrix)
  # } 
  #If there are one or more columns with zero variance, then retain one of the columns, 
  #and then deduplicate the columns correlated with any others from the retained copy number matrix
  else { #sum(csd == 0) >= 1
    csd_0 <- which(csd == 0)
    col_csd_0 <- new_copy_number_matrix[, csd_0[1], drop = F]
    new_copy_number_matrix <- new_copy_number_matrix[, -csd_0, drop = F]
    new_copy_number_matrix <- deduplicate_cnm(new_copy_number_matrix, threshold)
    new_copy_number_matrix <- data.frame(col_csd_0, new_copy_number_matrix)
  }
  
  cat(ncol(copy_number_matrix) - ncol(new_copy_number_matrix), " genomes have been removed from the copy number matrix.\n")
  removed_genomes <- colnames(copy_number_matrix)[! colnames(copy_number_matrix) %in% colnames(new_copy_number_matrix)]
  cat("The removed genomes: ", removed_genomes, "\n")
  cat(ncol(new_copy_number_matrix), " genomes were retained in the copy number matrix.\n")
  
  return (new_copy_number_matrix)
}

Filter_CNM <- function(cnm, tags_count_file) {
#filter the copynumber matrix according to the vsearch result (delete the tags which are not included in the sample)

  tags_count <- read.table(tags_count_file, sep = "\t", header = T, row.names = 1, comment="")

  #file.remove(tags_count_file) #Delete intermediate result files
  
  idx <- rownames(tags_count) %in% rownames(cnm)
  tags_count <- subset(tags_count, idx)

  idx <- rownames(cnm) %in% rownames(tags_count)
  cnm <- subset(cnm, idx)

  cnm <- cnm[, !duplicated(as.list(as.data.frame(cnm))), drop = F]
  
  # cnm <- unique(cnm, MARGIN=2)
  # cnm <- cnm[, -1] # The first column of the copy number matrix is the tag names

  cnm <- cnm[, which(colSums(cnm) > 0), drop = F]
  
  cnm <- check_corr_in_cnm(cnm)

  tags_count <- tags_count[order(rownames(tags_count)), , drop = FALSE]
  cnm <- cnm[order(rownames(cnm)), , drop = FALSE]
  # When drop = FALSE, the result will maintain its original dimension.

  result <- list(tags_count = tags_count, cnm = cnm)

  return (result)
}


Strain_Level_Profiling <- function(tags_count_matrix, cnm_matrix) {
  
  if(!identical(rownames(tags_count_matrix), rownames(cnm_matrix))) {
    stop("The tags of reads count matrix and that of copy number matrix is not identical")
  }

  predicted_abundance_matrix <- rmscols(tags_count_matrix, cnm_matrix)

  idx <- rowSums(predicted_abundance_matrix) > 0
  predicted_abundance_matrix <- subset(predicted_abundance_matrix, idx)

  predicted_abundance_matrix <- predicted_abundance_matrix / colSums(predicted_abundance_matrix)

  return (predicted_abundance_matrix)
}


Check_file_exists <- function(sample_list_file, species_list_file, cnm_matrix_dir = NULL, cnm_file = NULL, mode = 0) {
  
  if (!file.exists(sample_list_file)) {
    stop("The sequence-files-list file does not exist. Please make sure the file exists and try again.")
  }
  
  sample_list <- read.table(sample_list_file, sep = "\t", header = F)

  lapply(sample_list[, 2], function(x) {
    if (!file.exists(x)) {
      stop(paste0("The sequence file \"", x, "\" does not exist. Please make sure the file exists and try again."))
    }
  })
  
  if(mode == 0) { # mode == 0, profiling based on the 2bGDB
    
     if(!file.exists(cnm_matrix_dir)) {
       stop("The copy number matrix directory does not exist. Please make sure the directory exists and try again.")
     }

  } else { # Profiling based on a customized copy number matrix
    
    if(!file.exists(cnm_file)) {
      stop("The copy number matrix file does not exist. Please make sure the file exists and try again.")
    }
      
  }
  
}


One_Sample_Pipeline <- function(sample_info, species_info, output_path, mode, cnm_matrix_dir, cnm) {
  # mode == 0: profiling based on the 2bGDB; else: profiling based on the customized copy number matrix

  sample_name <- sample_info[1]
  sample_fa <- sample_info[2]

  if(mode == 0) {
    cnm <- Merge_Copynumber_Matrix(species_info, cnm_matrix_dir)
  } else { 
    cnm <- cnm
  }
  
  new_sample_fa <- paste0(output_path, "/new_", sample_name, ".fa")
  Rename_Fasta(sample_name, sample_fa, new_sample_fa)

  similarity <- 1
  output_tag_path <- paste0(output_path, "/", sample_name, ".BcgI.tag")
  tags_count_file <- paste0(output_path, "/", sample_name, "_", similarity, "_tags_count.txt")
  Vsearch(cnm, new_sample_fa, similarity, output_tag_path, tags_count_file)
 
  matrix <- Filter_CNM(cnm, tags_count_file)

  tags_count <- matrix$tags_count
  cnm <- matrix$cnm
  write.table(tags_count, paste0(output_path, "/filtered_", sample_name, "_reads_count_matrix.txt"), sep = "\t", quote = F, row.names = T, col.names = NA)
  write.table(cnm, paste0(output_path, "/filtered_", sample_name, "_copy_number_matrix.txt"), sep = "\t", quote = F, row.names = T, col.names = NA)

  strain_abd_matrix <- Strain_Level_Profiling(tags_count, cnm)

  out_file <- paste0(output_path, "/", sample_name, "_strain_level_abundance.txt")
  write.table(strain_abd_matrix, out_file, sep = "\t", row.names = T, col.names = NA, quote = F)
  
  return (data.frame(sample_name = sample_name, out_file = out_file))
}


Read_Profiling_Matrix <- function(profile_path) {
  result <- read.table(profile, sep = "\t", header = T)
  return (result)
}


Merge_Profiling_Matrix <- function(all_profiles_path, sample_name_list) { 
# merge the strain-level abundance matrix of all samples

  na_count <- sum(!complete.cases(all_profiles_path$out_file))
  
  if(na_count == length(all_profiles_path)) {
    stop("All samples are ineligible for strain-level profiling.")
  }
  
  if(length(all_profiles_path) == 1) {
    result <- Read_Profiling_Matrix(all_profiles_path$out_file[1])
  }
  else {
    # Read all rows where out_file is not NA
    valid_profiles <- all_profiles_path[complete.cases(all_profiles_path$out_file), ]
    # print("....")
    # print(valid_profiles)
    # print("....")
    # Read files corresponding to rows where out_file is not NA, and merge them together
    result <- data.frame()
    for (i in 1:nrow(valid_profiles)) {
      profile_data <- read.table(valid_profiles$out_file[i], header = TRUE, row.names = NULL)
      if (nrow(result) == 0) {
        result <- profile_data
      } else {
        result <- merge(result, profile_data, by = 1, all = TRUE)
        result[is.na(result)] <- 0
      }
    }


    # Check for sample_names that need to be added to the data frame
    missing_samples <- setdiff(sample_name_list, colnames(result))
    # print(missing_samples)
    # Add missing sample_names and set their values to 0
    for (sample in missing_samples) {
      result[, sample] <- 0
    }
  }
  
  rownames(result) <- result[, 1]
  result <- result[, -1]
  return (result)
}


Sample_List_Pipeline <- function(sample_list_file, species_list_file, output_path, cnm_matrix_dir = NULL, cnm_file = NULL, mode = 0) {
  
  if(!file.exists(output_path)) {
    dir.create(output_path)
  }

  Check_file_exists(sample_list_file, species_list_file, cnm_matrix_dir, cnm_file, mode);

  sample_list <- read.table(sample_list_file, sep = "\t", header = F)

  lapply(sample_list[, 1], Rename_sample)

  if (mode == 0){ # mode == 0, profiling based on the 2bGDB
    
    species_list <- read.table(species_list_file, sep = "\t", header = F, comment.char="")
    species_list <- species_list[, 1]
    cnm <- Merge_Copynumber_Matrix(species_list, cnm_matrix_dir)
  
  } else { # Profiling based on a customized copy number matrix
    
    cnm <- read.table(cnm_file, sep = "\t", header = T, row.names = 1)
    
  }

  profile_list <- apply(sample_list, 1, function(x) tryCatch({One_Sample_Pipeline(x, species_list, output_path, mode, cnm_matrix_dir, cnm)}, error=function(err) { NA }))
  
  # print(profile_list)

  profile_list <- do.call(rbind, profile_list)
 
  # print(profile_list)
  write.table(profile_list, paste0(output_path, "/tmp_profile_list.txt"), sep = "\t", quote = F, row.names = T, col.names = NA) 

  sample_name_list <- sample_list[, 1]
  abd_matrix <- Merge_Profiling_Matrix(profile_list, sample_name_list)
  
  write.table(abd_matrix, paste0(output_path, "/strain_level_abd.txt"), sep = "\t", quote = F, row.names = T, col.names = NA)
}



