Env <-Sys.getenv("Strain2bFunc")
if(nchar(Env)<1){
  cat('Please set the environment variable \"Strain2bFunc\" to the directory\n')
}
print(Env)

source(paste0(Env, "/Scripts/strain2b/composition.R"))


Read_Copynumber_Matrix <- function(species, cnm_matrix_dir) {
  
  if(length(species) == 0) {
    result = NULL
  }
  else {
    cnm_file <- paste0(cnm_matrix_dir, cnms_file, sep = "/")
    
    if (!file.exists(cnm_file)) {
      stop(paste0("The copy number matrix of ", sepcies, " does not exist. Please make sure the file exists and try again."))
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
    
    for (i in 1:nrow(valid_profiles)) {
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
  
  return (result)
}

Rename_sample <- function(sample_name) {
# Since sequence names for vsearch need to include file names, 
# but the file names cannot contain ".", 
# the "." in file names need to be replaced with "_"
  sample_name <- gsub("\\.", "_", sample_name)
  
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
                " --threads 0") # --threads: number of threads to use, zero for all cores (0)
  
  system(cmd, intern = TRUE)
}


Filter_CNM <- function(cnm, tags_count_file) { 
#filter the copynumber matrix according to the vsearch result (delete the tags which are not included in the sample)
  
  tags_count <- read.table(tags_count_file, sep = "\t", header = T, row.names = 1, comment="")
  
  idx <- rownames(tags_count) %in% rownames(cnm)
  tags_count <- subset(tags_count, idx)
  
  idx <- which(rownames(cnm) %in% rownames(tags_count))
  cnm <- subset(cnm, idx)

  cnm <- unique(cnm, MARGIN=2)
  cnm <- cnm[, -1] # The first column of the copy number matrix is the tag names
  
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

  return (result)
}


Check_file_exists <- function(sample_list_file, species_list_file, cnm_matrix_dir = NULL, cnm_file = NULL, mode = 0) {
  
  if (!file.exists(sample_list_file)) {
    stop("The sequence-files-list file does not exist. Please make sure the file exists and try again.")
  }
  
  sample_list <- read.table(sample_list_file, sep = "\t", header = F)
  
  apply(sample_list[, 2], 1, function(x) {
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


One_Sample_Pipeline <- function(sample_info, species_info, output_path, mode, cnm) {
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
  
  matrix <- Filter_CNM(cnm, tags_count_file, sample_name)

  tags_count <- matrix$tags_count
  cnm <- matrix$cnm
  
  strain_abd_matrix <- Strain_Level_Profiling(tags_count, cnm)

  out_file <- paste0(output_path, "/", sample_name, "_strain_level_abundance.txt")
  write.table(strain_abd_matrix, out_file, sep = "\t", row.names = T, col.names = NA, quote = F)
  
  return (data.frame(sample_name = sample_name, out_file = out_file))
}


Read_Profiling_Matrix <- function(profile_path) {
  result <- read.table(profile, sep = "\t", header = T)
  return (result)
}


Merge_Profiling_Matrix <- function(sample_list, all_profiles_path) { 
# merge the strain-level abundance matrix of all samples
  
  na_count <- sum(!complete.cases(all_profiles_path$out_file))
  
  if(na_count == length(all_profiles_path)) {
    stop("All samples are ineligible for strain-level profiling.")
  }
  
  if(length(all_profiles_path) == 1) {
    result <- Read_Profiling_Matrix(all_profiles_path[1])
  }
  else {
    # Read all rows where out_file is not NA
    valid_profiles <- all_profiles_path[complete.cases(all_profiles_path$out_file), ]
    
    # Read files corresponding to rows where out_file is not NA, and merge them together
    result <- data.frame()
    for (i in 1:nrow(valid_profiles)) {
      profile_data <- read.table(valid_profiles$out_file[i], header = TRUE) 
      if (nrow(result) == 0) {
        result <- profile_data
      } else {
        result <- merge(result, profile_data, by = 1, all = TRUE)
        result[is.na(result)] <- 0
      }
    }
    
    # Rename columns to sample_name
    colnames(result)[1] <- "sample_name"
    
    # Check for sample_names that need to be added to the data frame
    missing_samples <- setdiff(sample_list[, 1], result$sample_name)
    
    # Add missing sample_names and set their values to 0
    for (sample in missing_samples) {
      result$sample <- 0
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
  lapply(sample_list, Rename_sample)
  
  if (mode == 0){ # mode == 0, profiling based on the 2bGDB
    
    species_list <- read.table(species_list_file, sep = "\t", header = F, comment.char="")
    cnm <- Merge_Copynumber_Matrix(species, cnm_matrix_dir)
  
  } else { # Profiling based on a customized copy number matrix
    
    cnm <- read.table(cnm_file, sep = "\t", header = T)
    
  }

  profile_list <- apply(sample_list, 1, function(x) tryCatch({One_Sample_Pipeline(x, species, output_path, cnm)}, error=function(err) { NA }))
  
  profile_list <- do.call(rbind, profile_list)
  
  abd_matrix <- Merge_Profiling_Matrix(sample_list, profile_list)
  
  write.table(abd_matrix, paste0(output_path, "/strain_level_abd.txt"), sep = "\t", quote = F, row.names = T, col.names = NA)
}



