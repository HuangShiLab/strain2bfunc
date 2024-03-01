Env <-Sys.getenv("Strain2bFunc")
if(nchar(Env)<1){
  cat('Please set the environment variable \"Strain2bFunc\" to the directory\n')
 }
print(Env)

source(paste0(Env, "/Scripts/strain2b/composition.R"))

Select_Species_by_Species <- function(Species) {
	return (Species[, 1])
}

Select_Species_by_Abd_Table <- function(species_abd, sample_name, threshold) {
	idx <- which(species_abd[, sample_name] >= threshold)
	result <- rownames(species_abd)[idx]
	#result <- substr(result, 4, nchar(result))
	result <- gsub(" ", "_", result)
	return (result)
}

Merge_Two_Copynumber_Matrix <- function(cnm1, cnm2) { # merge two copynumber matrix
	result <- merge(cnm1, cnm2, by = 1, all = TRUE)
	result[is.na(result)] <- 0
	return (result)
}

Read_Copynumber_Matrix <- function(species) {
	#cnm_matrix_dir <- "/lustre1/g/aos_shihuang/Strain2b/databases/40w_GTDB/CNM/"
	#cnm_matrix_dir = "/lustre1/g/aos_shihuang/Strain2b/databases/CNM_0.01/"
	#cnm_matrix_dir = "/lustre1/g/aos_shihuang/Strain2b/databases/CNM_unique/"
	cnm_matrix_dir = "/lustre1/g/aos_shihuang/Strain2b/databases/CNM_0.001/"
	cnms_file <- list.files(path = cnm_matrix_dir, pattern = paste(species, ".CNM.xls", sep = ""))
	if(length(cnms_file) == 0) {
		result = NULL
	}
	else {
		result = read.table(paste(cnm_matrix_dir, cnms_file, sep = ""), sep = "\t", header = T)
	}
	return (result)
}

Merge_Copynumber_Matrix <- function(all_species) { # merge the copynumber matrix of all species
	if(length(all_species) == 0) {
		result <- NULL
	}
	else if(length(all_species) == 1) {
		result <- Read_Copynumber_Matrix(all_species[1])
	}
	else {
		result <- Read_Copynumber_Matrix(all_species[1])
		for (i in 2:length(all_species)) {
			cnm <-  Read_Copynumber_Matrix(all_species[i])
			result <- Merge_Two_Copynumber_Matrix(result, cnm)
		}
	}
	rownames(result) <- result[, 1]
	#result <- result[, -1]
	return (result)
}

Rename_Fasta <- function(sample_name, sample_fa, new_sample_fa) {
	fasta <- read.fasta(sample_fa, forceDNAtolower = F)	
	sample_name <- gsub("\\.", "_", sample_name)
	names(fasta) <- paste(sample_name, 1:length(fasta), sep = ".")
	write.fasta(sequences = fasta,  names = names(fasta), file.out = new_sample_fa)
}

Vsearch <- function(cnm, new_sample_fa, similarity, output_tag_path, tags_count_file) {
	sink(output_tag_path)
	for (tag in rownames(cnm)) {
		cat(paste0(">", tag, "\n"))
		cat(paste0(tag, "\n"))
	}
	sink()
	cmd <- paste("vsearch --usearch_global ", new_sample_fa, " --db ", output_tag_path, " --id ", similarity, " --iddef 4 --strand both -otutabout " , tags_count_file, " --threads 20", sep = "")
	system(cmd, intern = TRUE)
}

Filter_CNM <- function(cnm, tags_count_file, sample_name) { #filter the copynumber matrix according to the vsearch result (delete the tags which are not included in the sample)
	tags_count <- read.table(tags_count_file, sep = "\t", header = T, row.names = 1, comment="")
	idx <- which(rownames(tags_count) %in% rownames(cnm))
        new_tags_count <- as.data.frame(tags_count[idx, ])
        rownames(new_tags_count) <- rownames(tags_count)[idx]
        colnames(new_tags_count) <- colnames(tags_count)
        tags_count <- as.matrix(new_tags_count)
        idx <- which(rownames(cnm) %in% rownames(tags_count))
	#new_cnm <- as.data.frame(cnm[idx, ])
        #rownames(new_cnm) <- rownames(cnm)[idx]
        #colnames(new_cnm) <- colnames(cnm)
	#cnm <- as.matrix(new_cnm)
        cnm <- cnm[idx, ]  ###
	cnm <- unique(cnm, MARGIN=2)
	cnm <- cnm[, -1]
	print(identical(rownames(tags_count), rownames(cnm)))
	result <- list(tags_count = tags_count, cnm = cnm)
        return (result)
}

Sort_tags <- function(Matrix) {
	idx <- order(rownames(Matrix))
  	new_matrix <- as.data.frame(Matrix[idx, ])
  	rownames(new_matrix) <- rownames(Matrix)[idx]
  	colnames(new_matrix) <- colnames(Matrix)
  	new_matrix <- as.matrix(new_matrix)
	return (new_matrix)
}

Strain_Level_Profiling <- function(tags_count_matrix, cnm_matrix) {
	tags_count_matrix <- Sort_tags(tags_count_matrix)
	cnm_matrix <- Sort_tags(cnm_matrix)

	print(identical(rownames(tags_count_matrix), rownames(cnm_matrix)))

	predicted_abundance_matrix <- rmscols(tags_count_matrix, cnm_matrix)

	idx <- which(predicted_abundance_matrix > 0)
  	rowname <- rownames(predicted_abundance_matrix)
  	colname <- colnames(predicted_abundance_matrix)
  	result <- predicted_abundance_matrix[idx, ]
  	result <- result / sum(result)
  	result <- as.data.frame(result)
  	rownames(result) <- rowname[idx]
  	colnames(result) <- colname
	return (result)
}

One_Sample_Pipeline <- function(sample_info, species_info, output_path, mode, cnm = NULL, threshold = 0.001) {
#for mode == 0, species_info is species abundance table, for mode == 1, species_info is species list
	sample_name <- sample_info[1]
	sample_fa <- sample_info[2]
	#print("1")
	if(mode == 0) {
	  species <- Select_Species_by_Abd_Table(species_info, sample_name, threshold)
	  species <- sub("^s__", "", species)
	  cnm <- Merge_Copynumber_Matrix(species)
	} else { # mode == 1
	  #species <- species_info
	  #cnm <- cnm
	}
 	#print("2")
	output_tag_path <- paste0(output_path, "/", sample_name, ".BcgI.tag")
	new_sample_fa <- paste0(output_path, "/new_", sample_name, ".fa")
	Rename_Fasta(sample_name, sample_fa, new_sample_fa)
  	#print("3")
	similarity <- 1
	tags_count_file <- paste0(output_path, "/", sample_name, "_", similarity, "_tags_count.txt")
	Vsearch(cnm, new_sample_fa, similarity, output_tag_path, tags_count_file)
  	#print("4")
	matrix <- Filter_CNM(cnm, tags_count_file, sample_name)
	#print("5")
	tags_count <- matrix$tags_count
	cnm <- matrix$cnm
	#write.table(cnm, paste0(output_path, "/", sample_name, "_cnm.xls"), sep = "\t", row.names = T, col.names = NA, quote = F)
	predicted_abundance_matrix <- Strain_Level_Profiling(tags_count, cnm)
 	#print("6")
	out_file <- paste0(output_path, "/", sample_name, "_strain_level_profiling.txt")
	write.table(predicted_abundance_matrix, out_file, sep = "\t", row.names = T, col.names = NA, quote = F)
	return (out_file)
}

Merge_Two_Profiling_Matrix <- function(abd1, abd2) {
        result <- merge(abd1, abd2, by = 1, all = TRUE)
        result[is.na(result)] <- 0
        return (result)
}

Read_Profiling_Matrix <- function(sample_path) {
        result <- read.table(sample_path, sep = "\t", header = T)
        return (result)
}

Merge_Profiling_Matrix <- function(all_samples_info) { # merge the copynumber matrix of all samples
        if(length(all_samples_info) == 0) {
                result <- NULL
        }
        else if(length(all_samples_info) == 1) {
                result <- Read_Profiling_Matrix(all_samples_info[1])
        }
        else {
                result <- Read_Profiling_Matrix(all_samples_info[1])
                for (i in 2:length(all_samples_info)) {
                        abd <-  Read_Profiling_Matrix(all_samples_info[i])
                        result <- Merge_Two_Profiling_Matrix(result, abd)
                }
        }
        rownames(result) <- result[, 1]
        result <- result[, -1]
        return (result)
}

Sample_List_Pipeline <- function(sample_list_file, species_file, output_path, mode = 0, threshold = 0.001) {
  	if(!file.exists(output_path)) {
 	   	dir.create(output_path)
 	 }
 	sample_list <- read.table(sample_list_file, sep = "\t", header = F)
  	sample_list[, 1] <- gsub("-", ".", sample_list[, 1])
  
	if(mode == 0) { #mode == 0
	  species_abd <- read.table(species_file, sep = "\t", header = T, comment.char="")
	  if(ncol(species_abd) < 8) {
	  	stop("Please confirm that there is at least on sample in the species abundance table file.")
	  }
	  else if(ncol(species_abd) == 8) {
	  	rowname <- species_abd$Species
	  	colname <- gsub("-", ".", colnames(species_abd))[8:ncol(species_abd)]
	  	species_abd <- as.data.frame(species_abd[, 8:ncol(species_abd)]) #the first 7 columns of the result of 2bRAD-M is classification information
	  	rownames(species_abd) <- rowname
	  	colnames(species_abd) <- colname
	  }
	  else {
		rownames(species_abd) <- species_abd$Species
	  	species_abd <- species_abd[, 8:ncol(species_abd)] #the first 7 columns of the result of 2bRAD-M is classification information
		colnames(species_abd) <- gsub("-", ".", colnames(species_abd))
		sample_list <- sample_list[order(sample_list[, 1]), ]
		species_abd <- species_abd[, order(colnames(species_abd))]
	  }
	  sample_names0 <- sample_list[, 1]
          sample_names <- colnames(species_abd)
	  if(!identical(sample_names0, sample_names)) {
	    stop("Please confirm that the sample names in the first column of the sample list file match 
	         those in the first row of the species abundance table file.")
	  }
	  profile_list <- apply(sample_list, 1, function(x) tryCatch({One_Sample_Pipeline(x, species_abd, output_path, mode, threshold)}, error=function(err) { NA }))
	  #profile_list <- mclapply(sample_list, 1, function(x) tryCatch({One_Sample_Pipeline(x, species_abd, output_path, mode, threshold)}, error=function(err) { NA }), mc.cores = 4)
	  #profile_list <- mclapply(sample_list, 1, function(x) One_Sample_Pipeline(x, species_abd, output_path, mode, threshold), mc.cores = 4)
	} else { #mode == 1
	  species_list <- read.table(species_file, sep = "\t", header = F, comment.char="")
	  species <- Select_Species_by_Species(species_list)
	  cnm <- Merge_Copynumber_Matrix(species)
	  profile_list <- apply(sample_list, 1, function(x) tryCatch({One_Sample_Pipeline(x, species, output_path, mode, cnm)}, error=function(err) { NA }))
	  #profile_list <- mclapply(sample_list, 1, function(x) tryCatch({One_Sample_Pipeline(x, species, output_path, mode, cnm)}, error=function(err) { NA }), mc.cores = 4)
	  #profile_list <- mclapply(sample_list, 1, function(x) One_Sample_Pipeline(x, species, output_path, mode, cnm), mc.cores = 4)
	}

	abd_matrix <- Merge_Profiling_Matrix(profile_list)
	write.table(abd_matrix, paste0(output_path, "/strain_level_abd.txt"), sep = "\t", quote = F, row.names = T, col.names = NA)
}



