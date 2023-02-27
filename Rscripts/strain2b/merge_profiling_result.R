Merge_Two_Profiling_Matrix <- function(abd1, abd2) {
        result <- merge(abd1, abd2, by = 1, all = TRUE)
        result[is.na(result)] <- 0
        return (result)
}

Read_Profiling_Matrix <- function(sample_path) {
	result = read.table(sample_path, sep = "\t", header = T)
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

all_samples_info <- read.table("profiling_list.txt", header = F)
#head(all_samples_info)
all_samples_info <- all_samples_info[, "V1"]
#head(all_samples_info)
merged_profiling_matrix <- Merge_Profiling_Matrix(all_samples_info)
rowsum <- rowSums(merged_profiling_matrix)
merged_profiling_matrix <- merged_profiling_matrix[which(rowsum > 0), ]
#merged_profiling_matrix <- t(merged_profiling_matrix)
write.table(merged_profiling_matrix, "profiling_result.txt", quote = F, sep = "\t", row.names = T, col.names = NA)
