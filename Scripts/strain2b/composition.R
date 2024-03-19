p <- c("optparse", "seqinr", "Matrix", "parallel")

usePackage <- function(p){
        if (!is.element(p, installed.packages()[,1]))
                install.packages(p, dep=TRUE, repos="http://cran.us.r-project.org/")
        suppressWarnings(suppressMessages(invisible(require(p, character.only=TRUE))))
}
invisible(lapply(p, usePackage))

rmscols <- function(read_count, copy_number_matrix, trim = 0, reltol = 1e-6, verbose = TRUE){

  C <- nrow(copy_number_matrix)
  G <- ncol(copy_number_matrix)

  result <- matrix(0, nrow = G, ncol = 1)
  rownames(result) <- colnames(copy_number_matrix)
  colnames(result) <- colnames(read_count)
  
  total_reads_count <- sum(read_count[,1])
  nc <- sum(read_count[,1] > 0)
  cat("Sample ", colnames(read_count)[1], " having ", total_reads_count, " reads mapped to ", nc, " clusters\n")

  if(total_reads_count > 0){

    w <- rep(1,C)

    strain_count <- constrLS(read_count[,1], copy_number_matrix, w, reltol, verbose)
  
    if(trim > 0){

      if(verbose) cat("Trimming residuals...\n")

      residuals[,1] <- as.numeric(read_count[,1] - copy_number_matrix %*% strain_count)
    
      rsd <- sd(residuals[,1])
  
      for(g in 1:G){
        idx <- which(copy_number_matrix[,g] > 0)
        qq <- quantile(residuals[idx,j], c(trim/2, 1-trim/2))
        idd <- which(residuals[idx,1] < qq[1] | residuals[idx,1] > qq[2])
        w[idx[idd]] <- 0
      }
  
      strain_count <- constrLS(read_count[,1], copy_number_matrix, w, reltol, verbose)
  
    }
  
    strain_abundance <- strain_count/sum(strain_count)
  
    result[,1] <- as.numeric(strain_abundance)
  }

  return(result)
}

# X --- copy number matrix
# y --- reads count matrix
# w --- weight vector, original value = 1
# b --- strain level abundance

### Local functions
constrLS <- function(y, X, w, reltol = 1e-6, verbose = FALSE){
  p <- ncol(X)
print(1)

  # Initial estimate
  if(verbose) cat("Initial estimate...\n")
print(2)

  W <- Diagonal(x = w)
  X <- as.matrix(X) # X must be matrix, can NOT be data.frame
print(3)
	t <- crossprod(W, X)
	print(class(t))
	print(dim(t))
print(3.1)
	t <- crossprod(X, crossprod(W, X))
	        print(class(t))
	print(dim(t))
print(3.2)
	t <- crossprod(X, y)
        print(class(t))
	print(dim(t))
print(3.3)

  theta0 <- Matrix::solve(crossprod(X, crossprod(W, X)), crossprod(X, y))
print(3.4)

  theta0 <- pmax(1e-10, theta0)
  names(theta0) <- colnames(X)
print(4)

  # Constrained estimate
  ctl <- list(factr = reltol/.Machine$double.eps)
print(5)

  if(verbose){
    cat("Constrained optimization...\n")
    ctl <- list(trace = 1, REPORT = 1, factr = reltol/.Machine$double.eps, maxit = 10000)
  }
print(6)

  lst <- optim(theta0, fn = objectFun, gr = grr, y, X, w, method = "L-BFGS-B", lower = rep(0, p), control = ctl)
print(7)

  return(lst$par)
}

grr <- function(b, y, X, w) { ## Gradient of objectFun

  r <- t(-2 * w^2 * X) %*% (y - X %*% b) + 0.00002*b

  return (r)
}

objectFun <- function(b, y, X, w){

  r <- w * (y - (X %*% b))

  return(sum(r^2) + 0.00001*sum(b^2))
}

