p <- c("optparse", "seqinr", "Matrix", "parallel", "glmnet")

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
  # print(1)

  # Initial estimate
  if(verbose) cat("Initial estimate...\n")
  # print(2)

  W <- Diagonal(x = w)
  X <- as.matrix(X) # X must be matrix, can NOT be data.frame
  # print(3)

  theta0 <- Matrix::solve(crossprod(X, crossprod(W, X)), crossprod(X, y))

  theta0 <- pmax(1e-10, theta0)
  names(theta0) <- colnames(X)
  # print(4)

  # Constrained estimate
  ctl <- list(factr = reltol/.Machine$double.eps)
  # print(5)

  if(verbose){
    cat("Constrained optimization...\n")
    ctl <- list(trace = 1, REPORT = 1, factr = reltol/.Machine$double.eps, maxit = 10000)
  }
  
  # print(6)

  set.seed(123)
  # 定义 alpha 值的网格
  alpha_values <- seq(0, 1, by = 0.1)

  # 初始化存储结果的变量
  best_alpha <- 0
  best_lambda <- 0
  best_cv_fit <- NULL
  min_cv_error <- Inf

  # 进行交叉验证以找到最佳的 alpha 和 lambda 值，并设置非负约束
  for (alpha in alpha_values) {
    cv_fit <- cv.glmnet(X, y, alpha = alpha, lower.limits = 0, intercept = FALSE)
    cv_error <- min(cv_fit$cvm)

    if (cv_error < min_cv_error) {
      min_cv_error <- cv_error
      best_alpha <- alpha
      best_lambda <- cv_fit$lambda.min
      best_cv_fit <- cv_fit
    }
  }

  # 输出最佳的 alpha 和 lambda 值
  cat("Best alpha value: ", best_alpha, "\n")
  cat("Best lambda value: ", best_lambda, "\n")

  alpha <- best_alpha
  lambda <- best_lambda
  lst <- optim(theta0, fn = objectFun, gr = NULL, y, X, w, alpha, lambda, method = "L-BFGS-B", lower = rep(0, p), control = ctl)

  #lst <- optim(theta0, fn = objectFun, gr = grr, y, X, w, method = "L-BFGS-B", lower = rep(0, p), control = ctl)
  
  # print(7)

  return(lst$par)
}

objectFun <- function(b, y, X, w, alpha, lambda){

  # Compute the residual sum of squares (RSS)
  rss <- sum((w * (y - X %*% b))^2)
  # Compute the L1 norm (lasso penalty)
  l1_penalty <- sum(abs(b)) 
  # Compute the L2 norm (ridge penalty)
  l2_penalty <- sum(b^2)
  
  # Combine the RSS and penalties
  loss <- rss + lambda * (alpha * l1_penalty + (1 - alpha) / 2 * l2_penalty)
  
  return (loss)
}


#grr <- function(b, y, X, w) { ## Gradient of objectFun

#  r <- t(-2 * w^2 * X) %*% (y - X %*% b) + 0.00002*b

#  return (r)
#}

#objectFun <- function(b, y, X, w){

#  r <- w * (y - (X %*% b))

#  return(sum(r^2) + 0.00001*sum(b^2))
#}

