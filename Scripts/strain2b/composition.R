p <- c("optparse", "seqinr", "Matrix")

usePackage <- function(p){
        if (!is.element(p, installed.packages()[,1]))
                install.packages(p, dep=TRUE, repos="http://cran.us.r-project.org/")
        suppressWarnings(suppressMessages(invisible(require(p, character.only=TRUE))))
}
invisible(lapply(p, usePackage))

rmscols <- function(Readcount.mat, Cpn.mat, trim = 0, reltol = 1e-6, verbose = TRUE){
  J <- ncol(Readcount.mat)
  C <- nrow(Cpn.mat)
  G <- ncol(Cpn.mat)
  beta.matrix <- matrix(0, nrow = G, ncol = J)
  rownames(beta.matrix) <- colnames(Cpn.mat)
  colnames(beta.matrix) <- colnames(Readcount.mat)
  if(trim > 0) residuals <- matrix(0, nrow = C, ncol = J)
  for(j in 1:J){
    tot.reads <- sum(Readcount.mat[,j])
    nc <- sum(Readcount.mat[,j] > 0)
    if(verbose) cat("De-convolving sample", colnames(Readcount.mat)[j],
                    "having", tot.reads, "reads mapped to", nc, "clusters...\n")
    if(tot.reads > 0){
      w <- rep(1,C)
      s.hat <- constrLS(Readcount.mat[,j], Cpn.mat, w, reltol, verbose)
      if(trim > 0){
        if(verbose) cat("Trimming residuals...\n")
        residuals[,j] <- as.numeric(Readcount.mat[,j] - Cpn.mat %*% s.hat)
        rsd <- sd(residuals[,j])
        for(g in 1:G){
          idx <- which(Cpn.mat[,g] > 0)
          qq <- quantile(residuals[idx,j], c(trim/2, 1-trim/2))
          idd <- which(residuals[idx,j] < qq[1] | residuals[idx,j] > qq[2])
          w[idx[idd]] <- 0
        }
        s.hat <- constrLS(Readcount.mat[,j], Cpn.mat, w, reltol, verbose)
      }
      beta.hat <- s.hat/sum(s.hat)
      beta.matrix[,j] <- as.numeric(beta.hat)
    }
  }
  return(beta.matrix)
}

# X --- copy number matrix
# y --- reads count matrix
# w --- weight vector, original value = 1
# b --- strain level abundance

### Local functions
constrLS <- function(y, X, w, reltol = 1e-6, verbose = FALSE){
  p <- ncol(X)

  # Initial estimate
  if(verbose) cat(" initial estimate...\n")
  W <- Diagonal(x = w)
  theta0 <- Matrix::solve(crossprod(X, crossprod(W, X)), crossprod(X, y))
  theta0 <- pmax(1e-10, theta0)
  names(theta0) <- colnames(X)
  
  # Constrained estimate
  ctl <- list(factr = reltol/.Machine$double.eps)
  if(verbose){
    cat("   constrained optimization...\n")
    ctl <- list(trace = 1, REPORT = 1, factr = reltol/.Machine$double.eps, maxit = 10000)
  }
  lst <- optim(theta0, fn = objectFun, gr = grr, y, X, w,
               method = "L-BFGS-B", lower = rep(0, p),
               control = ctl)
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

