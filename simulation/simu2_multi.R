source("./R/SimuFunc.R")
source("./R/methods.R")
source("./R/CARLIS.R")
library(qvalue)

nrep = 100
m  = 50000
info.str = c("fair", "moderate", "strong")
pi1 = 0.8
pipi = c(pi1, (1-pi1-0.05)/2, (1-pi1-0.05)/2, 0.05)
A00 = A11 = A22 = 0.8
A33 = 0.7
mu1 = 2
mu2 = 2
mu3 = c(1.5, 1.5, 1.5)
alphas = 0.05

methods <- c("CARLIS1", "CARLIS2", "CARLIS3", "CARLIS12", "CARLIS13", "CARLIS23", "CARLIS123")
data.obj <- list()

A = matrix(c(A00, (1-A00)/3, (1-A00)/3, (1-A00)/3,
           (1-A11)/3, A11, (1-A11)/3, (1-A11)/3,
           (1-A22)/3, (1-A22)/3, A22, (1-A22)/3,
           (1-A33)/3, (1-A33)/3, (1-A33)/3, A33), 4, 4, byrow=TRUE)



fdp.carlis1 <- matrix(nrow = nrep, ncol = length(alphas)); pp.carlis1 <- matrix(nrow = nrep, ncol = length(alphas))
fdp.carlis2 <- matrix(nrow = nrep, ncol = length(alphas)); pp.carlis2 <- matrix(nrow = nrep, ncol = length(alphas))
fdp.carlis3 <- matrix(nrow = nrep, ncol = length(alphas)); pp.carlis3 <- matrix(nrow = nrep, ncol = length(alphas))
fdp.carlis12 <- matrix(nrow = nrep, ncol = length(alphas)); pp.carlis12 <- matrix(nrow = nrep, ncol = length(alphas))
fdp.carlis13 <- matrix(nrow = nrep, ncol = length(alphas)); pp.carlis13 <- matrix(nrow = nrep, ncol = length(alphas))
fdp.carlis23 <- matrix(nrow = nrep, ncol = length(alphas)); pp.carlis23 <- matrix(nrow = nrep, ncol = length(alphas))
fdp.carlis123 <- matrix(nrow = nrep, ncol = length(alphas)); pp.carlis123 <- matrix(nrow = nrep, ncol = length(alphas))

for (i in 1:nrep){
  cat(paste0("Replication ", i, "\n"))
  data.obj[[i]] <- SimuData2_multi(
    m  = m,
    pi = pipi,
    A = A,
    info.str = info.str,
    mu1 = mu1,
    mu2 = mu2,
    mu3 = mu3
  )
  pa = data.obj[[i]]$pa
  pb = data.obj[[i]]$pb
  x = data.obj[[i]]$x
  theta1 = data.obj[[i]]$theta1
  theta2 = data.obj[[i]]$theta2
  truth <- theta1 * theta2
  
  x1 <- x[[1]]
  x2 <- x[[2]]
  x3 <- x[[3]]
  x12 <- Cauchy(rbind(x1, x2))
  x13 <- Cauchy(rbind(x1, x3))
  x23 <- Cauchy(rbind(x2, x3))
  x123 <- Cauchy(rbind(x1, x2, x3))
  
  # carlis1
  res.carlis1 <- CARLIS(pa, pb, x1)
  padj.carlis1 <- res.carlis1$fdr
  
  # carlis2
  res.carlis2 <- CARLIS(pa, pb, x2)
  padj.carlis2 <- res.carlis2$fdr
  
  # carlis3
  res.carlis3 <- CARLIS(pa, pb, x3)
  padj.carlis3 <- res.carlis3$fdr
  
  # carlis12
  res.carlis12 <- CARLIS(pa, pb, x12)
  padj.carlis12 <- res.carlis12$fdr
  
  # carlis13
  res.carlis13 <- CARLIS(pa, pb, x13)
  padj.carlis13 <- res.carlis13$fdr
  
  # carlis23
  res.carlis23 <- CARLIS(pa, pb, x23)
  padj.carlis23 <- res.carlis23$fdr
  
  # carlis123
  res.carlis123 <- CARLIS(pa, pb, x123)
  padj.carlis123 <- res.carlis123$fdr
  
  j = 0
  for(q in alphas){
    j = j + 1 
    
    fdp.carlis1[i,j] <- sum(padj.carlis1 <= q & !truth)/max(sum(padj.carlis1 <= q), 1)
    pp.carlis1[i,j] <- sum(padj.carlis1 <= q & truth)/sum(truth)
    
    fdp.carlis2[i,j] <- sum(padj.carlis2 <= q & !truth)/max(sum(padj.carlis2 <= q), 1)
    pp.carlis2[i,j] <- sum(padj.carlis2 <= q & truth)/sum(truth)
    
    fdp.carlis3[i,j] <- sum(padj.carlis3 <= q & !truth)/max(sum(padj.carlis3 <= q), 1)
    pp.carlis3[i,j] <- sum(padj.carlis3 <= q & truth)/sum(truth)
    
    fdp.carlis12[i,j] <- sum(padj.carlis12 <= q & !truth)/max(sum(padj.carlis12 <= q), 1)
    pp.carlis12[i,j] <- sum(padj.carlis12 <= q & truth)/sum(truth)
    
    fdp.carlis13[i,j] <- sum(padj.carlis13 <= q & !truth)/max(sum(padj.carlis13 <= q), 1)
    pp.carlis13[i,j] <- sum(padj.carlis13 <= q & truth)/sum(truth)
    
    fdp.carlis23[i,j] <- sum(padj.carlis23 <= q & !truth)/max(sum(padj.carlis23 <= q), 1)
    pp.carlis23[i,j] <- sum(padj.carlis23 <= q & truth)/sum(truth)
    
    fdp.carlis123[i,j] <- sum(padj.carlis123 <= q & !truth)/max(sum(padj.carlis123 <= q), 1)
    pp.carlis123[i,j] <- sum(padj.carlis123 <= q & truth)/sum(truth)
    
  }
}