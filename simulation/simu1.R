library(JUMP)
library(STAREG)
library(adaFilter)
library(JMirror)
library(AdaPTGMM)
library(CARLIS)
library(CoCoNuT)
library(ReAD)
source("./R/SimuFunc.R")
source("./R/methods.R")

nrep = 1  
m  = 10000
pi = c(0.66, 0.01, 0.10, 0.01, 0.10, 0.01, 0.01, 0.10)

A = matrix(c(0.66, 0.01, 0.10, 0.01, 0.10, 0.01, 0.01, 0.10,
             0.10, 0.57, 0.10, 0.01, 0.10, 0.01, 0.01, 0.10,
             0.10, 0.01, 0.66, 0.01, 0.10, 0.01, 0.01, 0.10,
             0.10, 0.01, 0.10, 0.57, 0.10, 0.01, 0.01, 0.10,
             0.10, 0.01, 0.10, 0.01, 0.66, 0.01, 0.01, 0.10,
             0.10, 0.01, 0.10, 0.01, 0.10, 0.57, 0.01, 0.10,
             0.10, 0.01, 0.10, 0.01, 0.10, 0.01, 0.57, 0.10,
             0.04, 0.01, 0.04, 0.01, 0.04, 0.01, 0.01, 0.84),
           8, 8, byrow=TRUE)

mu1 = 2
mu2 = 2
mu3 = 1
alphas = seq(0.01,0.1,0.01)
alphas = c(5e-5, 1e-4, 5e-4, 1e-3, 5e-3, 1e-2, 5e-2, 1e-1)

methods <- c("JUMP", "AdaFilter", "JM", "STAREG", "ReAD", "AdaPTGMM", "CoCoNuT", "CARLIS")
data.obj <- list()

fdp.adafilter <- matrix(nrow = nrep, ncol = length(alphas)); pp.adafilter <- matrix(nrow = nrep, ncol = length(alphas))
fdp.jump <- matrix(nrow = nrep, ncol = length(alphas)); pp.jump <- matrix(nrow = nrep, ncol = length(alphas))
fdp.stareg <- matrix(nrow = nrep, ncol = length(alphas)); pp.stareg <- matrix(nrow = nrep, ncol = length(alphas))
fdp.jm <- matrix(nrow = nrep, ncol = length(alphas)); pp.jm <- matrix(nrow = nrep, ncol = length(alphas))
fdp.read <- matrix(nrow = nrep, ncol = length(alphas)); pp.read <- matrix(nrow = nrep, ncol = length(alphas))
fdp.adaptgmm <- matrix(nrow = nrep, ncol = length(alphas)); pp.adaptgmm <- matrix(nrow = nrep, ncol = length(alphas))
fdp.coconut <- matrix(nrow = nrep, ncol = length(alphas)); pp.coconut <- matrix(nrow = nrep, ncol = length(alphas))
fdp.carlis <- matrix(nrow = nrep, ncol = length(alphas)); pp.carlis <- matrix(nrow = nrep, ncol = length(alphas))
            
for (i in 1:nrep){
  cat(paste0("Replication ", i, "\n"))
  data.obj[[i]] <- SimuData1(m  = m,
                            pi = pi,
                            A  = A,
                            mu1 = mu1,
                            mu2 = mu2,
                            mu3 = mu3)
  
  pa = data.obj[[i]]$pa
  pb = data.obj[[i]]$pb
  x = data.obj[[i]]$x
  theta1 = data.obj[[i]]$theta1
  theta2 = data.obj[[i]]$theta2
  zeta = data.obj[[i]]$zeta
  
  truth <- theta1*theta2
  

  # STAREG
  # tic <- Sys.time()
  res.eb <- stareg(pa, pb)
  # toc <- Sys.time()
  padj.eb <- res.eb$fdr

  # ReAD
  res.hmm <- RepLIS(pa, pb)
  padj.hmm <- res.hmm$fdr

  # AdaPT-GMM
  formulas <- "splines::ns(x, df = 5)"
  res.adaptgmm <- adapt_gmm(x = as.data.frame(x), pvals = pmax(pa, pb), alphas = alphas,
                            beta_formulas = formulas, model_type = "neural", nclasses= c(5),
                            nfits = 5, niter_fit = 3, masking_shape = "comb")

  # AdaFilter
  res.adafilter <- adaFilter(cbind(pa, pb), r = 2, type.I.err = "FDR")
  
  # CoCoNuT
  res.coconut <- coconut(pa, pb, x)

  # CARLIS
  res.carlis <- CARLIS(pa, pb, x)
  padj.carlis <- res.carlis$fdr
  
  j = 0
  for(q in alphas){
    j = j + 1
    # JUMP
    jump.obj <- JUMP(pa, pb, q)
    jump.thr <- jump.obj$jump.thr
    p.max <- jump.obj$p.max
    fdp.jump[i,j] <- sum(p.max <= jump.thr & !truth)/max(sum(p.max <= jump.thr), 1)
    pp.jump[i,j]  <- sum(p.max <= jump.thr & truth)/sum(truth)

    # STAREG
    fdp.stareg[i,j] <- sum(padj.eb <= q & !truth)/max(sum(padj.eb <= q), 1)
    pp.stareg[i,j] <- sum(padj.eb <= q & truth)/sum(truth)

    # ReAD
    fdp.read[i,j] <- sum(padj.hmm <= q & !truth)/max(sum(padj.hmm <= q), 1)
    pp.read[i,j] <- sum(padj.hmm <= q & truth)/sum(truth)

    # AdaPT-GMM
    fdp.adaptgmm[i,j] <- sum(res.adaptgmm$qvals <= q & !truth)/max(sum(res.adaptgmm$qvals <= q), 1)
    pp.adaptgmm[i,j] <- sum(res.adaptgmm$qvals <= q & truth)/sum(truth)
    
    # AdaFilter
    fdp.adafilter[i,j] <- sum(res.adafilter$adjusted.p <= q & !truth)/max(sum(res.adafilter$adjusted.p <= q), 1)
    pp.adafilter[i,j] <- sum(res.adafilter$adjusted.p <= q & truth)/sum(truth)
    
    # Joint mirror
    res.jm <- JointMirror.R(cbind(pa, pb), rank.Mode = "EmptyPoset",
                            trgt.fdr.level = q)
    rej.jm <- rep(0, m)
    rej.jm[res.jm$selected] = 1
    fdp.jm[i,j] <- sum(rej.jm & !truth)/max(sum(rej.jm), 1)
    pp.jm[i,j] <- sum(rej.jm & truth)/sum(truth)

    # coconut
    fdp.coconut[i,j] <- sum(res.coconut$radj <= q & !truth)/max(sum(res.coconut$radj <= q), 1)
    pp.coconut[i,j] <- sum(res.coconut$radj <= q & truth)/sum(truth)

    # carlis
    fdp.carlis[i,j] <- sum(padj.carlis <= q & !truth)/max(sum(padj.carlis <= q), 1)
    pp.carlis[i,j] <- sum(padj.carlis <= q & truth)/sum(truth)
  }
}

