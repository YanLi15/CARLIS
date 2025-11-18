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

nrep = 100
m  = 10000
info.strs = c("uninformative", "fair", "moderate", "strong")
pi1 = 0.8
pi = c(pi1, (1-pi1-0.05)/2, (1-pi1-0.05)/2, 0.05)
A00 = A11 = A22 = 0.8
# A33s = c(0.6, 0.7, 0.8)
A33s = 0.7
mu1s = 2
mu2s = 2
mu3s = 1
q = 0.05

methods <- c("JUMP", "AdaFilter", "JM", "STAREG", "ReAD", "AdaPTGMM", "CoCoNuT", "CARLIS")
results <- c()
data.obj <- list()

for (A33 in A33s){
  A = matrix(c(A00, (1-A00)/3, (1-A00)/3, (1-A00)/3,
               (1-A11)/3, A11, (1-A11)/3, (1-A11)/3,
               (1-A22)/3, (1-A22)/3, A22, (1-A22)/3,
               (1-A33)/3, (1-A33)/3, (1-A33)/3, A33), 4, 4, byrow=TRUE)
  for (mu1 in mu1s){
    for (mu2 in mu2s){
      for (mu3 in mu3s) {
        for (info.str in info.strs){
          fdp.adhocBH <- c(); pp.adhocBH <- c()
          fdp.maxpbh <- c(); pp.maxpbh <- c()
          # fdp.marr <- c(); pp.marr <- c()
          fdp.radjust <- c(); pp.radjust <- c()
          fdp.adafilter <- c(); pp.adafilter <- c()
          fdp.jump <- c(); pp.jump <- c()
          fdp.stareg <- c(); pp.stareg <- c()
          fdp.jm <- c(); pp.jm <- c()
          fdp.read <- c(); pp.read <- c()
          fdp.adaptgmm <- c(); pp.adaptgmm <- c()
          fdp.coconut <- c(); pp.coconut <- c()
          fdp.carlis <- c(); pp.carlis <- c()
          
          for (i in 1:nrep){
            cat(paste0("Replication ", i, "\n"))
            data.obj[[i]] <- SimuData2(m  = m,
                                       pi = pi,
                                       A  = A,
                                       info.str = info.str,
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

            # JUMP
            jump.obj <- JUMP(pa, pb, q)
            jump.thr <- jump.obj$jump.thr
            p.max <- jump.obj$p.max
            fdp.jump <- c(fdp.jump, sum(p.max <= jump.thr & !truth)/max(sum(p.max <= jump.thr), 1))
            pp.jump  <- c(pp.jump, sum(p.max <= jump.thr & truth)/sum(truth))

            # STAREG
            # tic <- Sys.time()
            res.eb <- stareg(pa, pb)
            # toc <- Sys.time()
            padj.eb <- res.eb$fdr
            fdp.stareg <- c(fdp.stareg, sum(padj.eb <= q & !truth)/max(sum(padj.eb <= q), 1))
            pp.stareg <- c(pp.stareg, sum(padj.eb <= q & truth)/sum(truth))

            # ReAD
            res.hmm <- RepLIS(pa, pb)
            padj.hmm <- res.hmm$fdr
            fdp.read <- c(fdp.read, sum(padj.hmm <= q & !truth)/max(sum(padj.hmm <= q), 1))
            pp.read <- c(pp.read, sum(padj.hmm <= q & truth)/sum(truth))
            
            # AdaPT-GMM
            formulas <- "splines::ns(x, df = 5)"
            res.adaptgmm <- adapt_gmm(x = as.data.frame(x), pvals = pmax(pa, pb), alphas = q,
                                      beta_formulas = formulas, model_type = "neural", nclasses= c(5),
                                      nfits = 5, niter_fit = 3, masking_shape = "comb")
            fdp.adaptgmm <- c(fdp.adaptgmm, sum(res.adaptgmm$qvals <= q & !truth)/max(sum(res.adaptgmm$qvals <= q), 1))
            pp.adaptgmm <- c(pp.adaptgmm, sum(res.adaptgmm$qvals <= q & truth)/sum(truth))

            # AdaFilter
            res.adafilter <- adaFilter(cbind(pa, pb), r = 2, type.I.err = "FDR", alpha = q)
            fdp.adafilter <- c(fdp.adafilter, sum(res.adafilter$decision & !truth)/max(sum(res.adafilter$decision), 1))
            pp.adafilter <- c(pp.adafilter, sum(res.adafilter$decision & truth)/sum(truth))

            # Joint mirror
            res.jm <- JointMirror.R(cbind(pa, pb), rank.Mode = "EmptyPoset",
                                    trgt.fdr.level = q)
            rej.jm <- rep(0, m)
            rej.jm[res.jm$selected] = 1
            fdp.jm <- c(fdp.jm, sum(rej.jm & !truth)/max(sum(rej.jm), 1))
            pp.jm <- c(pp.jm, sum(rej.jm & truth)/sum(truth))

            # coconut
            res.coconut <- coconut(pa, pb, x)
            fdp.coconut <- c(fdp.coconut, sum(res.coconut$radj <= q & !truth)/max(sum(res.coconut$radj <= q), 1))
            pp.coconut <- c(pp.coconut, sum(res.coconut$radj <= q & truth)/sum(truth))

            # carlis
            res.carlis <- CARLIS(pa, pb, x)
            padj.carlis <- res.carlis$fdr
            fdp.carlis <- c(fdp.carlis, sum(padj.carlis <= q & !truth)/max(sum(padj.carlis <= q), 1))
            pp.carlis <- c(pp.carlis, sum(padj.carlis <= q & truth)/sum(truth))
          }
          fdr <- c(mean(fdp.adhocBH), mean(fdp.maxpbh), mean(fdp.radjust), mean(fdp.jump), mean(fdp.adafilter), mean(fdp.jm), mean(fdp.stareg), mean(fdp.read), mean(fdp.adaptgmm), mean(fdp.coconut), mean(fdp.carlis))
          power <- c(mean(pp.adhocBH), mean(pp.maxpbh), mean(pp.radjust), mean(pp.jump), mean(pp.adafilter), mean(pp.jm), mean(pp.stareg), mean(pp.read), mean(pp.adaptgmm), mean(pp.coconut), mean(pp.carlis))
          
          results <- rbind(results, data.frame("method" = methods, "mu1" = rep(mu1, 11), "mu2" = rep(mu2, 11), 
                                               "mu3" = rep(mu3, 11), "info.str" = rep(info.str, 11), "A33" = rep(A33, 11), "pi1" = rep(pi1, 11),
                                               "fdr" = fdr, "power" = power))
        }
      }
    }
  }
}  
