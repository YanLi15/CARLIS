CARLIS <- function(pa, pb, x){
  Rcpp::sourceCpp("src/em_carlis.cpp")
  
  pvals.cutoff = 1e-15
  
  pa[pa == 0] <- min(min(pa[pa != 0]), pvals.cutoff)
  pb[pb == 0] <- min(min(pb[pb != 0]), pvals.cutoff)
  x[x == 0] <- min(min(x[x != 0]), pvals.cutoff)
  
  pi0_pa <- min(pi0est(pa)$pi0, 0.999)
  pi0_pb <- min(pi0est(pb)$pi0, 0.999)
  pi0_x <- min(pi0est(x)$pi0, 0.999)
  
  res <- em_carlis(pa, pb, x, pi0_pa, pi0_pb, pi0_x)
  
  return(res)
}  