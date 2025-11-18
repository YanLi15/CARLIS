# CARLIS
 Covariate-assisted replicability analysis for GWAS

## Installation
```R
## Install dependency packages if necessary
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("qvalue")

## Install CARLIS
install.packages("devtools")
devtools::install_github("YanLi15/CARLIS")

## Load CARLIS
library(CARLIS)
```

## An numeric example
We illustrate the usage of CARLIS for replicability analysis of two large-scale multiple testing problem assisted by two covariates using simulated data under the setting of Simulation II in the manuscript.
```R
## Pre-specify the number of hypotheses, m, the prior probabilities of the joint hidden states of two primary studies, pi's, the transition matrix of the primary studies, the informativeness level of the auxiliary covariate and the alternative settings
m  = 10000;
pi = c(0.8, 0.075, 0.075, 0.05)
A00 = A11 = A22 = 0.7; A33 = 0.8
A = matrix(c(A00, (1-A00)/3, (1-A00)/3, (1-A00)/3,
            (1-A11)/3, A11, (1-A11)/3, (1-A11)/3,
            (1-A22)/3, (1-A22)/3, A22, (1-A22)/3,
            (1-A33)/3, (1-A33)/3, (1-A33)/3, A33), 4, 4, byrow=TRUE)
info.str = c("strong", "strong")
mu1 = 2; mu2 = 2; mu3_1 = 1, mu3_2 = 1

## Generate the hidden states and corresponding p-values in two primary studies 
s <- c()
  s[1] <- sample(0:3, 1, prob = pi)
  for (j in 2:m){
    s[j] <- sample(0:3, 1, prob = A[s[j-1]+1,])
  }
states1 = rep(0, m)
states1[c(which(s == 2), which(s == 3))] = 1
states2 = rep(0, m)
states2[c(which(s == 1), which(s == 3))] = 1

xa <- rnorm(m, mean = mu1 * states1, sd = 1)
xb <- rnorm(m, mean = mu2 * states2, sd = 1)
  
pa <- 1 - pnorm(xa)
pb <- 1 - pnorm(xb)

## Generate the hidden states and two corresponding 'strongly' informative auxiliary covariates
truth = states1 * states2
rand.s_1 = sample(0:1, m, replace = TRUE, prob = c(0.9, 0.1))
states3_1 = rep(0, m)
for(j in 1:m){
   if(rand.s_1[j]==0)
      states3_1[j] = truth[j]
   else
      states3_1[j] = sample(0:1, 1, prob = c(sum(pi[1:2]),sum(pi[3:4])))
}
stat3_1 = rnorm(m, states3_1*mu3_1, 1)
x3_1 = 1 - pnorm(stat3_1, mean = 0, sd = 1)

rand.s_2 = sample(0:1, m, replace = TRUE, prob = c(0.9, 0.1))
states3_2 = rep(0, m)
for(j in 1:m){
   if(rand.s_2[j]==0)
      states3_2[j] = truth[j]
   else
      states3_2[j] = sample(0:1, 1, prob = c(sum(pi[1:2]),sum(pi[3:4])))
}
stat3_2 = rnorm(m, states3_1*mu3_2, 1)
x3_2 = 1 - pnorm(stat3_2, mean = 0, sd = 1)

## Generate the integrated auxiliary covariate via the Cauchy combination rule
Cauchy <- function(Pvals){
  CCT <- colMeans(tan((0.5-Pvals)*pi))
  combined.pval <- 0.5-(atan(CCT))/pi
  return(combined.pval)
}
x <- Cauchy(rbind(x3_1, x3_2))

## Replicability analysis
library(CARLIS)
alpha <- 0.05
rep.obj <- CARLIS(pa, pb, x)
rep.snps <- which(rep.obj$fdr <= alpha)
```

## Data and reproducibility

All the R functions and scripts to reproduce the numeric simulation results in the manuscript are contained in “simulation” folder.

The data used in the real data analysis can be downloaded from the links provided in the *Data availability* section in the manuscript, and the replicability analysis code is provided in the “real_data” folder for illustration of the usage.
