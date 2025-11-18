 SimuData1 <- function(m   = 10000,
                     pi  = rep(0.125, 8),
                     A   = 0.6 * diag(8) + 0.05,
                     mu1 = 2,
                     mu2 = 2,
                     mu3 = 2){
  h <- c()
  h[1] <- sample(0:7, 1, prob = pi)
  for (j in 2:m){
    h[j] <- sample(0:7, 1, prob = A[h[j-1]+1,])
  }
  
  states1 = rep(0, m)
  states1[c(which(h == 4), which(h == 5), which(h == 6), which(h == 7))] = 1
  states2 = rep(0, m)
  states2[c(which(h == 2), which(h == 3), which(h == 6), which(h == 7))] = 1
  states3 = rep(0, m)
  states3[c(which(h == 1), which(h == 3), which(h == 5), which(h == 7))] = 1
  
  xa <- rnorm(m, mean = mu1 * states1, sd = 1)
  xb <- rnorm(m, mean = mu2 * states2, sd = 1)
  xx <- rnorm(m, mean = mu3 * states2, sd = 1)
  
  pa <- 1 - pnorm(xa)
  pb <- 1 - pnorm(xb)
  x  <- 1 - pnorm(xx)
  
  return(list(
    pa = pa,
    pb = pb,
    x  = x,
    theta1 = states1,
    theta2 = states2,
    zeta = states3))
}


# Four-state Markov chain
SimuData2 <- function(m  = 10000,
                      pi = c(0.25, 0.25, 0.25, 0.25),
                      A = 0.6 * diag(4) + 0.1,
                      info.str = "moderate",
                      mu1 = 2,
                      mu2 = 2,
                      mu3 = 1.5){
  set.seed(NULL)
  
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
  
  truth = states1 * states2
  if(info.str == "uninformative"){
    states3 = sample(0:1, m, replace = TRUE, prob = c(sum(pi[1:2]),sum(pi[3:4])))
  }else if(info.str == "inverse"){
    states3 = as.numeric(!truth)
  }else if(info.str == "same"){
    states3 = truth
  }else if(info.str == "weak"){
    rand.s = sample(0:1, m, replace = TRUE, prob = c(0.3, 0.7))
    states3 = rep(0, m)
    for(j in 1:m){
      if(rand.s[j]==0)
        states3[j] = truth[j]
      else
        states3[j] = sample(0:1, 1, prob = c(sum(pi[1:2]),sum(pi[3:4])))
    }
  }else if(info.str == "fair"){
    rand.s = sample(0:1, m, replace = TRUE, prob = c(0.5, 0.5))
    states3 = rep(0, m)
    for(j in 1:m){
      if(rand.s[j]==0)
        states3[j] = truth[j]
      else
        states3[j] = sample(0:1, 1, prob = c(sum(pi[1:2]),sum(pi[3:4])))
    }
  }else if(info.str == "moderate"){
    rand.s = sample(0:1, m, replace = TRUE, prob = c(0.7, 0.3))
    states3 = rep(0, m)
    for(j in 1:m){
      if(rand.s[j]==0)
        states3[j] = truth[j]
      else
        states3[j] = sample(0:1, 1, prob = c(sum(pi[1:2]),sum(pi[3:4])))
    }
  }else{
    rand.s = sample(0:1, m, replace = TRUE, prob = c(0.9, 0.1))
    states3 = rep(0, m)
    for(j in 1:m){
      if(rand.s[j]==0)
        states3[j] = truth[j]
      else
        states3[j] = sample(0:1, 1, prob = c(sum(pi[1:2]),sum(pi[3:4])))
    }
  }
  stat3 = rnorm(m, states3*mu3, 1)
  p3 = 1 - pnorm(stat3, mean = 0, sd = 1)
  
  
  return(list(
    pa = pa,
    pb = pb,
    x  = p3,
    theta1 = states1,
    theta2 = states2,
    zeta = states3
  ))
}

# Four-state Markov chain
SimuData2_multi <- function(m  = 10000,
                            pi = c(0.25, 0.25, 0.25, 0.25),
                            A = 0.6 * diag(4) + 0.1,
                            info.str = c("moderate", "moderate", "moderate"),
                            mu1 = 2,
                            mu2 = 2,
                            mu3 = c(1.5, 1.5, 1.5)){
  set.seed(NULL)
  
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
  
  truth = states1 * states2
  
  p3 <- list()
  for(i in 1:length(info.str)){
    if(info.str[i] == "uninformative"){
      states3 = sample(0:1, m, replace = TRUE, prob = c(sum(pi[1:2]),sum(pi[3:4])))
    }else if(info.str[i] == "weak"){
      rand.s = sample(0:1, m, replace = TRUE, prob = c(0.3, 0.7))
      states3 = rep(0, m)
      for(j in 1:m){
        if(rand.s[j]==0)
          states3[j] = truth[j]
        else
          states3[j] = sample(0:1, 1, prob = c(sum(pi[1:2]),sum(pi[3:4])))
      }
    }else if(info.str[i] == "fair"){
      rand.s = sample(0:1, m, replace = TRUE, prob = c(0.5, 0.5))
      states3 = rep(0, m)
      for(j in 1:m){
        if(rand.s[j]==0)
          states3[j] = truth[j]
        else
          states3[j] = sample(0:1, 1, prob = c(sum(pi[1:2]),sum(pi[3:4])))
      }
    }else if(info.str[i] == "moderate"){
      rand.s = sample(0:1, m, replace = TRUE, prob = c(0.7, 0.3))
      states3 = rep(0, m)
      for(j in 1:m){
        if(rand.s[j]==0)
          states3[j] = truth[j]
        else
          states3[j] = sample(0:1, 1, prob = c(sum(pi[1:2]),sum(pi[3:4])))
      }
    }else{
      rand.s = sample(0:1, m, replace = TRUE, prob = c(0.9, 0.1))
      states3 = rep(0, m)
      for(j in 1:m){
        if(rand.s[j]==0)
          states3[j] = truth[j]
        else
          states3[j] = sample(0:1, 1, prob = c(sum(pi[1:2]),sum(pi[3:4])))
      }
    }
    stat3 = rnorm(m, states3*mu3[i], 1)
    p3[[i]] = 1 - pnorm(stat3, mean = 0, sd = 1)
  }
  
  return(list(
    pa = pa,
    pb = pb,
    x  = p3,
    theta1 = states1,
    theta2 = states2
  ))
}

