SimulateData <- function(DOB, P.RSV, n.miss=0) {  #P.RSV is subjects x 365
  age <- 1:NCOL(P.RSV)
  n <- length(DOB)
  Age.RSV <- rep(NA, n)
  Obs.RSV <- rep(0, n)
  ind.obs <- which(runif(n) <= rowSums(P.RSV)); Obs.RSV[ind.obs] <- 1
  Age.RSV[ind.obs] <- apply(X = P.RSV[ind.obs,], MARGIN = 1, function(p){sample(x = age, size = 1, replace = F, prob = p/sum(p))})
  if (n.miss > 0) {
    Age.RSV[ind.obs[sample(x = length(ind.obs), size = n.miss, replace = F)]] <- NA
  }
  return(list(DOB=DOB, Age.RSV=Age.RSV, RSV.1year=Obs.RSV))
}

