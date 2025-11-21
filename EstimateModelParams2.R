Get.Data <- function(degree=3, n.internal.knots=10, knots=NULL, Lambda, Date.RSV, RSV.status, DOB, Lambda.null=F, equidistant.knots=F, Cov=NULL, include.int=F) {
  if (!is.null(Cov)) {
    ind.use <- which(!is.na(rowSums(Cov)) & !is.na(RSV.status))
    if (include.int) {
      Cov <- cbind(rep(1,length(ind.use)),apply(X = Cov, MARGIN = 2, function(x){x[ind.use]-mean(x[ind.use])}))
    } else {
      Cov <- cbind(apply(X = Cov, MARGIN = 2, function(x){x[ind.use]-mean(x[ind.use])}))
    }
    Date.RSV <- Date.RSV[ind.use]; RSV.status <- RSV.status[ind.use]; DOB <- DOB[ind.use]; rm(ind.use)
  }
  if (Lambda.null) {Lambda$lambda <- rep(1,length(Lambda$lambda))}
  ind.use <- which(!is.na(RSV.status)); Date.RSV <- Date.RSV[ind.use]; RSV.status <- RSV.status[ind.use]; DOB <- DOB[ind.use]
  out <- list()
  n <- length(DOB)
  DOB <- as.Date(DOB); Date.RSV <- as.Date(Date.RSV); Lambda$date <- as.Date(Lambda$date)
  a.max <- 365
  out$Ind.obs <- which(!is.na(Date.RSV)); out$Ind.NoRSV <- which(RSV.status==0); out$Ind.Miss <- which(RSV.status==1 & is.na(Date.RSV))
  out$Age.Infection <- as.numeric(Date.RSV[out$Ind.obs] - DOB[out$Ind.obs])
  out$DOB <- DOB
  if (is.null(knots)) {
    if (!equidistant.knots) {
      knots <- c(rep(0,degree+1),quantile(x = out$Age.Infection, (1:n.internal.knots)/(n.internal.knots+1)),rep(a.max,degree+1))
    } else {
      knots <- c(rep(0,degree+1),seq(0,a.max,length=n.internal.knots+2)[-c(1,n.internal.knots+2)],rep(a.max,degree+1))
    }
  }
  out$B.obs <- splines::splineDesign(knots = knots, x = out$Age.Infection, ord = degree+1)
  out$knots <- knots
  out$degree <- degree
  
  Integrals <- sapply(1:ncol(out$B.obs), function(j){  #A x J 
    ind.j <- rep(0,length(knots)-degree-1)
    ind.j[j] <- 1
    sapply(1:a.max,function(x){
      max(0,IntegrateBs::ibs(x=x, knots = c(rep(0,degree+1),knots), ord = degree+1, coef = c(rep(0,degree+1),ind.j)) - IntegrateBs::ibs(x=x-1, knots = c(rep(0,degree+1),knots), ord = degree+1, coef = c(rep(0,degree+1),ind.j)))
    })
  })
  
  t.0 <- min(DOB)
  Lambda.func <- rep(NA,as.numeric(max(DOB)-t.0)+1+a.max)
  for (i in 1:length(Lambda.func)) {
    date.i <- t.0 + i - 1
    time.i <- as.numeric(date.i - Lambda$date)
    ind.i.max <- max(which(time.i>=0))
    if (time.i[ind.i.max]==0) {Lambda.func[i] <- Lambda$lambda[ind.i.max]}
    lambda.1 <- Lambda$lambda[ind.i.max]; lambda.2 <- Lambda$lambda[ind.i.max+1]
    Lambda.func[i] <- lambda.1 + (lambda.2 - lambda.1)/(as.numeric(Lambda$date[ind.i.max+1]-Lambda$date[ind.i.max]))*time.i[ind.i.max]
  }
  Tmp.mat <- sapply(1:n, function(i) {  #A x n
    tmp.i <- as.numeric(DOB[i]-t.0)+1
    if (i %in% out$Ind.obs) {
      L <- rep(0,a.max)
      ind.i <- (tmp.i+1):(tmp.i+out$Age.Infection[which(out$Ind.obs==i)])
      L[1:length(ind.i)] <- Lambda.func[ind.i]
      return(L)
    }
    return(Lambda.func[(tmp.i+1):(tmp.i+a.max)])
  })
  out$B.tilde <- t(Tmp.mat)%*%Integrals
  out$Cov <- Cov
  return(out)
}

Optimize.W <- function(Data, lambda=0, Alpha.start=NULL, max.iter=1e4, tol=1e-6) {
  B.obs <- Data$B.obs
  colSum.B.tilde.obs <- colSums(Data$B.tilde[Data$Ind.obs,])
  colSum.B.tilde.noRSV <- colSums(Data$B.tilde[Data$Ind.NoRSV,])
  B.tilde.miss <- Data$B.tilde[Data$Ind.Miss,]
  J <- NCOL(B.obs)
  if (lambda>0) {
    M <- t(sapply(1:(J-2),function(j){out <- rep(0,J); out[j] <- 1; out[j+1] <- -2; out[j+2] <- 1; return(out)}))
    M <- t(M)%*%M
  } else {
    M <- NULL
  }
  
  #Starting point#
  if (is.null(Alpha.start)) {
    tmp <- optimise(f = Compute.Likelihood, interval = c(0.01,10), maximum = T, B.obs=B.obs, colSum.B.tilde.obs=colSum.B.tilde.obs, B.tilde.miss=B.tilde.miss, colSum.B.tilde.noRSV=colSum.B.tilde.noRSV, lambda=lambda, M=M)
    Alpha <- rep(tmp$maximum,J)
  } else {
    Alpha <- Alpha.start
  }
  
  #Run algorithm#
  for (i in 1:max.iter) {
    G <- Compute.Grad(Alpha = Alpha, B.obs = B.obs, colSum.B.tilde.obs = colSum.B.tilde.obs, B.tilde.miss = B.tilde.miss, colSum.B.tilde.noRSV = colSum.B.tilde.noRSV, lambda = lambda, M = M)
    if (Check.KKT(Alpha = Alpha, G = G, tol = tol)) {return(list(Alpha=Alpha,n.iter=i,out=0))}
    H <- Compute.Hess(Alpha = Alpha, B.obs = B.obs, B.tilde.miss = B.tilde.miss, lambda = lambda, M = M)
    dir.i <- Sub.problem(Alpha.0 = Alpha, G = G, H = H, max.iter = max.iter, tol = tol)$Alpha - Alpha
    cons.i <- optimise(f = Compute.Likelihood, interval = c(0,1), maximum = T, B.obs=B.obs, colSum.B.tilde.obs=colSum.B.tilde.obs, B.tilde.miss=B.tilde.miss, colSum.B.tilde.noRSV=colSum.B.tilde.noRSV, lambda=lambda, M=M, Alpha.0=Alpha, dir=dir.i)
    Alpha <- Alpha + cons.i$maximum*dir.i
  }
  return(list(Alpha=Alpha,n.iter=i,out=1))
}

Compute.Likelihood <- function(Alpha, B.obs, colSum.B.tilde.obs, B.tilde.miss, colSum.B.tilde.noRSV, lambda=0, M=NULL, Alpha.0=NULL, dir=NULL) {
  if (length(Alpha)==1) {
    if (is.null(dir)) {
      Alpha <- rep(Alpha,NCOL(B.obs))
    } else {
      Alpha <- Alpha*dir + Alpha.0
    }
  }
  out <- sum(log(c(B.obs%*%Alpha))) - sum(Alpha*colSum.B.tilde.obs) - sum(Alpha*colSum.B.tilde.noRSV)
  out <- out + sum(log(1 - exp(-c(B.tilde.miss%*%Alpha))))
  if (is.null(M)) {return(out)}
  return(out - lambda/2*sum(Alpha*c(M%*%Alpha)))
}

Compute.Grad <- function(Alpha, B.obs, colSum.B.tilde.obs, B.tilde.miss, colSum.B.tilde.noRSV, lambda=0, M=NULL) {
  out <- c(t(B.obs)%*%(1/c(B.obs%*%Alpha))) - colSum.B.tilde.obs - colSum.B.tilde.noRSV
  tmp <- 1/( exp(c(B.tilde.miss%*%Alpha)) - 1 )
  out <- out + c( t(B.tilde.miss)%*%tmp )
  if (is.null(M)) {return(out)}
  return(out - lambda*c(M%*%Alpha))
}

Compute.Hess <- function(Alpha, B.obs, B.tilde.miss, lambda=0, M=NULL) {
  if (lambda>0 & is.null(M)) {
    J <- length(Alpha)
    M <- t(sapply(1:(J-2),function(j){out <- rep(0,J); out[j] <- 1; out[j+1] <- -2; out[j+2] <- 1; return(out)}))
    M <- t(M)%*%M
  }
  tmp <- 1/c(B.obs%*%Alpha)^2
  out <- -t(B.obs*tmp)%*%B.obs
  tmp <- exp(c(B.tilde.miss%*%Alpha)) / (exp(c(B.tilde.miss%*%Alpha)) - 1)^2
  out <- out - t(B.tilde.miss*tmp)%*%B.tilde.miss
  if (is.null(M)) {return(out)}
  return(out - lambda*M)
}

Sub.problem <- function(Alpha.0, G, H, max.iter=1e4, tol=1e-6) {
  Alpha <- Alpha.0 - c(solve(H,G))
  if (all(Alpha>=0)) {return(list(Alpha=Alpha,n.iter=0,out=0))}
  if (all(Alpha<=0)) {
    Alpha <- Alpha.0
  } else {
    Alpha[Alpha<=0] <- min(Alpha[Alpha>0])
  }
  K <- length(Alpha)
  for (i in 1:max.iter) {
    for (j in 1:K) {
      x.mj <- Alpha; x.mj[j] <- 0
      Alpha[j] <- max(0, -(G[j] + sum(H[,j]*(x.mj - Alpha.0)))/H[j,j])
    }
    if (Check.KKT(Alpha = Alpha, G = c(H%*%(Alpha-Alpha.0)) + G)) {return(list(Alpha=Alpha,n.iter=i,out=0))}
  }
  return(list(Alpha=Alpha,n.iter=i,out=1))
}

Check.KKT <- function(Alpha, G, tol=1e-6) {
  if (all(abs(G)<tol)) {return(T)}
  if (all( abs(G)<tol | (G<=0 & Alpha==0) )) {return(T)}
  return(F)
}

Choose.Lambda <- function(Data, n.folds=5, n.lambda=100, max.lambda=10) {
  Fold.obs <- Get.folds(n = length(Data$Ind.obs), n.folds = n.folds)
  Fold.miss <- Get.folds(n = length(Data$Ind.Miss), n.folds = n.folds)
  Fold.NoRSV <- Get.folds(n = length(Data$Ind.NoRSV), n.folds = n.folds)
  lambda <- c(0,exp(seq(log(0.1),log(max.lambda),length=n.lambda)))
  Like.out <- matrix(NA,nrow=length(lambda),ncol=n.folds)
  n <- nrow(Data$B.tilde)

  for (k in 1:n.folds) {
    n.k <- 0
    Data.k <- Data
    ind.k <- which(Fold.obs==k); n.k <- n.k + length(ind.k)
    Data.k$Ind.obs <- Data.k$Ind.obs[-ind.k]; Data.k$B.obs <- Data.k$B.obs[-ind.k,]
    ind.k <- which(Fold.miss==k); Data.k$Ind.Miss <- Data.k$Ind.Miss[-ind.k]; n.k <- n.k + length(ind.k)
    ind.k <- which(Fold.NoRSV==k); Data.k$Ind.NoRSV <- Data.k$Ind.NoRSV[-ind.k]; n.k <- n.k + length(ind.k)
    Alpha.start <- NULL
    for (i in 1:length(lambda)) {
      out.i <- Optimize.W(Data = Data.k, lambda = (n-n.k)*lambda[i], max.iter = 1e3, tol = 1e-6, Alpha.start = Alpha.start)
      Alpha.start <- out.i$Alpha
      Like.out[i,k] <- 1/n.k*Compute.Likelihood(Alpha = out.i$Alpha, B.obs = Data$B.obs[which(Fold.obs==k),], colSum.B.tilde.obs = colSums(rbind(Data$B.tilde[Data$Ind.obs[which(Fold.obs==k)],])), B.tilde.miss = rbind(Data$B.tilde[Data$Ind.Miss[which(Fold.miss==k)],]), colSum.B.tilde.noRSV = colSums(rbind(Data$B.tilde[Data$Ind.NoRSV[which(Fold.NoRSV==k)],])), M = NULL)   
    }
  }
  lambda.opt <- lambda[which.max(rowSums(Like.out))]
  out <- Optimize.W(Data = Data, lambda = n*lambda.opt, max.iter = 1e3, tol = 1e-6)
  return(list(lambda=lambda,Like.all=Like.out,Like=rowSums(Like.out),lambda.opt=lambda.opt,Alpha.opt=out$Alpha))
}

Plot.Weight <- function(Alpha, Data, units=c("days","weeks","months"), add=F, col.add="blue", ...) {
  units <- match.arg(units, c("days","weeks","months"))
  x <- seq(0, 365, by=0.1)
  y <- as.numeric(c(splines::splineDesign(knots = Data$knots, x = x, ord = Data$degree+1)%*%Alpha))
  if (units=="days") {
    if (!add) {
      plot(x, y, xlab="Age (days)", ylab="Weight (norm. to have max 1)", type="l")
    } else {
      lines(x, y, col=col.add, ...)
    }
  }
  if (units=="weeks") {
    if (!add) {
      plot(x/7, y*7, xlab="Age (weeks)", ylab="Weight (norm. to have max 1)", type="l")
    } else {
      lines(x/7, y*7, col=col.add, ...)
    }
  }
  if (units=="months") {
    if (!add) {
      plot(x/30, y*30, xlab="Age (months)", ylab="Weight (norm. to have max 1)", type="l")
    } else {
      lines(x/30, y*30, col=col.add, ...)
    }
  }
}

Get.folds <- function(n, n.folds) {
  out <- rep(1:n.folds, each=ceiling(n/n.folds))[1:n]
  return(out[order(runif(n))])
}

Var.Pred <- function(x=NULL, Alpha, Data, lambda=0, M=NULL, units=c("days","weeks","months")) {
  units <- match.arg(units, c("days","weeks","months"))
  Var.Alpha <- solve(-Compute.Hess(Alpha=Alpha, B.obs=Data$B.obs, B.tilde.miss=Data$B.tilde[Data$Ind.Miss,], lambda=lambda, M=M))
  if (is.null(x)) {x <- seq(0, 365, by=0.1)}
  X <- splines::splineDesign(knots = Data$knots, x = x, ord = Data$degree+1)
  W.hat <- as.numeric(c(splines::splineDesign(knots = Data$knots, x = x, ord = Data$degree+1)%*%Alpha))
  s.What <- sqrt(rowSums((X%*%Var.Alpha)*X)); rm(X)
  x.plot <- x/ifelse(units=="days",1,ifelse(units=="weeks",7,30))
  plot(x.plot, W.hat/max(W.hat), xlab=paste0("Age (",units,")"), ylab="Weight (normalized to have max 1)", type="l", ylim=c(0,max(1.96*s.What/max(W.hat)+W.hat/max(W.hat))))
  for (i in 1:length(x.plot)) {
    lines(rep(x.plot[i],100), seq(-1.96*s.What[i]/max(W.hat)+W.hat[i]/max(W.hat),1.96*s.What[i]/max(W.hat)+W.hat[i]/max(W.hat),length=100), col="grey")
  }
  lines(x.plot, W.hat/max(W.hat), lwd=2)
}

Var.Pred.Cov <- function(x=NULL, Alpha, Beta=NULL, Data, lambda=0, units=c("days","weeks","months"), plot.it=T) {
  units <- match.arg(units, c("days","weeks","months"))
  tmp <- solve(-Hess.full(Alpha = Alpha, Beta = Beta, Data = Data, lambda = lambda))
  Var.Total <- tmp%*%(-Hess.full(Alpha = Alpha, Beta = Beta, Data = Data, lambda = 0))%*%tmp
  Eff.Params <- -sum(diag(Hess.full(Alpha = Alpha, Beta = Beta, Data = Data, lambda = 0)  %*% tmp))
  Var.Alpha <- Var.Total[1:length(Alpha),1:length(Alpha)]
  if (!is.null(Beta)) {
    Var.Beta <- Var.Total[(length(Alpha)+1):(length(Alpha)+length(Beta)),(length(Alpha)+1):(length(Alpha)+length(Beta))]
  } else {
    Var.Beta <- NULL
  }
  if (is.null(x)) {x <- seq(0, 365, by=0.1)}
  X <- splines::splineDesign(knots = Data$knots, x = x, ord = Data$degree+1)
  k <- 1
  if (units=="weeks") {k <- 1/7}
  if (units=="months") {k <- 1/30}
  W.hat <- as.numeric(c(splines::splineDesign(knots = Data$knots, x = x, ord = Data$degree+1)%*%Alpha))/k
  s.What <- sqrt(rowSums((X%*%Var.Alpha)*X))/k; rm(X)
  x.plot <- x*k
  if (plot.it) {
    plot(x.plot, W.hat/max(W.hat), xlab=paste0("Age (",units,")"), ylab="Weight (normalized to have max 1)", type="l", ylim=c(0,max(1.96*s.What/max(W.hat)+W.hat/max(W.hat))))
    for (i in 1:length(x.plot)) {
      lines(rep(x.plot[i],100), seq(-1.96*s.What[i]/max(W.hat)+W.hat[i]/max(W.hat),1.96*s.What[i]/max(W.hat)+W.hat[i]/max(W.hat),length=100), col="grey")
    }
    lines(x.plot, W.hat/max(W.hat), lwd=2)
  }
  if (!is.null(Beta)) {return(list( Eff.Params=Eff.Params, x.plot=x.plot, W.hat=W.hat, se.W.hat=s.What, Beta=Beta, se.Beta=sqrt(diag(Var.Beta)), Var.Beta=Var.Beta ))}
  return(list( Eff.Params=Eff.Params, x.plot=x.plot, W.hat=W.hat, se.W.hat=s.What ))
}


#####With random effect#####
Update.phi <- function(Alpha, B.obs, B.tilde.obs, B.tilde.miss, B.tilde.NoRSV, phi.try=c(0,seq(0.1,10,by=0.1))^2, plot.it=F) {
  ind.0 <- which(phi.try==0)
  if (length(ind.0)) {
    ll.0 <- Compute.Likelihood(Alpha = Alpha, B.obs = B.obs, colSum.B.tilde.obs = colSums(B.tilde.obs), B.tilde.miss = B.tilde.miss, colSum.B.tilde.noRSV = colSums(B.tilde.NoRSV), lambda = 0, M = NULL)
    phi.try <- phi.try[-ind.0]
  } else {
    ll.0 <- -Inf
  }
  a.try <- 1/phi.try
  tmp <- c(B.tilde.obs%*%Alpha)
  ll1 <- length(tmp)*(a.try+1)*log(a.try) - (a.try+1)*sapply(a.try,function(a){sum(log(a+tmp))}) + sum(log(c(B.obs%*%Alpha)))
  tmp <- c(B.tilde.NoRSV%*%Alpha)
  ll2 <- length(tmp)*a.try*log(a.try) - a.try*sapply(a.try,function(a){sum(log(a+tmp))})
  tmp <- c(B.tilde.miss%*%Alpha)
  ll3 <- sapply(a.try,function(a){sum(log( 1 - 1/(1+tmp/a)^a ))})
  if (plot.it) {
    if (length(ind.0)) {
      plot(c(0,sqrt(phi.try)), c(ll.0,ll1+ll2+ll3), xlab=expression(sqrt(phi)), ylab="log-likelihood")
    } else {
      plot(sqrt(phi.try), ll1+ll2+ll3, xlab=expression(sqrt(phi)), ylab="log-likelihood")
    }
  }
  return(c(0,phi.try)[which.max(c(ll.0,ll1+ll2+ll3))])
}

####With additional covariates####
Optimize.W.Cov <- function(Data, lambda=0, Alpha.start=NULL, beta.start=NULL, max.iter=1e4, tol=1e-6, tol.beta=1e-6) {
  B.tilde.obs <- Data$B.tilde[Data$Ind.obs,]
  B.tilde.miss <- Data$B.tilde[Data$Ind.Miss,]
  B.tilde.NoRSV <- Data$B.tilde[Data$Ind.NoRSV,]
  X.obs <- Data$Cov[Data$Ind.obs,]; X.miss <- Data$Cov[Data$Ind.Miss,]; X.NoRSV <- Data$Cov[Data$Ind.NoRSV,]
  
  if (is.null(Alpha.start) | is.null(beta.start)) {
    Alpha <- Optimize.W(Data = Data, lambda = lambda, Alpha.start = NULL, max.iter = max.iter, tol = tol)$Alpha
    beta <- c(1,rep(0,NCOL(X.obs)-1))
  } else {
    Alpha <- Alpha.start
    beta <- beta.start
  }
  
  for (i in 1:max.iter) {
    beta.new <- Update.beta(beta.start = beta, X.obs = X.obs, X.NoRSV = X.NoRSV, X.miss = X.miss, aBtilde.obs = c(B.tilde.obs%*%Alpha), aBtilde.NoRSV = c(B.tilde.NoRSV%*%Alpha), aBtilde.miss = c(B.tilde.miss%*%Alpha))
    delta <- max(abs(beta.new-beta))
    if (delta<tol.beta) {
      beta <- beta.new
      return(list(Alpha=Alpha,beta=beta,n.iter=i,out=0,delta=delta))
    }
    beta <- beta.new
    Data.i <- Data; Data.i$B.tilde <- exp(c(Data.i$Cov%*%beta))*Data.i$B.tilde
    Alpha <- Optimize.W(Data = Data.i, lambda = lambda, Alpha.start = Alpha, max.iter = max.iter, tol = tol)$Alpha
  }
  return(list(Alpha=Alpha,beta=beta,n.iter=i,out=1,delta=delta))
}

LogLike.beta <- function(beta, X.obs, X.NoRSV, X.miss, aBtilde.obs, aBtilde.NoRSV, aBtilde.miss) {
  X.obs <- cbind(X.obs); X.NoRSV <- cbind(X.NoRSV); X.miss <- cbind(X.miss)
  tmp <- c(X.obs%*%beta)
  l1 <- sum(tmp - exp(tmp)*aBtilde.obs)
  l2 <- -sum(exp(c(X.NoRSV%*%beta))*aBtilde.NoRSV)
  l3 <- sum(log( 1 - exp( -exp(c(X.miss%*%beta))*aBtilde.miss ) ))
  return(l1+l2+l3)
}

Grad.LogLike.beta <- function(beta, X.obs, X.NoRSV, X.miss, aBtilde.obs, aBtilde.NoRSV, aBtilde.miss) {
  X.obs <- cbind(X.obs); X.NoRSV <- cbind(X.NoRSV); X.miss <- cbind(X.miss)
  l1 <- t(X.obs)%*%( 1 - exp(c(X.obs%*%beta))*aBtilde.obs )
  l2 <- -t(X.NoRSV)%*%( exp(c(X.NoRSV%*%beta))*aBtilde.NoRSV )
  tmp <- c(X.miss%*%beta)
  l3 <- t(X.miss)%*%( exp(tmp)*aBtilde.miss/( exp(exp(tmp)*aBtilde.miss) - 1 ) )
  return(c(l1+l2+l3))
}

Update.beta <- function(beta.start=NULL, X.obs, X.NoRSV, X.miss, aBtilde.obs, aBtilde.NoRSV, aBtilde.miss) {
  if (is.null(beta.start)) {beta.start <- c(1,rep(0,NCOL(X.obs)-1))}
  out <- optim(par = beta.start, fn = LogLike.beta, gr = Grad.LogLike.beta, method = "BFGS", control = list(fnscale=-1), X.obs=X.obs, X.NoRSV=X.NoRSV, X.miss=X.miss, aBtilde.obs=aBtilde.obs, aBtilde.NoRSV=aBtilde.NoRSV, aBtilde.miss=aBtilde.miss)
  return(out$par)
}

Hess.full <- function(Alpha, Beta=NULL, Data, lambda=0) {
  Cov <- Data$Cov
  if (is.null(Cov) | is.null(Beta)) {
    return(Compute.Hess(Alpha = Alpha, B.obs = Data$B.obs, B.tilde.miss = Data$B.tilde[Data$Ind.Miss,], lambda = lambda))
  }
  J <- length(Alpha)
  Data$B.tilde <- exp(c(Cov%*%Beta))*Data$B.tilde
  p <- length(Beta)
  aB <- c(Data$B.tilde%*%Alpha)
  
  H <- matrix(NA,nrow=J+p,ncol=J+p)
  H[1:J,1:J] <- Compute.Hess(Alpha = Alpha, B.obs = Data$B.obs, B.tilde.miss = Data$B.tilde[Data$Ind.Miss,], lambda = lambda)
  H[(J+1):(J+p),(J+1):(J+p)] <- Hess.Beta(X.obs = Cov[Data$Ind.obs,], X.NoRSV = Cov[Data$Ind.NoRSV,], X.Miss = Cov[Data$Ind.Miss,], aB.obs = aB[Data$Ind.obs], aB.NoRSV = aB[Data$Ind.NoRSV], aB.Miss = aB[Data$Ind.Miss])
  H[(J+1):(J+p),1:J] <- Hess.BetaAlpha(X.obs = Cov[Data$Ind.obs,], X.NoRSV = Cov[Data$Ind.NoRSV,], X.Miss = Cov[Data$Ind.Miss,], Btilde.obs = Data$B.tilde[Data$Ind.obs,], Btilde.NoRSV = Data$B.tilde[Data$Ind.NoRSV,], Btilde.Miss = Data$B.tilde[Data$Ind.Miss,], Alpha = Alpha)
  H[1:J,(J+1):(J+p)] <- t(H[(J+1):(J+p),1:J])
  return(H)
}

Hess.Beta <- function(X.obs, X.NoRSV, X.Miss, aB.obs, aB.NoRSV, aB.Miss) {
  H1 <- -t(X.obs*aB.obs)%*%X.obs
  tmp <- aB.Miss/(exp(aB.Miss) - 1); tmp <- tmp - tmp^2*exp(aB.Miss)
  H2 <- t(X.Miss*tmp)%*%X.Miss
  H3 <- -t(X.NoRSV*aB.NoRSV)%*%X.NoRSV
  return(H1 + H2 + H3)
}

Hess.BetaAlpha <- function(X.obs, X.NoRSV, X.Miss, Btilde.obs, Btilde.NoRSV, Btilde.Miss, Alpha) {
  H1 <- -t(X.obs)%*%Btilde.obs
  aB.Miss <- c(Btilde.Miss%*%Alpha)
  tmp <- 1/(exp(aB.Miss) - 1) - aB.Miss*exp(aB.Miss)/(exp(aB.Miss) - 1)^2
  H2 <- t(X.Miss*tmp)%*%Btilde.Miss
  H3 <- -t(X.NoRSV)%*%Btilde.NoRSV
  return(H1 + H2 + H3)
}

Get.Distribution <- function(Alpha, beta, DOB, x=rep(0,length(beta)), Data, Lambda) {
  a.max <- 365
  mult <- exp(sum(beta*x))
  knots <- Data$knots; degree <- Data$degree
  B <- splines::splineDesign(knots = knots, x = 1:a.max, ord = degree+1)
  DOB <- as.Date(DOB); Lambda$date <- as.Date(Lambda$date)
  
  t.0 <- DOB
  Lambda.func <- rep(NA,1+a.max)
  for (i in 1:length(Lambda.func)) {
    date.i <- t.0 + i - 1
    time.i <- as.numeric(date.i - Lambda$date)
    ind.i.max <- max(which(time.i>=0))
    if (time.i[ind.i.max]==0) {Lambda.func[i] <- Lambda$lambda[ind.i.max]}
    lambda.1 <- Lambda$lambda[ind.i.max]; lambda.2 <- Lambda$lambda[ind.i.max+1]
    Lambda.func[i] <- lambda.1 + (lambda.2 - lambda.1)/(as.numeric(Lambda$date[ind.i.max+1]-Lambda$date[ind.i.max]))*time.i[ind.i.max]
  }
  
  #First part#
  Prob.a <- c(B%*%Alpha)*mult*Lambda.func[2:length(Lambda.func)]
  
  #Second part#
  Integrals <- sapply(1:ncol(B), function(j){  #A x J 
    ind.j <- rep(0,length(knots)-degree-1)
    ind.j[j] <- 1
    sapply(1:a.max,function(x){
      max(0,IntegrateBs::ibs(x=x, knots = c(rep(0,degree+1),knots), ord = degree+1, coef = c(rep(0,degree+1),ind.j)) - IntegrateBs::ibs(x=x-1, knots = c(rep(0,degree+1),knots), ord = degree+1, coef = c(rep(0,degree+1),ind.j)))
    })
  })
  Tmp.mat <- cbind(rep(1,a.max))%*%rbind(Lambda.func[2:length(Lambda.func)])
  Tmp.mat[upper.tri(x = Tmp.mat, diag = F)] <- 0
  Prob.a <- Prob.a*exp(-mult*c(Tmp.mat%*%Integrals%*%Alpha)); Prob.a <- Prob.a/sum(Prob.a)
  A.vec <- 1:a.max
  mu.a <- sum(Prob.a*A.vec); sd.a <- sqrt( sum(Prob.a*A.vec^2) - mu.a^2 )
  return(list(A=A.vec,Prob.a=Prob.a,mu.a=mu.a,sd.a=sd.a))
}

Combine.Prob.Months <- function(prob) {
  if (!is.list(prob)) {
    m <- cut(1:length(prob), 12, labels=F)
    return(sapply(unique(m),function(i){sum(prob[m==i])}))
  }
  return(sapply(prob,function(p){return(Combine.Prob.Months(p))}))
}

Est.Prob.Step <- function(w, c=rep(1,length(w)), Lambda, DOB, Age.RSV, a.max=365, Month.use=NULL, return.OneYearProb=F) {
  if (!is.null(Month.use)) {
    ind <- which(sapply(strsplit(as.character(DOB),"-"),function(x){as.numeric(x[2])}) %in% Month.use)
    DOB <- DOB[ind]; Age.RSV <- Age.RSV[ind]
  }
  ind.use <- which(!is.na(Age.RSV)); DOB <- as.Date(DOB[ind.use]); Age.RSV <- Age.RSV[ind.use]
  out.return <- list()
  out.return$DOB <- DOB; out.return$Age.RSV <- Age.RSV
  
  n <- length(DOB)
  t.0 <- min(DOB)
  Lambda$date <- as.Date(Lambda$date)
  Lambda.func <- rep(NA,as.numeric(max(DOB)-t.0)+1+a.max)
  for (i in 1:length(Lambda.func)) {
    date.i <- t.0 + i - 1
    time.i <- as.numeric(date.i - Lambda$date)
    ind.i.max <- max(which(time.i>=0))
    if (time.i[ind.i.max]==0) {Lambda.func[i] <- Lambda$lambda[ind.i.max]}
    lambda.1 <- Lambda$lambda[ind.i.max]; lambda.2 <- Lambda$lambda[ind.i.max+1]
    Lambda.func[i] <- lambda.1 + (lambda.2 - lambda.1)/(as.numeric(Lambda$date[ind.i.max+1]-Lambda$date[ind.i.max]))*time.i[ind.i.max]
  }
  Tmp.mat <- sapply(1:n,function(i){  #A x n
    tmp.i <- as.numeric(DOB[i]-t.0)+1
    return(Lambda.func[(tmp.i+1):(tmp.i+a.max)])
  })
  cc <- cut(x = 1:a.max, 12, labels = F)
  Tmp.mat <- w*t(sapply(unique(cc),function(i){colMeans(Tmp.mat[cc==i,])}))  #A x n
  Prob <- sapply(1:12,function(a){
    if (a==1) {return(1 - exp(-Tmp.mat[a,]))}
    return((1 - exp(-Tmp.mat[a,]))*exp(-colSums(rbind(Tmp.mat[1:(a-1),]))))
  })  #n x A
  if (return.OneYearProb) {return(mean(rowSums(Prob)))}
  Prob <- sweep(x = Prob, MARGIN = 2, STATS = c, FUN = "*")/rowSums(Prob)
  return(colMeans(Prob))
}

VarExp.Model <- function(Alpha, knots, degree=3, Beta, DOB.start, DOB.end, Lambda, Z, Var.Beta=NULL) {
  a.max <- 365
  DOB <- as.Date(DOB.start) + 1:(as.Date(DOB.end) - as.Date(DOB.start) + 1)
  Z <- Z[!is.na(rowSums(Z)),]
  BZ <- exp(c(Z%*%Beta))
  if (!is.null(Var.Beta)) {BZ <- BZ/(1+rowSums((Z%*%Var.Beta)*Z))}
  BZ.unique <- unique(BZ); p.BZ <- sapply(BZ.unique,function(b){sum(BZ==b)})/length(BZ); rm(BZ)
  n <- length(DOB)
  Lambda$date <- as.Date(Lambda$date)
  
  Integrals <- sapply(1:length(Alpha), function(j){  #A x J 
    ind.j <- rep(0,length(knots)-degree-1)
    ind.j[j] <- 1
    sapply(1:a.max,function(x){
      max(0,IntegrateBs::ibs(x=x, knots = c(rep(0,degree+1),knots), ord = degree+1, coef = c(rep(0,degree+1),ind.j)) - IntegrateBs::ibs(x=x-1, knots = c(rep(0,degree+1),knots), ord = degree+1, coef = c(rep(0,degree+1),ind.j)))
    })
  })
  IntAlpha <- c(Integrals%*%Alpha)  #Length A
  
  t.0 <- min(DOB)
  Lambda.func <- rep(NA,as.numeric(max(DOB)-t.0)+1+a.max)
  for (i in 1:length(Lambda.func)) {
    date.i <- t.0 + i - 1
    time.i <- as.numeric(date.i - Lambda$date)
    ind.i.max <- max(which(time.i>=0))
    if (time.i[ind.i.max]==0) {Lambda.func[i] <- Lambda$lambda[ind.i.max]}
    lambda.1 <- Lambda$lambda[ind.i.max]; lambda.2 <- Lambda$lambda[ind.i.max+1]
    Lambda.func[i] <- lambda.1 + (lambda.2 - lambda.1)/(as.numeric(Lambda$date[ind.i.max+1]-Lambda$date[ind.i.max]))*time.i[ind.i.max]
  }
  if (return.Lambda) {return(Lambda.func)}
  Tmp.mat <- sapply(1:n,function(i){  #A x n
    tmp.i <- as.numeric(DOB[i]-t.0)+1
    return(Lambda.func[(tmp.i+1):(tmp.i+a.max)])
  })
  
  S.diff <- Tmp.mat*IntAlpha  #A x n
  S.cumulative <- rbind(rep(0,n),apply(S.diff[-nrow(S.diff),],2,cumsum))  #A x n
  tmp <- lapply(BZ.unique,function(b){
    probs <- exp(-b*S.cumulative)*( 1 - exp(-b*S.diff) )
    Prob.1 <- colSums(probs) #Length n 
    probs <- t(t(probs)/Prob.1)
    ER <- colSums(probs*(1:nrow(probs)))
    ER2 <- colSums(probs*(1:nrow(probs))^2)
    return(list(E1=Prob.1,V1=Prob.1*(1-Prob.1),ER=ER,VR=ER2-ER^2))
  })
  ER.B <- cbind(c(sapply(tmp,function(x){x$E1})%*%p.BZ), c(sapply(tmp,function(x){x$ER})%*%p.BZ))
  VER.B <- colMeans(cbind(c(sapply(tmp,function(x){x$E1^2})%*%p.BZ)-ER.B[,1]^2, c(sapply(tmp,function(x){x$ER^2})%*%p.BZ)-ER.B[,2]^2))
  ER.B <- apply(ER.B,2,var)
  Error.VBZ <- c( mean(c(sapply(tmp,function(x){x$V1})%*%p.BZ)), mean(c(sapply(tmp,function(x){x$VR})%*%p.BZ)) )
  
  return(list(PVE.1=100*c(ER.B[1],VER.B[1])/(ER.B[1]+VER.B[1]+Error.VBZ[1]), PVE.R=100*c(ER.B[2],VER.B[2])/(ER.B[2]+VER.B[2]+Error.VBZ[2])))
}