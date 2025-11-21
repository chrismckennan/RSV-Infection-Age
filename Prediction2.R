PredictCondDensity <- function(Lambda, DOB, IDs=NULL, Age.RSV=NULL, a.max=365, knots, degree=3, Alpha, Cov=NULL, Beta=NULL, Mu=NULL, return.Lambda=F, Month=NULL, Month.use=NULL, Cond.Prob=T) {
  if (is.list(DOB)) {
    if (is.null(IDs)) {return(lapply(1:length(DOB),function(j){PredictCondDensity(Lambda = Lambda[[j]], DOB = DOB[[j]], Age.RSV = NULL, a.max = a.max, knots = knots, degree = degree, Alpha = Alpha, Cov = Cov[[j]], Beta = Beta, Mu = Mu, return.Lambda = return.Lambda, Month = Month, Month.use = Month.use, Cond.Prob = Cond.Prob)}))}
    if (!is.null(IDs)) {return(lapply(1:length(DOB),function(j){PredictCondDensity(Lambda = Lambda[[j]], DOB = DOB[[j]], Age.RSV = NULL, a.max = a.max, knots = knots, degree = degree, Alpha = Alpha, Cov = Cov[[j]], Beta = Beta, Mu = Mu, return.Lambda = return.Lambda, Month = Month, Month.use = Month.use, Cond.Prob = Cond.Prob, IDs = IDs[[j]])}))}
  }
  if (is.null(Age.RSV)) {Age.RSV <- rep(NA,length(DOB))}
  out.return <- list()
  if (!is.null(Month.use)) {
    ind <- which(sapply(strsplit(as.character(DOB),"-"),function(x){as.numeric(x[2])}) %in% Month.use)
    DOB <- DOB[ind]; Age.RSV <- Age.RSV[ind]
    if (!is.null(IDs)) {IDs <- IDs[ind]}
    if (!is.null(Cov)) {Cov <- cbind(cbind(Cov)[ind,])}
    if (!is.null(Month)) {Month <- Month[ind]}
  }
  DOB <- as.Date(DOB); Lambda$date <- as.Date(Lambda$date)
  if (!is.null(Cov) & !is.null(Beta)) {
    ind.use <- which(!is.na(rowSums(cbind(Cov))))
    DOB <- DOB[ind.use]
    Age.RSV <- Age.RSV[ind.use]
    Cov <- cbind(cbind(Cov)[ind.use,])
    if (!is.null(IDs)) {IDs <- IDs[ind.use]}
  }
  if (!is.null(Month)) {out.return$Month <- Month[ind.use]}
  out.return$DOB <- DOB; out.return$Age.RSV <- Age.RSV
  if (!is.null(IDs)) {out.return$IDs <- IDs}
  n <- length(DOB)
  
  Integrals <- sapply(1:length(Alpha), function(j){  #A x J 
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
  if (return.Lambda) {return(Lambda.func)}
  Tmp.mat <- sapply(1:n,function(i){  #A x n
    tmp.i <- as.numeric(DOB[i]-t.0)+1
    return(Lambda.func[(tmp.i+1):(tmp.i+a.max)])
  })
  S <- c(t(Tmp.mat)%*%Integrals%*%Alpha) #n

  if (!is.null(Cov) & !is.null(Beta)) {
    Cov <- Cov - cbind(rep(1,NROW(Cov)))%*%rbind(Mu)
    tmp <- c(cbind(Cov)%*%Beta)
    c <- exp(tmp)
  } else {
    c <- rep(1,n)
  }
  Denom <- 1 - exp(-c*S)

  #Numerator#
  tmp <- exp(-c*t(apply(X = c(Integrals%*%Alpha) * Tmp.mat, MARGIN = 2, cumsum)))  #X is A x n, output is n x A
  if (Cond.Prob) {
    out.return$Prob <- (cbind(rep(1,n),tmp[,-a.max]) - tmp)/Denom
  } else {
    out.return$Prob <- (cbind(rep(1,n),tmp[,-a.max]) - tmp)
  }
  return(out.return)
}