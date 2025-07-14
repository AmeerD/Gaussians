ARconstructor <- function(p, theta) {
  if (length(theta) == 1) {
    return(diag(rep(theta[1], p)))
  } else {
    Sig <- matrix(0, nrow=p, ncol=p)
    covs <- ts.extend::ARMA.autocov(p, ar=theta[-1])
    
    for (i in 1:p) {
      Sig[abs(row(Sig) - col(Sig)) == i-1] <- covs[i]
    }
    
    return(theta[1]*Sig)
  }
}

AR1constructor <- function(p, rho) {
  Qtemp <- diag(p)*(1+rho^2)
  for (i in 1:(p-1)) {
    Qtemp[i,i+1] <- Qtemp[i+1,i] <- -rho
  }
  return(solve(Qtemp))
}

buildAR1 <- function(par, evar) {
  sigma2 <- par[1]
  ar <- par[2]
  dlmModARMA(ar=ar, sigma2=sigma2, dV=evar)
}

### Q Alg related functions

Qalg <- function(X, r, Q, Sigt) {
  # Assume that X is an nxp matrix
  p <- ncol(X) 
  
  W <- mvtnorm::rmvnorm(r, sigma=Sigt, checkSymmetry=F)
  Xaug <- rbind(X,W)
  
  return(Q %*% Xaug)
}

## Use this version for constant diagonal Sigt
QalgCDiag <- function(X, r, Q, sig) {
  # Assume that X is an nxp matrix
  p <- ncol(X) 
  
  W <- matrix(rnorm(r*p, sd=sqrt(sig)), nrow=r, ncol=p)
  Xaug <- rbind(X,W)
  
  return(Q %*% Xaug)
}

## THESE FUNCTIONS ASSUME CONSTANT DIAGONAL SIGT AND ZERO MEAN AND AR1 PROCESS
QmxnX <- function(X, rho, cor, p, Q=NULL, ARvar=1, vard=1) {
  eigGamma <- eigen(ARvar*AR1constructor(p, rho)) #eigen(ARconstructor(p, c(ARvar, rho)))
  eigDelta <- eigen(cor)
  
  QR <- eigGamma$vectors %x% eigDelta$vectors
  AB <- diag(diag(eigGamma$values) %x% diag(eigDelta$values))
  
  if (min(AB) < 0) {
    return(-Inf)
  } else {
    return(sum(dnorm(t(QR) %*% X, sd=sqrt(AB), log=T)))
  }
}

QmxnX1 <- function(X1, rho, cor, p, Q, ARvar=1, vard=1) {
  q1 <- Q[1,1]
  eigGamma <- eigen(ARvar*AR1constructor(p, rho)) #eigen(ARconstructor(p, c(ARvar, rho)))
  eigDelta <- eigen(cor)
  
  QR <- eigGamma$vectors %x% eigDelta$vectors
  AB <- diag(diag(eigGamma$values) %x% diag(eigDelta$values))
  
  if (min(AB) < 0) {
    print(min(AB))
    return(-Inf)
  } else {
    return(sum(dnorm(t(QR) %*% X1, sd=sqrt(q1*q1*AB+(1-q1^2)*vard), log=T)))
  }
}

QmxnX2 <- function(X2, rho, cor, p, Q, ARvar=1, vard=1) {
  q2 <- Q[2,1]
  eigGamma <- eigen(ARvar*AR1constructor(p, rho)) #eigen(ARconstructor(p, c(ARvar, rho)))
  eigDelta <- eigen(cor)
  
  QR <- eigGamma$vectors %x% eigDelta$vectors
  AB <- diag(diag(eigGamma$values) %x% diag(eigDelta$values))
  
  if (min(AB) < 0) {
    print(min(AB))
    return(-Inf)
  } else {
    return(sum(dnorm(t(QR) %*% X2, sd=sqrt(q2*q2*AB+(1-q2^2)*vard), log=T)))
  }
}

QmxnX2X1 <- function(X2, X1, rho, cor, p, Q, ARvar=1, vard=1) {
  q1 <- Q[1,1]
  q2 <- Q[2,1]
  eigGamma <- eigen(ARvar*AR1constructor(p, rho)) #eigen(ARconstructor(p, c(ARvar, rho)))
  eigDelta <- eigen(cor)
  
  QR <- eigGamma$vectors %x% eigDelta$vectors
  AB <- diag(diag(eigGamma$values) %x% diag(eigDelta$values))
  
  if (min(AB) < 0) {
    print(min(AB))
    return(-Inf)
  } else {
    return(sum(dnorm(t(QR) %*% X2, mean=(q1*q2)*((AB-vard)/(q1*q1*(AB-vard)+vard)) * (t(QR) %*% X1),
                     sd=pmax(sqrt(vard*AB/(q1*q1*(AB-vard)+vard)), 0.1), log=T)))
  }
}

QmxnJoint <- function(X1, X2, rho, cor, p, Q, ARvar=1, vard=1) {
  QXQXT <- Q[,1] %*% t(Q[,1])
  eigGamma <- eigen(ARvar*AR1constructor(p, rho)) #eigen(ARconstructor(p, c(ARvar, rho)))
  eigDelta <- eigen(cor)
  
  QR <- eigGamma$vectors %x% eigDelta$vectors
  AB <- diag(eigGamma$values) %x% diag(eigDelta$values)
  S <- QR %*% AB %*% t(QR)
  
  if (min(AB) < 0) {
    print(min(AB))
    return(-Inf)
  } else {
    return(dmvnorm(c(X1,X2),sigma=QXQXT %x% S + (diag(2)-QXQXT) %x% (vard*diag(dim(S)[1])), log=T))
  }
}
