library(dplyr)
library(datathin)
library(mvtnorm)
library(cmdstanr)
library(argparse)
library(dlm)

source("MXNfunctions.R")
options(dplyr.summarise.inform = FALSE)
mxnmod <- cmdstanr::cmdstan_model("mxn_Q.stan")

## -----------------------------------------
## Load any command line arguments
## -----------------------------------------
parser <- ArgumentParser()
parser$add_argument("--simname", default = "null",
                    help = "name of simulation")
parser$add_argument("--nreps", type = "double", default = 50,
                    help = "number of replicates for each set of params")
parser$add_argument("--n_node", type = "double", default = 10,
                    help = "number of nodes in the network")
args <- parser$parse_args()
jobid <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
print(jobid)

nreps <- args$nreps

## -----------------------------------------
## Null simulation
## -----------------------------------------

set.seed(jobid)

filename <- paste("res/", args$simname, args$n_node, args$eps, "_", jobid, ".txt", sep="")

# Set the simulation parameters 
# Here we are assuming that Q=[[sqrt(eps), sqrt(1-eps)],[sqrt(1-eps), -eps]], r=1, and Sigt=diag(sigd)
n_time <- 100 
n_node <- args$n_node
rho <- 0.9

Delta <- diag(n_node)
if (args$simname != "null") {
  Delta[1,2] <- Delta[2,1] <- 0.5
}
Gamma <- AR1constructor(n_time, rho)
Sig <- Gamma %x% Delta

sigd <- 1

Qtemp <- diag(n_time)*(1+rho^2)
for (i in 1:(n_time-1)) {
  Qtemp[i,i+1] <- Qtemp[i+1,i] <- -rho
}
evecQAR <- eigen(Qtemp)$vectors
evalQARseq <- cos((1:n_time)*pi/(n_time+1))

q1 <- 0.5^0.25
Q <- matrix(c(q1, sqrt(1-q1^2), sqrt(1-q1^2), -q1), nrow=2, byrow=T)

for (i in 1:nreps) {
  ## Generate and decompose the data
  X <- rmvnorm(1, sigma = Sig)
  split <- QalgCDiag(X, 1, Q, sigd) #datathin(X, "mvgaussian", arg=Sigt)
  X1 <- split[1,]
  X2 <- split[2,]
  
  ## Naive method
  Xmat <- matrix(X, nrow=n_node)
  ar1 <- arima(as.vector(t(cbind(Xmat, matrix(NA, nrow=nrow(Xmat), ncol=ncol(Xmat))))), 
               order=c(1,0,0), include.mean=F)
  
  ## Pick out the largest absolute entry in the sample covariance matrix
  naivecor <- cov(t(Xmat))
  diag(naivecor) <- NA
  znaive <- as.vector(which(abs(naivecor) == max(abs(naivecor), na.rm=T), arr.ind=T)[1,])
  
  ### Maximize the likelihood under the null and the alternative
  naivedatH0 <- list(n_node=n_node, n_time=n_time, X=as.vector(X), Y=rep(0,length(X)), mod=0, zero=znaive,
                     evecQAR=evecQAR, seq=evalQARseq, QX=Q[,1])
  naivedatH1 <- list(n_node=n_node, n_time=n_time, X=as.vector(X), Y=rep(0,length(X)), mod=0, zero=c(0,0),
                     evecQAR=evecQAR, seq=evalQARseq, QX=Q[,1])
  
  ### Check if inserting a zero into the sample correlation matrix will be a valid starting point for H0
  initlist <- list(list(rho=ar1$coef[1], cor=cor(t(Xmat))))
  
  naiveH1 <- mxnmod$optimize(data = naivedatH1, seed = jobid, refresh=250, init=initlist, 
                             tol_rel_grad=1, #tol_param=1e-10, tol_rel_obj=10,
                             algorithm="lbfgs", iter=3000, sig_figs=8)
  
  initlist <- list(list(rho=naiveH1$mle("rho"), cor=matrix(naiveH1$mle("cor"), nrow=n_node)))
  H0test <- matrix(naiveH1$mle("cor"), nrow=n_node) #cor(t(Xmat))
  H0test[znaive[1], znaive[2]] <- H0test[znaive[2], znaive[1]] <- 0
  if (min(eigen(H0test, only.values=T)$values) > 0) {
    initlist[[1]]$cor <- H0test
  } else {
    H0test.eig <- eigen(H0test)
    H0temp <- H0test.eig$vectors %*% diag(pmax(H0test.eig$values, 0.00001)) %*% t(H0test.eig$vectors)
    H0temp <- diag(1/sqrt(diag(H0temp))) %*% H0temp %*% diag(1/sqrt(diag(H0temp)))
    diag(H0temp) <- rep(1, n_node)
    initlist <- list(list(rho=naiveH1$mle("rho"), cor=H0temp))
  }
  naiveH0 <- mxnmod$optimize(data = naivedatH0, seed = jobid, refresh=250, init=initlist, 
                             tol_rel_grad=1, #tol_param=1e-10, tol_rel_obj=10,
                             algorithm="lbfgs", iter=3000, sig_figs=8)
  
  
  LRTnaive <- tryCatch(-2*(naiveH0$mle("ll") - naiveH1$mle("ll")), error=function(err) NA)
  
  ## Fission
  Xmat <- matrix(X1, nrow=n_node)
  
  ## Pick out the largest absolute entry in the sample covariance matrix
  fissioncor <- cov(t(Xmat))
  diag(fissioncor) <- NA
  zfission <- as.vector(which(abs(fissioncor) == max(abs(fissioncor), na.rm=T), arr.ind=T)[1,])
  
  ## Marginal and Conditional tests
  Xmat <- matrix(X2, nrow=n_node)
  ar1 <- arima(as.vector(t(cbind(Xmat, matrix(NA, nrow=nrow(Xmat), ncol=ncol(Xmat))))), 
               order=c(1,0,0), include.mean=F)
  initlist <- list(list(rho=ar1$coef[1], cor=cor(t(Xmat))  ))
  margdatH0 <- list(n_node=n_node, n_time=n_time, X=X2, Y=rep(0,length(X)), mod=2, zero=zfission,
                    evecQAR=evecQAR, seq=evalQARseq, QX=Q[,1])
  margdatH1 <- list(n_node=n_node, n_time=n_time, X=X2, Y=rep(0,length(X)), mod=2, zero=c(0,0),
                    evecQAR=evecQAR, seq=evalQARseq, QX=Q[,1])
  conddatH0 <- list(n_node=n_node, n_time=n_time, X=X2, Y=X1, mod=3, zero=zfission,
                    evecQAR=evecQAR, seq=evalQARseq, QX=Q[,1])
  conddatH1 <- list(n_node=n_node, n_time=n_time, X=X2, Y=X1, mod=3, zero=c(0,0),
                    evecQAR=evecQAR, seq=evalQARseq, QX=Q[,1])
  
  margH1 <- mxnmod$optimize(data = margdatH1, seed = jobid, refresh=250, init=initlist, 
                            tol_rel_grad=1, #tol_param=1e-10, tol_rel_obj=10,
                            algorithm="lbfgs", iter=3000, sig_figs=8)
  condH1 <- mxnmod$optimize(data = conddatH1, seed = jobid, refresh=250, init=initlist, 
                            tol_rel_grad=1, #tol_param=1e-10, tol_rel_obj=10,
                            algorithm="lbfgs", iter=3000, sig_figs=8)
  
  # initlist <- list(list(rho=margH1$mle("rho"), cor=matrix(margH1$mle("cor"), nrow=n_node)))
  H0test <- matrix(margH1$mle("cor"), nrow=n_node) #cor(t(Xmat))
  H0test[zfission[1], zfission[2]] <- H0test[zfission[2], zfission[1]] <- 0
  if (min(eigen(H0test, only.values=T)$values) > 0) {
    initlist <- list(list(rho=margH1$mle("rho"), cor=H0test))
  } else {
    H0test.eig <- eigen(H0test)
    H0temp <- H0test.eig$vectors %*% diag(pmax(H0test.eig$values, 0.00001)) %*% t(H0test.eig$vectors)
    H0temp <- diag(1/sqrt(diag(H0temp))) %*% H0temp %*% diag(1/sqrt(diag(H0temp)))
    diag(H0temp) <- rep(1, n_node)
    initlist <- list(list(rho=margH1$mle("rho"), cor=H0temp))
  }
  margH0 <- mxnmod$optimize(data = margdatH0, seed = jobid, refresh=250, init=initlist, 
                            tol_rel_grad=1, #tol_param=1e-10, tol_rel_obj=10,
                            algorithm="lbfgs", iter=3000, sig_figs=8)
  
  # initlist <- list(list(rho=condH1$mle("rho"), cor=matrix(condH1$mle("cor"), nrow=n_node)))
  H0test <- matrix(condH1$mle("cor"), nrow=n_node) #cor(t(Xmat))
  H0test[zfission[1], zfission[2]] <- H0test[zfission[2], zfission[1]] <- 0
  if (min(eigen(H0test, only.values=T)$values) > 0) {
    initlist <- list(list(rho=condH1$mle("rho"), cor=H0test))
  } else {
    H0test.eig <- eigen(H0test)
    H0temp <- H0test.eig$vectors %*% diag(pmax(H0test.eig$values, 0.00001)) %*% t(H0test.eig$vectors)
    H0temp <- diag(1/sqrt(diag(H0temp))) %*% H0temp %*% diag(1/sqrt(diag(H0temp)))
    diag(H0temp) <- rep(1, n_node)
    initlist <- list(list(rho=condH1$mle("rho"), cor=H0temp))
  }
  condH0 <- mxnmod$optimize(data = conddatH0, seed = jobid, refresh=250, init=initlist, 
                            tol_rel_grad=1, #tol_param=1e-10, tol_rel_obj=10,
                            algorithm="lbfgs", iter=3000, sig_figs=8)
  
  ## Marginal test
  LRTmarg <- tryCatch(-2*(margH0$mle("ll") - margH1$mle("ll")), error=function(err) NA)
  
  ## Conditional test
  LRTcond <- tryCatch(-2*(condH0$mle("ll") - condH1$mle("ll")), error=function(err) NA)
  
  ## Write results to text file
  write(c("Naive", i, znaive, LRTnaive), file = filename, append=TRUE, ncolumns = 5)
  write(c("Marginal", i, zfission, LRTmarg), file = filename, append=TRUE, ncolumns = 5)
  write(c("Conditional", i, zfission, LRTcond), file = filename, append=TRUE, ncolumns = 5)
  
  print("***")
  print(paste0("Iteration: ", i, " complete"))
  print(paste0(LRTnaive, "; ", LRTmarg, "; ", LRTcond, "; "))
  print("***")
}
