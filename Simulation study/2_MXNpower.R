library(dplyr)
library(datathin)
library(mvtnorm)
library(cmdstanr)
library(argparse)

source("MXNfunctions.R")
options(dplyr.summarise.inform = FALSE)
mxnmod <- cmdstanr::cmdstan_model("mxn_Q.stan")

## -----------------------------------------
## Load any command line arguments
## -----------------------------------------
parser <- ArgumentParser()
parser$add_argument("--nreps", type = "double", default = 50,
                    help = "number of replicates for each set of params")
parser$add_argument("--n_node", type = "double", default = 10,
                    help = "number of nodes in the network")
parser$add_argument("--eps", type = "double", default = 0.5,
                    help = "fission tuning parameter")
args <- parser$parse_args()
jobid <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
print(jobid)

nreps <- args$nreps
eps <- args$eps

## -----------------------------------------
## Power simulation
## -----------------------------------------

set.seed(jobid)

filename <- paste("res/power", args$n_node, args$eps, "_", jobid, ".txt", sep="")

# Set the simulation parameters 
# Here we are assuming that K=2, and Sigt=I
n_time <- 100 
n_node <- args$n_node
rho <- 0.9
Delta <- diag(n_node)
Gamma <- ARconstructor(n_time, c(1, rho))

Qtemp <- diag(n_time)*(1+rho^2)
for (i in 1:(n_time-1)) {
  Qtemp[i,i+1] <- Qtemp[i+1,i] <- -rho
}
evecQAR <- eigen(Qtemp)$vectors
evalQARseq <- cos((1:n_time)*pi/(n_time+1))

q1 <- eps^0.25
Q <- matrix(c(q1, sqrt(1-q1^2), sqrt(1-q1^2), -q1), nrow=2, byrow=T)

corvals <- seq(-0.9, 0.9, length.out=19)

for (cval in corvals) {
  Delta[1,2] <- Delta[2,1] <- cval
  
  Sig <- Gamma %x% Delta
  
  sigd <- 1
  
  for (i in 1:nreps) {
    ## Generate and thin the data
    X <- rmvnorm(1, sigma = Sig)
    split <- QalgCDiag(X, 1, Q, sigd) #datathin(X, "mvgaussian", arg=Sigt)
    X1 <- split[1,]
    X2 <- split[2,]
    
    ## Pick out the largest absolute entry in the X1 sample covariance matrix
    Xmat <- matrix(X1, nrow=n_node)
    fissioncor <- cov(t(Xmat))
    diag(fissioncor) <- NA
    zfission <- as.vector(which(abs(fissioncor) == max(abs(fissioncor), na.rm=T), arr.ind=T)[1,])
    
    ## Conduct the conditional test if zfission is the entry 1,2
    if ((zfission[1] == 2) & (zfission[2] == 1)) {
      Xmat <- matrix(X2, nrow=n_node)
      ar1 <- arima(as.vector(t(cbind(Xmat, matrix(NA, nrow=nrow(Xmat), ncol=ncol(Xmat))))), 
                   order=c(1,0,0), include.mean=F)
      initlist <- list(list(rho=ar1$coef[1], cor=cor(t(Xmat))  ))
      
      conddatH0 <- list(n_node=n_node, n_time=n_time, X=X2, Y=X1, mod=3, zero=zfission,
                        evecQAR=evecQAR, seq=evalQARseq, QX=Q[,1])
      conddatH1 <- list(n_node=n_node, n_time=n_time, X=X2, Y=X1, mod=3, zero=c(0,0),
                        evecQAR=evecQAR, seq=evalQARseq, QX=Q[,1])
      
      condH1 <- mxnmod$optimize(data = conddatH1, seed = jobid, refresh=250, init=initlist, 
                                tol_rel_grad=1, #tol_param=1e-10, tol_rel_obj=10,
                                algorithm="lbfgs", iter=3000, sig_figs=8)
      
      initlist <- list(list(rho=condH1$mle("rho"), cor=matrix(condH1$mle("cor"), nrow=n_node)))
      H0test <- matrix(condH1$mle("cor"), nrow=n_node)
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
      
      LRTcond <- tryCatch(-2*(condH0$mle("ll") - condH1$mle("ll")), error=function(err) NA)
    } else {
      LRTcond <- NA
    }
    
    ## Write results to text file
    write(c(cval, eps, i, zfission, LRTcond), file = filename, append=TRUE, ncolumns = 6)
  }
  
  print("***")
  print(paste0("Correlation: ", cval, " complete"))
  print("***")
}






