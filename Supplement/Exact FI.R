library(dplyr)
library(tidyr)
library(datathin)
library(mvtnorm)
library(numDeriv)
library(purrr)
library(ggplot2)
library(latex2exp)

ARconstructor <- function(theta, p) {
  if (length(theta) == 0) {
    return(diag(rep(1, p)))
  } else {
    Sig <- matrix(0, nrow=p, ncol=p)
    covs <- ts.extend::ARMA.autocov(p, ar=theta)
    
    for (i in 1:p) {
      Sig[abs(row(Sig) - col(Sig)) == i-1] <- covs[i]
    }
    
    return(Sig)
  }
}

AR1constructor <- function(theta, p) {
  Qtemp <- diag(p)*(1+theta^2)
  Qtemp[abs(row(Qtemp) - col(Qtemp)) == 1] <- -theta
  return(solve(Qtemp))
}

EXCHconstructor <- function(theta, p) {
  Sig <- (1-theta)*diag(p) + matrix(theta, nrow=p, ncol=p)
  return(Sig)
}

## banded=T assumes symmetry!!!
gradSig <- function(SigFun, pars, p, banded=F) {
  npars <- length(pars)
  
  SigFunCoord <- function(pars, p, xcoord, ycoord) {
    return(SigFun(pars, p)[xcoord, ycoord])
  }
  
  dQ <- replicate(npars, matrix(NA, p, p), simplify=F) #array(NA, dim=c(npars, p, p))
  if (banded) {
    for (i in 1:p) {
      temp <- grad(SigFunCoord, pars, p=p, xcoord=1, ycoord=i)
      for (k in 1:npars) {
        dQ[[k]][abs(row(dQ[[k]]) - col(dQ[[k]])) == (i-1)] <- temp[k]
      }
    }
  } else {
    for (i in 1:p) {
      for (j in 1:p) {
        temp <- grad(SigFunCoord, pars, p=p, xcoord=i, ycoord=j)
        for (k in 1:npars) {
          dQ[[k]][i,j] <- temp[k]
        }
      }
    }
  }
  
  return(dQ)
}

FIexact <- function(SigFun, pars, p, q1, Sigp, banded=F) {
  npars <- length(pars) 
  
  Sig <- SigFun(pars, p)
  Prec <- solve(Sig)
  Prec1 <- solve((q1^2)*Sig + (1-q1^2)*Sigp)
  
  dSig <- gradSig(SigFun, pars, p, banded)
  
  IX <- matrix(NA, nrow=npars, ncol=npars)
  IX1 <- matrix(NA, nrow=npars, ncol=npars)
  
  for (i in 1:npars) {
    for (j in 1:i) {
      IX[i,j] <- IX[j,i] <- 0.5*sum(diag(Prec %*% dSig[[i]] %*% Prec %*% dSig[[j]]))
      IX1[i,j] <- IX1[j,i] <- 0.5*(q1^4)*sum(diag(Prec1 %*% dSig[[i]] %*% Prec1 %*% dSig[[j]]))
      
    }
  }
  
  return(list(IX,IX1))
}

# AR(1) Experiment
## Goal: Allocate 50% of the FI in X about rho into X1 by choosing c
ptest <- 100
gamma <- 0.5
q1 <- 0.5^(1/4)
## Grid of rho and c values
AR1pars <- expand_grid(
  rho = seq(-0.9, 0.9, by=0.005),
  c = seq(0.01, 6, by=0.005)
) %>% 
  rowwise %>% 
  mutate(
    FI.X = 0.5*sum(((2*rho+2*cos((1:ptest)*pi/(ptest+1)))/(1+rho^2+2*rho*cos((1:ptest)*pi/(ptest+1))))^2),
    FI.X1 = ((q1^4)/2)*sum(((2*rho+2*cos((1:ptest)*pi/(ptest+1)))/((q1*q1/(1+rho^2+2*rho*cos((1:ptest)*pi/(ptest+1)))+(1-q1*q1)*c)*(1+rho^2+2*rho*cos((1:ptest)*pi/(ptest+1)))^2))^2) ) 

FI.AR1 <- AR1pars %>% 
  mutate(prop = FI.X1/FI.X) %>%
  ggplot(aes(x=rho, y=c, fill=prop)) +
  geom_tile() +
  geom_line(data=AR1pars %>% mutate(prop = FI.X1/FI.X) %>%
              group_by(rho) %>% filter(abs(prop - 0.6) == min(abs(prop - 0.6))),
            aes(x=rho, y=c)) +
  geom_line(data=AR1pars %>% mutate(prop = FI.X1/FI.X) %>%
              group_by(rho) %>% filter(abs(prop - 0.4) == min(abs(prop - 0.4))),
            aes(x=rho, y=c)) +
  scale_fill_gradient2(low="red", high="blue", midpoint=0.5, lim=c(0,1)) +
  xlab(TeX("\\phi")) + ylab(TeX(r"($\sigma'^2$)")) +
  theme(legend.position="bottom") +
  labs(fill="") 
FI.AR1
# ggtitle("Fisher information about AR(1) parameter allocated to X1")

ggsave("FIar1.png", FI.AR1, width=4, height=5)

phitilde <- 0.25
sigtilde.AR1 <- AR1pars %>%
  ungroup %>% 
  filter(between(rho, phitilde-0.001, phitilde+0.001)) %>% 
  mutate(prop=FI.X1/FI.X) %>% 
  filter(abs(prop - 0.5) == min(abs(prop - 0.5))) %>% 
  pull(c)
sigrange.AR1 <- AR1pars %>% 
  ungroup %>% 
  filter(c == sigtilde.AR1) %>% 
  mutate(prop = FI.X1/FI.X) %>% 
  filter(between(prop, 0.4, 0.6))

AR1tilde <- AR1pars %>% 
  mutate(prop = FI.X1/FI.X) %>%
  ggplot(aes(x=rho, y=c, fill=prop)) +
  geom_tile() +
  geom_segment(x=phitilde, y=-Inf, xend=phitilde, yend=sigtilde.AR1, linetype=2, linewidth=0.25) +
  geom_segment(x=-Inf, y=sigtilde.AR1, 
               xend=pull(sigrange.AR1 %>% filter(rho == min(rho)), rho), yend=sigtilde.AR1, linetype=2, linewidth=0.25) +
  geom_segment(x=pull(sigrange.AR1 %>% filter(rho == min(rho)), rho), y=sigtilde.AR1, 
               xend=pull(sigrange.AR1 %>% filter(rho == max(rho)), rho), yend=sigtilde.AR1) +
  geom_point(x=phitilde, y=sigtilde.AR1) +
  annotate("text", y=sigtilde.AR1, x=-Inf, label=TeX(r"($\hat{\sigma}'^2$)"), hjust=1.2) +
  annotate("text", y=-Inf, x=phitilde, label=TeX(r"($\tilde{\phi}=0.25$)"), vjust=1.2) +
  scale_fill_gradient2(low="red", high="blue", midpoint=0.5, lim=c(0,1)) +
  xlab(TeX("\\phi")) + ylab(TeX(r"($\sigma'^2$)")) +
  theme(legend.position="bottom") +
  labs(fill="") +
  coord_cartesian(clip="off")

ggsave("Setar1.png", AR1tilde, width=4, height=5)

# Exchangeable Experiment
## Goal: Allocate 50% of the FI in X about rho into X1 by choosing c
ptest <- 100
gamma <- 0.5
q1 <- 0.5^(1/4)
## Grid of rho and c values
EXCHpars <- expand_grid(
  rho = seq(0.01, 0.99, by=0.005),
  c = seq(0.01, 6, by=0.005)
) %>% 
  mutate(FI.X = (1/(2*(1-rho)^2))*((ptest/(1-rho+rho*ptest))^2-2*ptest/(1-rho+rho*ptest)+ptest),
         A = q1*q1*(1-rho)+(1-q1^2)*c,
         FI.X1 = (ptest*(q1^4)/(2*A^2))*(ptest*((A+rho*q1^2)/(A+rho*ptest*q1^2))^2-2*(A+rho*q1^2)/(A+rho*ptest*q1^2)+1)) %>% 
  select(-A)


FI.EXCH <- EXCHpars %>% 
  mutate(prop = FI.X1/FI.X) %>%
  ggplot(aes(x=rho, y=c, fill=prop)) +
  geom_tile() +
  geom_line(data=EXCHpars %>% mutate(prop = FI.X1/FI.X) %>%
              group_by(rho) %>% filter(abs(prop - 0.6) == min(abs(prop - 0.6))),
            aes(x=rho, y=c)) +
  geom_line(data=EXCHpars %>% mutate(prop = FI.X1/FI.X) %>%
              group_by(rho) %>% filter(abs(prop - 0.4) == min(abs(prop - 0.4))),
            aes(x=rho, y=c)) +
  scale_fill_gradient2(low="red", high="blue", midpoint=0.5, lim=c(0,1)) +
  xlab(TeX("\\phi")) + ylab(TeX(r"($\sigma'^2$)")) +
  theme(legend.position="bottom") +
  labs(fill="") 
  # ggtitle("Fisher information about AR(1) parameter allocated to X1")
FI.EXCH

ggsave("FIexch.png", FI.EXCH, width=4, height=5)

phitilde <- 0.25
sigtilde.exch <- EXCHpars %>%
  ungroup %>% 
  filter(between(rho, phitilde-0.001, phitilde+0.001)) %>% 
  mutate(prop=FI.X1/FI.X) %>% 
  filter(abs(prop - 0.5) == min(abs(prop - 0.5))) %>% 
  pull(c)
sigrange.exch <- EXCHpars %>% 
  ungroup %>% 
  filter(c == sigtilde.exch) %>% 
  mutate(prop = FI.X1/FI.X) %>% 
  filter(between(prop, 0.4, 0.6))

EXCHtilde <- EXCHpars %>% 
  mutate(prop = FI.X1/FI.X) %>%
  ggplot(aes(x=rho, y=c, fill=prop)) +
  geom_tile() +
  geom_segment(x=phitilde, y=-Inf, xend=phitilde, yend=sigtilde.exch, linetype=2, linewidth=0.25) +
  geom_segment(x=-Inf, y=sigtilde.exch, 
               xend=pull(sigrange.exch %>% filter(rho == min(rho)), rho), yend=sigtilde.exch, linetype=2, linewidth=0.25) +
  geom_segment(x=pull(sigrange.exch %>% filter(rho == min(rho)), rho), y=sigtilde.exch, 
               xend=pull(sigrange.exch %>% filter(rho == max(rho)), rho), yend=sigtilde.exch) +
  geom_point(x=phitilde, y=sigtilde.exch) +
  annotate("text", y=sigtilde.exch, x=-Inf, label=TeX(r"($\hat{\sigma}'^2$)"), hjust=1.2) +
  annotate("text", y=-Inf, x=phitilde, label=TeX(r"($\tilde{\phi}=0.25$)"), vjust=1.2) +
  scale_fill_gradient2(low="red", high="blue", midpoint=0.5, lim=c(0,1)) +
  scale_x_continuous(breaks=c(0,0.5,0.75,1)) +
  xlab(TeX("\\phi")) + ylab(TeX(r"($\sigma'^2$)")) +
  theme(legend.position="bottom") +
  labs(fill="") +
  coord_cartesian(clip="off")
EXCHtilde

ggsave("Setexch.png", EXCHtilde, width=4, height=5)

# ggsave("fisher.png", fisher, width=4, height=5)



