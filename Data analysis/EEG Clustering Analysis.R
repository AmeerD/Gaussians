library(dplyr)
library(Matrix)
library(datathin)
library(mvtnorm)
library(dlm)
library(patchwork)
library(reshape)
library(ggplot2)
library(ggdendro)
library(dendextend)
library(latex2exp)
library(eegUtils)

source("MXNfunctions.R")
options(dplyr.summarise.inform = FALSE)

set.seed(2024)

eeg <- readRDS("EEG data.RData")

eeg1 <- eeg %>% 
  filter(Trial == 0) %>%
  # filter(!(Channel %in% c("X", "Y", "nd"))) %>%
  group_by(Channel) %>% 
  mutate(Value = Value - mean(Value),
         Value = Value/sd(Value)) %>%
  ungroup %>%
  arrange(Time, Channel)

rm(eeg)

n_time <- length(unique(eeg1$Time))
n_node <- length(unique(eeg1$Channel))

X <- eeg1$Value
n <- length(X)
covX <- cov(t(matrix(X, nrow=n_node))) 

## SPLIT THE DATA

q1 <- 0.5^0.25
sigd <- 1
Q <- matrix(c(q1, sqrt(1-q1^2), sqrt(1-q1^2), -q1), nrow=2, byrow=T)
split <- QalgCDiag(matrix(X, nrow=1), 1, Q, sigd)
X1 <- split[1,]
X2 <- split[2,]

## VALIDATION STEP
### Rescaling is challening since we can't guarantee positive definiteness - set all negative eigenvalues to 0.001
eigX1 <- eigen((cov(t(matrix(X1, nrow=n_node))) - (1-q1^2)*diag(n_node))/(q1^2))
covX1 <- eigX1$vectors %*% diag(pmax(eigX1$values, 0.01)) %*% t(eigX1$vectors)
### We work with correlation matrices since after standardizing we know that Delta has unit diagonal
corX1 <- diag(1/sqrt(diag(covX1))) %*% covX1 %*% diag(1/sqrt(diag(covX1))) 
rownames(corX1) <- unique(eeg1$Channel)
ar1X1 <- dlmMLE(as.vector(t(cbind(matrix(X1, nrow=n_node), 
                                  matrix(NA, nrow=n_node, ncol=n_time)))), 
                c(1,0.5), buildAR1, evar=(1-q1^2))

## Apply hierarchical clustering to the correlation matrix
hc <- hclust(as.dist(1-corX1))
plot(hc, hang=-1)

## Use X2|X1 to identify the number of clusters to cut from the tree
height <- c(1:n_node)
cuts <- cutree(hc, k=height)
ll <- rep(NA, length(height))
nonzeros <- rep(NA, length(height))

corclust <- function(cormat, clust) {
  for (i in 1:length(clust)) {
    cormat[i,which(clust != clust[i])] <- 0
  }
  return(cormat)
}

for (h in height) {
  cortemp <- corclust(corX1, cuts[,h])
  
  nonzeros[h] <- sum((cortemp[lower.tri(cortemp)]) == 0)/length(cortemp[lower.tri(cortemp)])

  ### Conditional method
  ll[h] <- QmxnX2X1(X2, X1, ar1X1$par[2], cortemp, n_time, Q=Q) 
  
  print(paste0("Done: ", h))
}

hhat <- height[which(ll == max(ll))[1]]
hhat


## Save plots
dendplot <- ggplot(as.ggdend(hc %>% as.dendrogram %>% set("branches_k_color", k=hhat))) +
  theme_dendro()  + 
  ylim(c(-0.15,NA))

ggsave("dendplot.pdf", dendplot, width=12, height=4)

cll1 <- data.frame(x=1:n_node, y=ll) %>%
  ggplot(aes(x=x, y=y)) +
  geom_hline(yintercept=ll[hhat]) +
  geom_vline(xintercept=hhat) +
  geom_line(colour="steelblue") +
  geom_point(colour="steelblue", size=1) +
  xlab("Number of clusters") +
  ylab(unname(TeX("CLL of $X^{(2)}|X^{(1)}$")))

ggsave("cll1.pdf", cll1, width=4, height=4)

dendcols <- unique(ggplot_build(dendplot)$data[1][[1]] %>% arrange(group) %>% pull(colour))[1:hhat]

brain <- electrode_locations(data.frame(electrode=unique(eeg1$Channel), cluster=as.factor(cuts[,hhat]))) %>%
  ggplot(aes(x=x, y=y)) + 
  geom_text(aes(label=electrode, colour=cluster), fontface="bold", size=2.5) +
  theme_minimal() +
  theme(legend.position="none")

ggsave("brain.pdf", brain, width=4, height=4)

