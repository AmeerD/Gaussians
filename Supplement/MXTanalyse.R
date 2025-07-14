library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(latex2exp)

load("Ameer's Notes/T Simulations/mxtresults.rda")

t1 <- fullres %>%
  mutate(LRStat = pchisq(as.numeric(LRStat), df=1, lower.tail=F)) %>%
  filter(!is.na(LRStat)) %>%
  mutate(df = factor(df, levels = c(4,10,100))) %>%
  ggplot(aes(sample=LRStat, colour=df)) +
  stat_qq(distribution=qunif) +
  geom_abline(intercept=0, slope=1) +
  theme(legend.position="bottom") +
  labs(colour="Degrees of Freedom") +
  xlim(c(0,1)) + ylim(c(0,1)) +
  xlab("Significance cutoff") + ylab("P-value quantiles")

t1

ggsave("Ameer's Notes/T Simulations/tdist_t1e.pdf", t1, width=4, height=4)


