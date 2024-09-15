library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

load("mxnresults.rda")
load("mxnpowerresults.rda")

nullres <- fullres %>%
  filter(grepl("null", sim)) 

f1 <- nullres %>%
  mutate(LRStat = pchisq(as.numeric(LRStat), df=1)) %>%
  filter(!is.na(LRStat)) %>%
  mutate(method = case_when(
    method == "Naive" ~ "Method (a)",
    method == "Marginal" ~ "Method (b)",
    TRUE ~ "Method (c)"
  )) %>%
  ggplot(aes(sample=LRStat, colour=method)) +
  stat_qq(distribution=qunif) +
  geom_abline(intercept=0, slope=1) +
  theme(legend.title=element_blank(),
        legend.position="bottom") +
  xlab("Significance cutoff") + ylab("Null hypotheses rejected")

f2 <- powerres %>%
  mutate(LRStat = as.numeric(LRStat)) %>%
  filter(!is.na(LRStat)) %>% 
  mutate(rej = LRStat > qchisq(0.95, 1)) %>%
  group_by(cor, eps) %>%
  summarise(cpower = sum(rej)/n()) %>%
  mutate(epsilon = paste0("(", eps, ",", 1-eps,")")) %>%
  ggplot(aes(x=cor, y=cpower, linetype=epsilon, shape=epsilon)) +
  geom_line(colour="#619CFF") +
  geom_point(colour="#619CFF") +
  scale_shape_manual(values=c(2, 3, 19)) +
  ylim(c(0,1)) + 
  xlab("Correlation between nodes 1 and 2") +
  ylab("Conditional power") +
  theme(legend.position="bottom")
  

f3 <- powerres %>%
  mutate(detected = !is.na(as.numeric(LRStat))) %>%
  group_by(cor, eps) %>% 
  summarise(dprob = sum(detected)/n()) %>%
  mutate(epsilon = paste0("(", eps, ",", 1-eps,")")) %>%
  ggplot(aes(x=cor, y=dprob, linetype=epsilon, shape=epsilon)) +
  geom_point(colour="#619CFF") +
  geom_line(colour="#619CFF") +
  scale_shape_manual(values=c(2, 3, 19)) +
  ylim(c(0,1)) +
  xlab("Correlation between nodes 1 and 2") +
  ylab("Detection probability") +
  theme(legend.position="bottom")


ggsave("t1e.pdf", f1, width=4, height=4)
ggsave("cpower.pdf", f2, width=4, height=4)
ggsave("detect.pdf", f3, width=4, height=4)
