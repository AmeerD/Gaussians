library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(latex2exp)

load("Ameer's Notes/Data Section - August 2024/mxnresults.rda")
load("Ameer's Notes/Data Section - August 2024/mxnpowerresults.rda")

nullres <- fullres %>%
  filter(grepl("null", sim)) 

f1 <- nullres %>%
  mutate(LRStat = pchisq(as.numeric(LRStat), df=1, lower.tail=F)) %>%
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
  xlim(c(0,1)) + ylim(c(0,1)) +
  xlab("Significance cutoff") + ylab("P-value quantiles")

f1

f2 <- powerres %>%
  mutate(LRStat = as.numeric(LRStat)) %>%
  filter(!is.na(LRStat)) %>% 
  mutate(rej = LRStat > qchisq(0.95, 1)) %>%
  group_by(cor, eps) %>%
  summarise(cpower = sum(rej)/n()) %>%
  mutate(epsilon=as.character(eps)) %>%
  ggplot(aes(x=cor, y=cpower, linetype=epsilon, shape=epsilon)) +
  geom_line(colour="#619CFF") +
  geom_point(colour="#619CFF") +
  scale_shape_manual(values=c(2, 3, 19), 
                     labels=unname(TeX(c("$0.25^{1/4}$", "$0.5^{1/4}$", "$0.75^{1/4}$")))) +
  scale_linetype_manual(values = c(1,3,5),
    labels=unname(TeX(c("$0.25^{1/4}$", "$0.5^{1/4}$", "$0.75^{1/4}$")))) +
  ylim(c(0,1)) + 
  xlab(unname(TeX(c("Correlation between nodes 1 and 2 ($\\omega$)")))) +
  ylab("Conditional power") +
  labs(shape=unname(TeX("$q_1$")), linetype=unname(TeX("$q_1$"))) +
  theme(legend.position="bottom")

f2

f3 <- powerres %>%
  mutate(detected = !is.na(as.numeric(LRStat))) %>%
  group_by(cor, eps) %>% 
  summarise(dprob = sum(detected)/n()) %>%
  mutate(epsilon=as.character(eps)) %>%
  ggplot(aes(x=cor, y=dprob, linetype=epsilon, shape=epsilon)) +
  geom_point(colour="#619CFF") +
  geom_line(colour="#619CFF") +
  scale_shape_manual(values=c(2, 3, 19), 
                     labels=unname(TeX(c("$0.25^{1/4}$", "$0.5^{1/4}$", "$0.75^{1/4}$")))) +
  scale_linetype_manual(values = c(1,3,5),
                        labels=unname(TeX(c("$0.25^{1/4}$", "$0.5^{1/4}$", "$0.75^{1/4}$")))) +
  ylim(c(0,1)) +
  xlab(unname(TeX(c("Correlation between nodes 1 and 2 ($\\omega$)")))) +
  ylab("Detection probability") +
  labs(shape=unname(TeX("$q_1$")), linetype=unname(TeX("$q_1$"))) +
  theme(legend.position="bottom")

f3


ggsave("Ameer's Notes/Data Section - August 2024/t1e.pdf", f1, width=4, height=4)
ggsave("Ameer's Notes/Data Section - August 2024/cpower.pdf", f2, width=4, height=4)
ggsave("Ameer's Notes/Data Section - August 2024/detect.pdf", f3, width=4, height=4)

