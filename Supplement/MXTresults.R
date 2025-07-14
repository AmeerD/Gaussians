library(dplyr)
library(tidyr)
library(purrr)

resdir <- "./res/"

resfiles <- list.files(resdir)

loadres <- function(fpath, cols) {
  res <- read.delim(fpath, sep=" ", header=F)
  colnames(res) <- cols
  res <- res %>% mutate(sim = strsplit(gsub(resdir, "", fpath), "_")[[1]][1])
  return(res)
}

fullres <- map(resfiles[grepl("^t", resfiles)], 
               ~loadres(paste(resdir, .x, sep=""),
                        c("df", "idx", "row", "col", "LRStat"))) %>%
  list_rbind()

save(fullres, file="mxtresults.rda")