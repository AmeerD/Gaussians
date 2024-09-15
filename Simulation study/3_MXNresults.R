library(dplyr)
library(tidyr)
library(purrr)

resdir <- "./res/"

resfiles <- list.files(resdir)

loadres <- function(fpath, cols) {
  res <- read.delim(fpath, sep=" ", header=F)
  # res <- bind_cols(
  #   res[c(TRUE,FALSE),],
  #   data.frame(LRStat = res[c(FALSE,TRUE),1])
  # )
  colnames(res) <- cols
  # print(head(res))
  res <- res %>% mutate(sim = strsplit(gsub(resdir, "", fpath), "_")[[1]][1])
  return(res)
}

fullres <- map(resfiles[!grepl("power", resfiles)], 
               ~loadres(paste(resdir, .x, sep=""),
                        c("method", "idx", "row", "col", "LRStat"))) %>%
  list_rbind()

save(fullres, file="mxnresults.rda")

powerres <- map(resfiles[grepl("power", resfiles)], 
                ~loadres(paste(resdir, .x, sep=""),
                         c("cor", "eps", "idx", "row", "col", "LRStat"))) %>%
  list_rbind()

# powerresold <- map(resfiles[grepl("power25_", resfiles)], 
#                 ~loadres(paste(resdir, .x, sep=""),
#                          c("cor", "idx", "row", "col", "LRStat"))) %>%
#   list_rbind() %>% mutate(eps = 0.5)

save(powerres, file="mxnpowerresults.rda")