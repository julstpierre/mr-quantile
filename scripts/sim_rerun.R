source("utils.R")
library(stringr)
library(dplyr)

sim_name <- c("main")
backup_dir <- paste0("../results/rds/", sim_name, "/backup")
all_files <- list.files(backup_dir)

tag_method <- function(res_list, method) {
  bind_rows(res_list) %>% 
    dplyr::mutate(method = method)
}

for (f in all_files){
  
  (scenario <- str_extract(f, "D_S\\d+"))
  (start_idx <- as.integer(str_match(f, "chunk_(\\d+)_to_")[,2]))
  (end_idx <- as.integer(str_match(f, "_to_(\\d+)_")[,2]))
  
  load(paste0("../data/simulations/", sim_name, "/", scenario, ".RData"))
  D <- get(scenario)
  chunk <- D[start_idx:end_idx]
  
  res_mr_quantreg <- 
    tag_method(lapply(chunk, fit_mr_quantreg), "mr_quantreg")
  
  res_all <-
    readRDS(paste0(backup_dir, "/", f)) %>%
    filter(method != "mr_quantreg") %>%
    bind_rows(res_mr_quantreg)
  
  saveRDS(res_all, paste0("../results/rds/", sim_name, "/", f))
}
