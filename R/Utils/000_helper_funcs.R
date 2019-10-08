#select the first plot of each experiment
select_first <- function(data){
  data <- data %>% filter(Plot_ID %in% c("FPWW0120001", "FPWW0180001", "FPWW0220001"))
}

#remove unmeasured replicates from design
remove_reps <- function(data, value, missing){
  remove <- ifelse(sum(is.na(data %>% pull(value)), na.rm = TRUE) > missing, TRUE, FALSE)
}

#remove unmeasured Experiments from design (wrapper for remove_reps)
##this function cleans datasets after match-joining model predictions and scoring data
remove_exps <- function(data, value, missing){
  data <- data %>% group_by(Exp) %>% nest() %>% 
    mutate(remove = purrr::map_chr(data, remove_reps, value = value, missing = missing)) %>% 
    filter(remove == FALSE) %>% dplyr::select(-remove) %>% unnest()
}

#extract number of unique Plot_IDs from datasets
check_fun <- function(data){
  length <- unique(data$Plot_ID) %>% length
}

#remove plots that were not measured for the whole duration of the experiment
remove_incomplete_series <- function(data) {
  complete_series <- data %>% filter(Lot != 6) %>% filter(!(Lot == 3 & Exp == "FPWW022"))
}
