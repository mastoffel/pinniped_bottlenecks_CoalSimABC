#' loads multiple data frames into list from xlsx and assigns the 
#' tab names from the spreadsheet to the list elements
#' 
#' @param data_path path to excel file with multiple data.frames 


load_seal_data <- function(data_path) {
  ### load all data
  library(readxl)
  # sheet numbers to load
  data_path <- "../data/seal_data_largest_clust_and_pop.xlsx"
  dataset_names <- excel_sheets(data_path)
  
  load_dataset <- function(dataset_names) {
    read_excel(data_path, sheet = dataset_names)
  }
  # load all datasets
  all_seals <- lapply(dataset_names, load_dataset)
  names(all_seals) <- dataset_names
  all_seals
}


#' creates the tbs file for using priors
#' @param N_sim number of values to draw for all parameters
#' @param N_pop population size
#' @param N_samp sample size
#' @param N_loc number of loci
#' @param model so far c("bottleneck", "neutral", "decline")



create_tbs_file <- function(N_sim, N_pop, N_samp, N_loc, model = c("bottleneck", "neutral", "decline")) {
  ## diploid pop size
  N0 <- round(runif(1, min = N_pop / 20, max = N_pop / 10), 0)
  # N0 <- N_pop
  
  ## mutation rate
  mu <- runif(1, min = 0.0005, max = 0.005)
  # mu <- 0.0005
  
  ## theta
  theta <- 4 * N0 * mu
  
  #### bottleneck end ####
  # time is always thought in GENERATIONS
  latest_end <- 1 # generations ago
  earliest_end <- 20 # generations ago
  # see ms manual. time is expressed in units of 4*N0
  end_bot <- runif(1, min = latest_end / (4*N0), max = earliest_end / (4*N0))
  
  #### bottleneck start ####
  latest_start <- 20 # generations ago
  earliest_start <- 100 # generations ago
  # see ms manual. time is expressed in units of 4*N0
  start_bot <- runif(1, min = latest_start / (4*N0), max = earliest_start / (4*N0))
  
  ## bottleneck population size 1 - 1000, expressed relative to N0
  N_bot <- round(runif(1, min = 1 / N0, max = 1000 / N0), 4)
  
  ## historical populaiton size 1 - 100 times as big as current
  N_hist <- round(runif(1, min = 1 , max = 1000), 0)
  
  if (model == "bottleneck") {
    out <- data.frame(theta, end_bot, N_bot, start_bot, N_hist)
  }
  
  if (model == "bottleneck") {
    out <- data.frame(theta, start_bot, N_hist)
  }
  
  if (model == "constant") {
    out <- data.frame(theta)
  }

  out
  
}





