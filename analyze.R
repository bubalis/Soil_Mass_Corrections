options(warn = -1) 
source('ESM_and_mass_correction_functions.R')
#library(constrKriging)



site_metadata <- read.csv('Source_data/site_metadata.csv')



get_max_depths <- function(data){
  return(data$sample_depths[which.max(nchar(data$sample_depths))])
}

sigmoid <- function(x){
  return(1/(1+ exp(x*-1)))
}

logit <- function(p){
  return(log(p/ (1-p)))
}

get_site_meta <- function(.study){
  
  meta <- site_metadata %>% filter(study == .study)
  if (nrow(meta) == 0){print(paste('Could not find metadata for', .study))}
  return(meta)
}


#' Read in data from the van-haden ESM paper
read.VanHaden.data <- function(){
  input_file_name <- "Source_Data/Example_datasets.xlsx"
  # Name of the sheet on the spreadsheet that contains the data
  input_file_sheet <- "a_temporal_paired"
  return(read.xlsx(input_file_name, sheet=input_file_sheet))
}



simulate_benchmark.kriged <- function(data, mass_of_quant, n_sims = 500){
  data <- calc.cumulative_masses(data)
  model <- kmMonotonic1D(c(0, data$Cum_Min_Soil_g_cm2/max(data$Cum_Min_Soil_g_cm2)),
                          c(0, data$Cum_SOC_g_cm2), covtype = 'matern3_2' )
  simulated <- simulate_process.kmMonotonic1D(model, n_sims, newdata = mass_of_quant / max(data$Cum_Min_Soil_g_cm2))
  
  return(simulated[1,])
}


sim_benchmark_empirical_errs <- function(data, yest, mass_of_quant, err_fact, nsims = 500){
  err_sd <- 1/sqrt(sum_weights_mass(data, mass_of_quant)) * err_fact
  
  return(yest + rnorm(nsims, 0, err_sd))
}

sim_benchmark_logit <- function( data, yest, mass_of_quant, nsims = 500, sd_logit = .15){
  lower_i <- data %>% filter(Cum_Min_Soil_g_cm2 < mass_of_quant) %>% 
    pull(Cum_Min_Soil_g_cm2) %>% which.max()
  upper_i <- lower_i + 1
  dif <- yest-data$Cum_SOC_g_cm2[lower_i]
  span <- data$Cum_SOC_g_cm2[upper_i] - data$Cum_SOC_g_cm2[lower_i]
  
  log_o <- logit(dif/span)
  
  
  
  return (sigmoid(rnorm(nsims, log_o, sd_logit))* span + data$Cum_SOC_g_cm2[lower_i])
  
  
}

sim_benchmark_proportionate_logit <- function(data, yest, mass_of_quant, nsims = 500, sd_logit = .55){
  lower_i <- data %>% filter(Cum_Min_Soil_g_cm2 < mass_of_quant) %>% 
    pull(Cum_Min_Soil_g_cm2) %>% which.max()
  upper_i <- lower_i + 1
  dif <- yest-data$Cum_SOC_g_cm2[lower_i]
  span <- data$Cum_SOC_g_cm2[upper_i] - data$Cum_SOC_g_cm2[lower_i]
  
  log_o <- logit(dif/span)
  distance_to_edge <- min(abs(c(yest-0, yest-1)))
  return (sigmoid(rnorm(nsims, log_o, sd_logit * distance_to_edge))* span + data$Cum_SOC_g_cm2[lower_i])
}




sim_all_benchmarks_empirical_errs <- function(FD_data, res_data, depth_of_est, nsims = 500,
                                              err_factor = .0025){
  out <- data.frame()
  FD_data <- calc.cumulative_masses(FD_data)
  max_depths <- res_data %>% filter(nchar(sample_depths) == max(nchar(res_data$sample_depths))) %>% 
    pull(sample_depths) %>% first()
  for (i in seq(1, nrow(distinct(res_data, ID, Ref_ID)))){
    
    .ID <- distinct(res_data, ID, Ref_ID)[i,] %>% pull(ID) %>% first()
    .Ref_ID <- distinct(res_data, ID, Ref_ID)[i,] %>% pull(Ref_ID) %>% first()
    sim_data <- filter(FD_data, ID == .ID) %>% calc.cumulative_masses()
    
    yest <- res_data %>% filter(sample_depths == max_depths &  ID == .ID & Ref_ID == .Ref_ID) %>% 
      pull(Cum_SOC_g_cm2) %>% first()
    mass_of_quant <- filter(FD_data, ID == .Ref_ID & Lower_cm == depth_of_est) %>% 
      pull(Cum_Min_Soil_g_cm2) %>% first()
    
    simmed <- sim_benchmark_empirical_errs(sim_data, yest, mass_of_quant, err_fact = err_factor,
                                           nsims = nsims)
    
  
    out <- rbind(out, data.frame(SOC_benchmark = simmed, ID = .ID, Ref_ID = .Ref_ID,
                                 simulation.N = seq(1:nsims)))
    
  }
  return (out)
  }

sim_all_benchmarks_logit<- function(FD_data, res_data, depth_of_est, nsims = 500, sd_logit =.15){
  out <- data.frame()
  FD_data <- calc.cumulative_masses(FD_data)
  max_depths <- res_data %>% filter(nchar(sample_depths) == max(nchar(res_data$sample_depths))) %>% 
    pull(sample_depths) %>% first()
  for (i in seq(1, nrow(distinct(res_data, ID, Ref_ID)))){
     
    .ID <- distinct(res_data, ID, Ref_ID)[i,] %>% pull(ID) %>% first()
    .Ref_ID <- distinct(res_data, ID, Ref_ID)[i,] %>% pull(Ref_ID) %>% first()
    sim_data <- filter(FD_data, ID == .ID) %>% calc.cumulative_masses()
    
    yest <- res_data %>% filter(sample_depths == max_depths &  ID == .ID & Ref_ID == .Ref_ID) %>% 
      pull(Cum_SOC_g_cm2) %>% first()
    mass_of_quant <- filter(FD_data, ID == .Ref_ID & Lower_cm == depth_of_est) %>% 
      pull(Cum_Min_Soil_g_cm2) %>% first()
    
    
    #simmed <- sim_benchmark_logit(sim_data, yest, mass_of_quant, nsims, sd_logit )
    simmed <- sim_benchmark_proportionate_logit(sim_data, yest, mass_of_quant, nsims)
    out <- rbind(out, data.frame(SOC_benchmark = simmed, ID = .ID, Ref_ID = .Ref_ID,
                                 simulation.N = seq(1:nsims)))
    }
  return (out)
}


sim_all_benchmarks_kriged <- function(FD_data, res_data, depth_of_est, n_sims =500){
  out <- data.frame()
  d <- distinct(res_data, ID, Ref_ID) 
  for (i in seq(1, nrow(distinct(res_data, ID, Ref_ID)))){
      
     .ID <- d[i,] %>% pull(ID)
     .Ref_ID <- d[i,] %>% pull(Ref_ID)
     print(.ID)
     mass_of_quant <- FD_data %>% filter(ID == .Ref_ID) %>% calc.cumulative_masses() %>%
       filter(Lower_cm == depth_of_est) %>% pull(Cum_Min_Soil_g_cm2)%>% first()
     sim_data <- filter(FD_data, ID == .ID)
     simmed <- simulate_benchmark.kriged(sim_data, mass_of_quant, n_sims)
     out <- rbind(out, data.frame(SOC_benchmark = simmed, ID = .ID, Ref_ID = .Ref_ID,
                                  simulation.N = seq(1:n_sims)))
     
  }
  return (out)
  
}

#'compare the ESM and Mass-correction methods on a dataset.
#'@param fd_data a dataframe: data in the format used by VanHaden
#'@param ESM_depths a list of vectors: depth combinations to use the ESM method on
#'@param MC_depths a list of vectors: depth combinations to use the Mass-correction methods on.
#'@param depth_of_estimate: numeric: the depth we are trying get an accurate estimate of SOC for.
compare.methods <- function(fd_data,  ESM_depths, MC_depths, split_core_single_SOC_depths = c(),
                            depth_of_estimate = 30, lat = NULL, lon = NULL,
                            muname = NULL){
  #fd_data <- fd_data %>% select(-c('study'))
  comparisons <- fd_data %>% select(ID, Rep, Ref_ID) %>% distinct()
  
  #fixed depth results
  single_depth_agged <- df_agg_to_depths(fd_data, c(depth_of_estimate))
  FD_res <-  single_depth_agged  %>% 
    mutate(Cum_SOC_g_cm2 = Lower_cm * BD_g_cm3 * SOC_pct/100,
           method = 'Single Depth',
           sample_depths = as.character(depth_of_estimate)
    )
  
  mass_corr_avg_res <- single_depth_agged %>% 
    group_run_MC(adjustment_factor = 'standard_correction',  quantification_depth = depth_of_estimate) %>% 
    mutate(method = "Mass Correction, SOC_corrected?", 
           sample_depths = as.character(depth_of_estimate))
  
  mass_corr_fowler_weights <- single_depth_agged %>% 
    group_run_MC(adjustment_factor = 17/32,  quantification_depth = depth_of_estimate) %>% 
    mutate(method = "Mass Correction, Fowler Weights", 
           sample_depths = as.character(depth_of_estimate))
  
  all_res <- rbind(FD_res, mass_corr_avg_res, mass_corr_fowler_weights)
  
  for (depth_vals in ESM_depths){
    print(depth_vals)
    sample <- df_agg_to_depths(fd_data %>% filter(ID != Ref_ID), depth_vals)
    initial <- df_agg_to_depths(fd_data %>% filter(ID == Ref_ID), c(depth_of_estimate))
    data <- rbind(sample, initial)
    ESM_res <- data %>% 
      calc_ESM_VanHaden(T, depth_vals) %>% 
      filter(Lower_cm == depth_of_estimate) %>% #only keep the estimate for the depth we are calculating
      mutate(sample_depths = paste(depth_vals, collapse = ', '),
             method = 'ESM') %>% filter(ID != Ref_ID)
    all_res <- rbind(all_res, ESM_res)
    
    if (length(depth_vals) ==2){
      
      print('Running exponential decay ESM')
      print(data )
      ESM_decay_res <- run_decay_ESM(data, depth_of_estimate) %>% #only keep the estimate for the depth we are calculating
        mutate(sample_depths = paste(depth_vals, collapse = ', '),
               method = 'decay_ESM', Type = 'decay_ESM')
      
      
      all_res  <- rbind(all_res, ESM_decay_res)
      }}
  
  for (depth_vals in MC_depths){
    print(depth_vals)
    to_depths <-  rbind(df_agg_to_depths(fd_data %>% filter(ID != Ref_ID), depth_vals),
                        df_agg_to_depths(fd_data %>% filter(ID == Ref_ID), c(depth_of_estimate)))
    
    mass_correction_res <- to_depths %>% 
      group_run_MC(adjustment_factor = 1, quantification_depth = depth_of_estimate) %>% 
      mutate(method = "Mass Correction", 
             sample_depths = paste(depth_vals, collapse = ', '))
    
    all_res <- rbind(all_res, mass_correction_res)
  }
  for (depth_vals in split_core_single_SOC_depths){
    res <- group_run_joined_sample(fd_data, depth_vals, depth_of_estimate) %>%
      mutate(method = 'ESM, Single OC test', sample_depths = paste(depth_vals, collapse = ', '))
    all_res <- rbind(all_res, res)
    
  }
  
  if (!all(is.na(c(lat, lon, muname)))){
    ssurgo_res <- run_SSurgo_Mass_Corr(lat = lat, lon = lon,
                                       muname = muname, data = fd_data,
                                       depth_of_estimate = depth_of_estimate)%>% 
                              mutate(sample_depths = as.character(sample_depths))
    
    all_res <- rbind(all_res, ssurgo_res)
  }
  
  all_res <-all_res %>% ungroup() %>%select(ID, Rep, Upper_cm, Lower_cm, 
                                            Cum_SOC_g_cm2,
                                            method, sample_depths) %>% merge(comparisons)
  
  baseline_data <- all_res %>%ungroup() %>% filter( method == 'Single Depth' & Ref_ID == ID) %>% 
    rename( Cum_SOC_g_cm2_baseline = Cum_SOC_g_cm2) %>%
    select(Ref_ID, Rep, Cum_SOC_g_cm2_baseline)
  
  res_changes <- all_res %>% filter(ID != Ref_ID) %>% merge(baseline_data) %>%
    mutate(soc_change = Cum_SOC_g_cm2 - Cum_SOC_g_cm2_baseline) %>% 
    select(ID, Rep, Ref_ID, Cum_SOC_g_cm2, method, sample_depths, Cum_SOC_g_cm2_baseline,
           soc_change)
  
  return (res_changes)
}


mut_Ref_ID <- function(id){
  return (paste('Y5', substring(id, 4), sep = '_'))
}


#'Calculate the mass change to depth of quantification
mass_change <- function(FD_data, quant_depth){
  cum_FD <- calc.cumulative_masses(FD_data) 
  
  return (cum_FD %>% filter(ID != Ref_ID & Lower_cm == quant_depth) %>% 
    merge(cum_FD %>% filter(Lower_cm == quant_depth) %>% select(ID, Cum_Min_Soil_g_cm2), 
          by.x = 'Ref_ID', by.y = 'ID') %>% 
    mutate(mass.change = Cum_Min_Soil_g_cm2.x -  Cum_Min_Soil_g_cm2.y) %>% 
    select(ID, mass.change))
  
}

comparison_summarize_sim <- function(res, FD_data, .ID, depth_of_estimate,
                                 simulation_function = sim_all_benchmarks_empirical_errs){
  res <- res %>% arrange(ID) %>% 
    tidyr::fill(Cum_SOC_g_cm2_baseline, .direction = 'down') %>% 
    mutate(soc_change = Cum_SOC_g_cm2 - Cum_SOC_g_cm2_baseline, Ref = .ID )
  
  simulated_benchmarks <- simulation_function( FD_data, res, depth_of_estimate)
  res <- res %>% merge(simulated_benchmarks, by = c('ID', 'Ref_ID')) %>% 
    mutate(dif_from_benchmark = Cum_SOC_g_cm2 - SOC_benchmark,
           soc_change_benchmark = SOC_benchmark- Cum_SOC_g_cm2_baseline) %>% 
    
    left_join(mass_change(FD_data, depth_of_estimate))    
  
  return(res)
}



comparison_summarize <- function(res, FD_data, .ID, depth_of_estimate){
  res <- res %>% arrange(ID) %>% 
    tidyr::fill(Cum_SOC_g_cm2_baseline, .direction = 'down') %>% 
    mutate(soc_change = Cum_SOC_g_cm2 - Cum_SOC_g_cm2_baseline, Ref = .ID )
  
  max_depths <- res %>% filter(nchar(sample_depths) == max(nchar(res$sample_depths))) %>% 
    pull(sample_depths) %>% first()
  res <- res %>% merge(res %>% filter(sample_depths == max_depths) %>%
                         select(ID, soc_change, Ref) %>% 
                         rename(soc_change_benchmark = soc_change)) %>% 
    mutate(dif_from_benchmark = soc_change - soc_change_benchmark) %>% 
    left_join(mass_change(FD_data, depth_of_estimate))          
  return(res)
}


run_comparison <- function(FD_data, 
                          ESM_depths, MC_depths,
                          split_core_single_SOC_depths = c(),
                          depth_of_estimate = 30, 
                          study_name, .site = NULL){
  
  site_meta <- get_site_meta(study_name)
  all_out <- data.frame()
  
  #if a specific site is not specified, and each comparison is already specified
  if (is.null(.site) & ('Ref_ID' %in% colnames(FD_data))){
      lat <- site_meta %>% pull(lat) %>% first()
      lon <- site_meta %>% pull(lon) %>% first()
      muname <- site_meta %>% pull(muname) %>% first()
      
      res <- compare.methods(FD_data, ESM_depths = ESM_depths, MC_depths = MC_depths, 
                             split_core_single_SOC_depths = split_core_single_SOC_depths,
                             lat = lat, lon = lon, muname = muname) %>% mutate(study_name = study_name)
    
      all_out <- comparison_summarize_sim(res, FD_data, 0, depth_of_estimate) %>% 
        mutate(weight = 1)
      
  }else{
    if ('site' %in% colnames(FD_data)){
      #if a site is listed in the dataframe, only run comparisons between plots at the same site
      
      for (site_name in unique(FD_data$site)){
        
        new_res <- run_comparison(FD_data %>% filter(site == site_name) %>% select(-c(site)),
                                  ESM_depths = ESM_depths,
                                  MC_depths = MC_depths,
                                  split_core_single_SOC_depths = split_core_single_SOC_depths,
                                  depth_of_estimate = depth_of_estimate,
                                  study_name = study_name, .site = site_name)$res
        all_out <- rbind(all_out, new_res)
      }
      
     
    }else{
      #if a site argument was passed, get location data for specific site
      if (!is.null(.site)){
      site_dat <- site_meta %>% filter(site == .site)}
      
      else{site_dat <- site_meta}
      
      lat <- site_dat %>% pull(lat) %>% first()
      lon <- site_dat %>% pull(lon) %>% first()
      muname <- site_dat %>% pull(muname) %>% first()
  
      
  for (.ID in unique(FD_data$ID)){ 
    sub_fd <- FD_data %>% mutate(Ref_ID = .ID)
    res <- compare.methods(sub_fd, 
                           ESM_depths = ESM_depths,
                          MC_depths = MC_depths, 
                          split_core_single_SOC_depths = split_core_single_SOC_depths,
                          depth_of_estimate = depth_of_estimate, 
                          lat = lat, lon = lon, muname = muname) %>%
      mutate(study_name = study_name) %>%
      comparison_summarize_sim(sub_fd, .ID, depth_of_estimate) %>% 
      mutate(weight = 1/length(unique(FD_data$ID)))
    
    all_out <- rbind(all_out, res)
    
    
  }
      
      
      }
  }
  
  summary<- all_out %>% group_by(method, sample_depths) %>% 
    summarize(RMSE_from_benchmark = sqrt(mean(dif_from_benchmark**2)), 
              Bias.from.benchmark = mean(dif_from_benchmark * ifelse(mass.change >0, 1, -1)))
  return (list(summary = summary, res = all_out))
}


#Run comparisons
filtered_FD <- read.VanHaden.data() %>% preprocess.fd.data()

vanHoden_res <- compare.methods(filtered_FD, 
                                ESM_depths = list( c(10, 30), c(10, 30, 50, 100),  c(30, 50),
                                                  c(10, 30, 50)),
                                MC_depths = list(c(10,30), c(30)))


#run the vanHoden data backwards, so that the "after" soil is the 
#reference for the "before" soil
reversed <- filtered_FD %>% mutate (Ref_ID = sapply(Ref_ID, mut_Ref_ID))                            
vanHoden_rev <- compare.methods(reversed, 
                                ESM_depths = list(c(10, 30, 50, 100),  c(10, 30), c(30, 50),
                                                  c(10, 30, 50)),
                                MC_depths = list(c(10,30), c(30)))

vanHoden_all <- rbind(vanHoden_res, vanHoden_rev)

vanHoden_all <- vanHoden_all %>%  merge(vanHoden_all %>% filter(method == 'Single Depth') 
                                        %>% select(ID, Rep, soc_change) 
                                        %>% rename(FD_error = soc_change)) %>% 
  mutate(relative_error = soc_change / abs(FD_error) ) %>% select( -c(FD_error)) %>% 
  mutate(error = soc_change)

vanHoden_all %>% group_by(method, sample_depths) %>% summarize(RMSE = sqrt(mean(soc_change**2)))


#compare the results of different methods on the fowler simulated data, and 
#variations on it
source("data_sim.R")
fowler_sim_data <- simulate.soil.Fowler(soc_profile)


t_dt <- fowler_sim_data$t_dt
soc_delt_actual <-  fowler_sim_data$true_delta_soc

soc_delt_actual <- soc_delt_actual %>% 
  mutate(true_change_soc =true_change_soc/100,
         Ref_ID = ifelse(grepl('t0', ID) & grepl('t2', Ref_ID), paste(Ref_ID, 'x', sep = '_'),
                         Ref_ID)
  ) %>% filter(!(grepl('t1', ID) & grepl('t2', Ref_ID))) %>% filter(Ref_ID != ID) %>%
  select(ID, true_change_soc)


t_dt <- t_dt %>% filter(scenario != 's1')

y0_rev <- t_dt %>% filter(period %in% c('t0', 't1')) %>%
  mutate(Ref_ID = paste(scenario, 't1', sep = '_'))

#reversed sample_dt, with t0 as the ending point, and t2 as the starting point
y0_rev2 <- t_dt %>% filter(period %in% c('t0', 't2')) %>%
  mutate(period = paste(period, 'x', sep = '_')) %>%
  mutate(Ref_ID = paste(scenario, 't2_x', sep = '_'),
         ID = paste(scenario, period, sep = '_'))


fowler.all <- data.frame()
for (sample_data in list(t_dt, y0_rev, y0_rev2))
   {
     ESM_depths <- list(c(10,20,30),  c(15,30), c(24,30), c(24,30, 36),
                        c(30, 36), c(15, 30,50), c(10, 30, 50), c(30, 50), 
                        c(10, 30, 50, 70), c(15, 30, 45, 60, 75))
     MC_depths <- list(c(10,20,30),  c(15,30), c(24,30), c(10, 30), c(30))
     fowler.sim     <- compare.methods(sample_data, ESM_depths, MC_depths)
     fowler.sim <- fowler.sim %>% 
            distinct(method, ID, sample_depths, .keep_all = T)  %>% merge(soc_delt_actual) %>%
            mutate(error = soc_change - true_change_soc)
     
     fowler.all <- rbind(fowler.all, fowler.sim)
     
}

fowler.all <- fowler.all %>% merge(fowler.all %>% filter(method == 'Single Depth') 
                                   %>% select(ID, Rep, error) 
                                   %>% rename(FD_error = error),
                                   by = c('ID', 'Rep')) %>% 
  mutate(relative_error = error / abs(FD_error) ) %>% select( -c(FD_error))

fowler.all <- fowler.all %>% distinct(ID, sample_depths, soc_change, .keep_all = TRUE)
fowler.all %>% group_by(method, sample_depths) %>% summarize(RMSE = sqrt(mean(error**2)))




time_series_comparison <- function(FD_data, ESM_depths, MC_depths,
                                   study_name, depth_of_estimate = 30){
  out <- run_comparison(FD_data, ESM_depths = ESM_depths,
                                MC_depths = MC_depths, study_name = study_name)
  
  
  res <- out$res
  summary <- out$summary
  
  
  #summarize the between-treatment comparisons
  #Compare the calculated differences to the benchmark (ESM on all available data)
  calced_difs <- res %>% filter(ID != 'all_t1') %>% group_by(method, sample_depths) %>% 
    reframe(changes = dif_matrix(soc_change))
  
  
  
  as.df <- as.data.frame(calced_difs$changes) %>% 
    mutate(method = calced_difs$method, sample_depths = calced_difs$sample_depths)
  calced_difs <- as.df
  
    max_depths <- as.df %>% 
    filter(nchar(sample_depths) == max(nchar(res$sample_depths))) %>% 
    pull(sample_depths) %>% first()
  
  benchmark_difs <- as.df %>% filter(sample_depths == max_depths) %>% 
    select(-c(method, sample_depths))
  
  
  differences_summary <- data.frame()
  
  dist <- distinct(calced_difs, method, sample_depths)
  for ( i in seq(1: nrow(dist))){
    .sample_depths <- dist[i, 'sample_depths']
    .method <- dist[i, 'method']
    subset <- calced_difs %>% filter(.sample_depths == sample_depths, .method == method
    )
    difs <- select(subset, -c(method, sample_depths)) - benchmark_difs
    perc_difs <- difs / benchmark_difs
    RMSE <- sqrt(mean(colMeans(difs **2, na.rm = T)))
    RMSE_perc <- sqrt(mean(colMeans(perc_difs **2, na.rm = T)))
    mean_absolute_perc_err <- mean(colMeans(abs(perc_difs), na.rm = T))
    mean_absolute_err <- mean(colMeans(abs(difs), na.rm = T))
    .bias <- mean(difs[upper.tri(difs)])
    
    differences_summary <- rbind(differences_summary,
                                 data.frame(method = .method, sample_depths = .sample_depths,
                                            RMSE_differences = RMSE, RMSE_differences_perc = RMSE_perc,
                                            mean_absolute_perc_err_dif = mean_absolute_perc_err,
                                            difs_err = mean_absolute_err,
                                            bias = .bias))
  }
  summary <- merge(summary, differences_summary, by = c('method', 'sample_depths') )
  
  return (list('summary' = summary, 'res' = res %>% mutate('study' = study_name)))}




sainju_FD <- load_study('Sainju_2013')
muname <- 'Dooley sandy loam, 0 to 4 percent slopes'
ESM_depths <- list(c(15, 30), c(15,30, 60), c(15, 30, 60, 90), c(30, 60), c(7.5, 15, 30, 60),
                  c(7.5, 15, 30, 60, 90, 120), c(7.5, 15,30))
MC_depths <- list(c(15, 30), c(30), c(7.5, 15, 30))
single_core_depths = list(c(15,30), c(15,60), c(30,60))

sainju_summary <- run_comparison(sainju_FD, ESM_depths = ESM_depths,
                                 MC_depths = MC_depths, 
                                 split_core_single_SOC_depths = single_core_depths,
                                 study_name = 'Sainju_2013'
                                            )


mishra_fd <- load_study('mishra_2010')
study_name = 'mishra_2010'
mishra_ESM <- list(c(10, 20, 30, 40), c(20, 30), c(10,20,30),
                    c(30, 40), c(20,30, 40))
mishra_MC <- list(c(20, 30), c(10, 20, 30), c(10, 30), c(30))
mishra_splitcore <- list(c(20,30), c(30,40) )
mishra_fd <- mishra_fd %>% mutate(ID = paste(ID, site, sep = '_'))

mishra_res <- run_comparison(mishra_fd,
        ESM_depths = mishra_ESM, MC_depths = mishra_MC, 
        split_core_single_SOC_depths = mishra_splitcore,
                                    study_name = study_name)

vanDoren_fd <- load_study('Van_doren_1986')

vanDoren_sp_4_t <- run_comparison(vanDoren_fd %>% filter(!grepl('t1', ID)) %>% select(-c(Ref_ID)), 
                                           ESM_depths = list(c(30,40), c(20,30),
                                                             c(20,30,40), c(30, 45),
                                                             c(10, 20,30),
                                                             c(10, 30),
                                                             c(30,35),
                                                             c(15,30), c(25, 30),
                                                             c(25,30, 35), c(20, 30, 45),
                                                             c(15, 30, 45),
                                                             c(5, 10, 15, 20, 25, 30, 35, 40, 45)),
                                           MC_depths = list(c(25,30), c(20, 30), c(15, 30),
                                                            c(10, 30), c(5,30), c(30)), 
                                  split_core_single_SOC_depths = list(c(20, 40), c(25,35), c(20,30), c(20, 45),
                                                                     c(25, 30), c(15,30), c(20, 35), 
                                                                     c(30,45)),
                                  study_name = 'Van_doren_1986')


devine_fd <- load_study('Devine_2014')

devine_res <-run_comparison(devine_fd, ESM_depths = list(c(5, 15,30, 50, 100),
                                                        c(5, 15, 30), 
                                                          c(15, 30), c(15,30, 50),
                                        c(15, 30, 50), c(15, 30, 50, 100), c(30,50)),
                                        MC_depths = list(c(5, 15, 30), c(30), c(15, 30)),
                            split_core_single_SOC_depths = list(c(15,30), c(15, 50), c(30,50)),
                                        study_name = 'Devine_2014'
                                        )

blanco_fd <- load_study('Blanco_Canqui(2008)')
blanco_res <- run_comparison(blanco_fd, 
                             ESM_depths = list(c(5, 10, 30), c(10, 30),  c(30, 50),
                                               c(30,60),
                                                   c(10, 30, 50), c(10, 30, 50, 60),
                                                   c(5, 10, 30, 50, 60)
                                                   ),
                      MC_depths = list(c(30), c(10, 30)),
                      split_core_single_MC_depths = list(c(10,30), c(30,50), c(10, 50)),
                     study_name = 'Blanco_Canqui(2008)')
  
  




venterea_fd <- load_study('Venterea_2006')

venterea_s4_t1 <- run_comparison(venterea_fd %>% filter(years == 2000) %>% select(-c(Ref_ID)),
                                           ESM_depths = list(c(5, 10, 20, 30, 45, 60),
                                                             c(30, 45), c(20, 30),
                                                             c(20, 30, 45), c(30, 60),
                                                             c(10, 30, 45), c(5, 10, 30),
                                                             c(10, 20, 30), c(5, 10, 20, 30),
                                                             c(5, 10, 20, 30, 45),
                                                             c(5, 10, 20, 30, 45, 60)
                                            
                                                             ),
                                           MC_depths = list(c(30), c(20, 30), c(10, 20, 30), 
                                                            c(5, 10, 20, 30)),
                                 split_core_single_SOC_depths = list(c(20,30), 
                                                                     c(20, 45),
                                                                     c(30, 45),
                                                                     ),
                                           study_name = 'Venterea_2006'
                                           )

venterea_s4_t2 <- run_comparison(venterea_fd %>% filter(years == 2005) %>% select(-c(Ref_ID)),
                                 ESM_depths = list(c(5, 10, 20, 30, 45, 60),
                                                   c(30, 45), c(20, 30),
                                                   c(20, 30, 45), c(30, 60),
                                                   c(10, 30, 45), c(5, 10, 30),
                                                   c(10, 20, 30), c(5, 10, 20, 30),
                                                   c(5, 10, 20, 30, 45),
                                                   c(5, 10, 20, 30, 45, 60)
                                                   
                                 ),
                                 MC_depths = list(c(30), c(20, 30), c(10, 20, 30), 
                                                  c(5, 10, 20, 30)),
                                 split_core_single_SOC_depths = list(c(20,30), c(20, 45),
                                                                     c(30,45)),
                                 study_name = 'Venterea_2006'
)



chatterjee_fd <- load_study('chatterjee_2009')
chj_res <- run_comparison(chatterjee_fd, 
                                ESM_depths = list(c(5, 10, 30), c(10, 30),  c(30, 50),
                                                  c(30,60),
                                                   c(10, 30, 50), c(10, 30, 50, 60),
                                                   c(5, 10, 30, 50, 60)),
                                   MC_depths = list(c(30), c(10, 30), c(5,10,30)),
                          split_core_single_SOC_depths = list(c(30,50), c(10,50)),
                                   study_name = 'chatterjee_2009'
                                   )
    

blanco_fd2 <- load_study('Blanco-canqui_2011')
blanc_res_2 <- run_comparison(blanco_fd2, ESM_depths = list(c(30,40), c(20,30), c(20,30, 40), c(30,60),
                                                   c(20,30,60), c(20, 30,40, 60),
                        
                                                   c(15, 20, 30,40, 60),
                                                   c(15,30), c(10, 30),
                                                   c(20, 30, 40, 60, 80),
                                                   c(15, 20, 30, 40, 60, 80),
                                                   c(2.5, 5, 10, 15, 20, 30, 40, 60, 80, 100)
  ),
  MC_depths = list(c(30), c(10, 30), c(20, 30), c(15,30)),
  split_core_single_SOC_depths = list(c(20,40), c(15, 40), c(20,30), c(30,40)),
   study = 'Blanco_Canqui_2011'
  )
  
  
poffenbarger_fd <- load_study('poffenbarger_2020')
poffenbarger_res <- run_comparison(poffenbarger_fd, 
                        ESM_depths = list(c(15, 30), 
                        c(30, 60), c(15,30, 60), c(15,30,60, 90) ),
                         MC_depths = list(c(15, 30), c(30)),
                        split_core_single_SOC_depths = list(c(15,30), c(15,60), c(30,60)),
                        study_name = 'poffenbarger_2020')
                         

yang_FD <- load_study('Yang_1999')    
yang_res <- run_comparison(yang_FD, 
                           ESM_depths = list(c(20,30), c(30,40), c(30,50),
                                             c(20,30, 40), c(20,30,40, 50),
                                             c(5, 20,30,40), c(5,20,30,50),
                                             c(5,20, 30,40,50,70,90),
                                             c(10,30)
                                             ),
                           MC_depths = list(c(5,30), c(10,30), c(20,30), c(30)),
                           split_core_single_SOC_depths = list(c(20, 40), c(20, 50), c(20,30),
                                                               c(30,40), c(30,50)),
                           study_name = 'Yang_1999')

out_s4t <- rbind(mishra_res$res, blanco_res$res, 
      devine_res$res, vanDoren_sp_4_t$res,
      sainju_summary$res, 
      
      venterea_s4_t1$res,
      venterea_s4_t2$res,
      chj_res$res, blanc_res_2$res,
      poffenbarger_res$res,
      yang_res$res)


summarize_res <- function(grouped_res){
  grouped_res %>% summarise(RMSE = sqrt(weighted.mean(dif_from_benchmark**2, weight)), 
                            RMPSE = sqrt(weighted.mean((dif_from_benchmark /SOC_benchmark )**2, weight )), 
                            bias = weighted.mean(dif_from_benchmark * ifelse(mass.change >0, 1, -1), weight), 
                            MAPE_change = weighted.mean(abs(dif_from_benchmark/soc_change_benchmark), weight),
                            Max_err =  max(abs(dif_from_benchmark)),
                            mean_abs_err = mean(abs(dif_from_benchmark)))
}

out_s4t <- out_s4t %>% distinct(ID, Ref_ID, method, sample_depths, simulation.N, .keep_all = T)

ag_res <- out_s4t %>% filter(!grepl('WL', ID) & ID!= "Forest Succession" & !grepl('WL', Ref_ID) & Ref_ID != "Forest Succession")

ag_res %>%  group_by(method, sample_depths) %>% summarize_res() %>% 
  filter(sample_depths %in% c('30, 40', '30, 45', '30, 50', '30, 60', '30' ))



ggplot(ag_res %>% filter(sample_depths %in% c('30, 40', '30, 45', '30, 50', '30, 60')), 
       aes(mass.change, dif_from_benchmark, color = as.factor(method))) + 
  geom_violin(alpha = .25) + geom_hline(yintercept  = 0) + 
  facet_wrap(~sample_depths)

ggplot(ag_res %>% filter(sample_depths %in% c('10, 30', '20, 30', '15, 30')), 
       aes(mass.change, dif_from_benchmark, color = as.factor(method))) + 
  geom_point(alpha = .25) + geom_hline(yintercept  = 0) + 
  facet_wrap(~sample_depths)


s<- ag_res %>% filter(sample_depths == '30, 45') %>% pull(study_name) %>% unique()
dat_to_plot <- ag_res %>% filter(
  study_name %in% s & ((sample_depths == '30, 45' & grepl('ESM', method)) | 
  (sample_depths == '30' & 
     method %in% c('Mass Correction', 'Mass Correction, SSurgo avg', 'Single Depth' ))))

ggplot(dat_to_plot, aes(x = mass.change, y = dif_from_benchmark*100, color = as.factor(method))) + geom_point(alpha = .5) + 
  geom_hline(yintercept = 0)


cross_fold_weights_res <- data.frame()
for (study in unique(out_s4t$study_name)){
  training_data <- out_s4t %>% filter(study_name != study)
  sumed_dat <- training_data %>% filter(!grepl("WL", Ref) &  !grepl('WL', ID) & !grepl('Forest', ID) & !grepl('Forest', Ref)) %>%
    filter(sample_depths == 30) %>% group_by(method) %>% 
    
    summarise(bias = weighted.mean(dif_from_benchmark * ifelse(mass.change >0, 1, -1), weight))
  
  single_d_bias <- sumed_dat %>% filter(method == 'Single Depth') %>% 
    pull('bias') %>% first()
  
  linear_MC_bias <- sumed_dat %>% filter(method == 'Mass Correction') %>% 
    pull('bias') %>% first()
   MC_weight <- abs(single_d_bias)/(abs(single_d_bias) + abs(linear_MC_bias))
   sd_weight <- 1-MC_weight
   testing_data <- out_s4t %>% filter(study == study_name) %>% 
     filter(!grepl("WL", Ref) &  !grepl('WL', ID) & !grepl('Forest', ID) & !grepl('Forest', Ref))
   testing_data %>% filter(method %in% c('Single Depth', 'Mass Correction') & sample_depths == '30') %>% 
     pivot_wider(id_cols = c(ID, Ref), names_from = method, values_from = soc_change) %>% 
     mutate(soc_change = `Mass Correction` * MC_weight + `Single Depth` * sd_weight ) %>% 
     merge(testing_data %>% select(ID, Ref, soc_change_benchmark, mass.change, weight) %>% distinct(ID, Ref, .keep_all = T)) %>%
     mutate(dif_from_benchmark = soc_change - soc_change_benchmark, study = study_name,
            method = 'Mass Correction, empirical weights') %>%
     rbind(cross_fold_weights_res ) -> cross_fold_weights_res
   
}








### calculate time-series differences on real data
ESM_depths <- list(c(5, 10, 15,20, 25, 30,35, 40, 45), c(30, 45), c(25, 30, 35), c(15,30, 45), 
                   c(20, 30, 40), c(30, 40),
                   c(10, 20, 30), c(15, 30), c(20, 30))
MC_depths <- list(c(30), c(15,30), c(20, 30), c(25, 30), c(10, 20, 30))
muname <- 'Wooster silt loam, 2 to 6 percent slopes'

vanDoren_res <- time_series_comparison(vanDoren_fd, ESM_depths, MC_depths, 
                                       study_name = 'van_doren_1986')




#sanford_res <- time_series_comparison(sanford_fd, 
#                                      ESM_depths = list(c(15, 30), c(15,30, 60), 
#                                                   c(15, 30, 60, 90), c(30, 60)),
#                                      MC_depths = list(c(15, 30), c(30)), 
#                                      study_name = 'Sanford_2012'
#)

venterea_res <- time_series_comparison(venterea_fd,  
                                    ESM_depths = list(
                                      c(30, 45), c(20, 30),
                                               c(20, 30, 45), c(30, 60),
                                               c(10, 30, 45), c(5, 10, 30),
                                               c(10, 20, 30), c(5, 10, 20, 30),
                                               c(5, 10, 20, 30, 45),
                                               c(5, 10, 20, 30, 45, 60)),
                          MC_depths = list(c(30), c(20, 30), c(10, 20, 30), 
                 c(5, 10, 20, 30), c(30)),
                
                 
                 study_name = 'Venterea_2006.csv')



FD_data <- vanDoren_fd %>% filter(grepl('t2', ID))

id <- "COM_NT"

dat <- rbind(df_agg_to_depths(FD_data, c(5,10,15,20,25,35,40,45)) %>% 
               mutate(Ref_ID = paste(ID, 'to30', sep = '_')), 
             df_agg_to_depths(FD_data, c(30)) %>% mutate(ID = paste(ID, 'to30', sep = '_'))
                                                        )


powerset = function(s){
  len = length(s)
  l = vector(mode="list",length=2^len) ; l[[1]]=numeric()
  counter = 1L
  for(x in 1L:length(s)){
    for(subset in 1L:counter){
      counter=counter+1L
      l[[counter]] = c(l[[subset]],s[x])
    }
  }
  return(l)
}

errors_empirical <- function(FD_data, depth_to_predict, spline_depths){
  dat <- rbind(df_agg_to_depths(FD_data, spline_depths) %>% 
          mutate(Ref_ID = paste(ID, 'to30', sep = '_')), 
        df_agg_to_depths(FD_data, c(depth_to_predict)) %>% mutate(ID = paste(ID, 'to30', sep = '_')) 
        %>% mutate(Ref_ID = ID)
        
  )
  cumulative_ESM <- calc_ESM_VanHaden(dat, T, c(depth_to_predict))
  out <- cumulative_ESM %>% merge(dat%>% calc.cumulative_masses() %>% ungroup() %>% filter(Ref_ID == ID) %>% 
                                   select(Ref_ID, Cum_SOC_g_cm2), by = 'Ref_ID')
  return(out$Cum_SOC_g_cm2.x - out$Cum_SOC_g_cm2.y)
}



all_errors_empirical <- function(FD_data, depth_to_predict){
  depth_vals <- FD_data$Lower_cm %>% unique() %>% setdiff(c(depth_to_predict))
  all_combos <- powerset(depth_vals)
  out <- data.frame()
  for (v_set in all_combos[2: length(all_combos)]){
    errs <- errors_empirical(FD_data, depth_to_predict, v_set)
    res <- data.frame(matrix(nrow = length(errs), ncol = length(depth_vals)))
    colnames(res) <- as.character(depth_vals)
    res[as.character(v_set)] = 1
    res <- mutate(res, err = errs, ID = unique(FD_data$ID))
    out <- rbind(out, res)
  }
  out[is.na(out)] <-0
  return(out)
}
