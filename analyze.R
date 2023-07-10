
source('ESM_and_mass_correction_functions.R')


#' Read in data from the van-haden ESM paper
read.VanHaden.data <- function(){
  input_file_name <- "Source_Data/Example_datasets.xlsx"
  # Name of the sheet on the spreadsheet that contains the data
  input_file_sheet <- "a_temporal_paired"
  return(read.xlsx(input_file_name, sheet=input_file_sheet))
}



#'compare the ESM and Mass-correction methods on a dataset.
#'@param fd_data a dataframe: data in the format used by VanHaden
#'@param ESM_depths a list of vectors: depth combinations to use the ESM method on
#'@param MC_depths a list of vectors: depth combinations to use the Mass-correction methods on.
#'@param depth_of_estimate: numeric: the depth we are trying get an accurate estimate of SOC for.
compare.methods <- function(fd_data,  ESM_depths, MC_depths,
                            depth_of_estimate = 30){
  
  
  
  comparisons <- fd_data %>% select(ID, Rep, Ref_ID) %>% distinct()
  
  #fixed depth results
  FD_res <- fd_data %>% df_agg_to_depths(c(depth_of_estimate)) %>% 
    mutate(Cum_SOC_g_cm2 = Lower_cm * BD_g_cm3 * SOC_pct/100,
           method = 'Single Depth',
           sample_depths = as.character(depth_of_estimate)
    )
  
  
  all_res <- FD_res
  
  
  for (depth_vals in ESM_depths){
    ESM_res <- fd_data %>% df_agg_to_depths(depth_vals) %>% 
      calc_ESM_VanHaden(T, depth_vals) %>% 
      filter(Lower_cm == depth_of_estimate) %>% #only keep the estimate for the depth we are calculating
      mutate(sample_depths = paste(depth_vals, collapse = ', '),
             method = 'ESM')
    all_res <- rbind(all_res, ESM_res)
  }
  
  for (depth_vals in MC_depths){
    to_depths <- fd_data %>% df_agg_to_depths(depth_vals)
    
    mass_correction_res <- to_depths %>% 
      group_run_MC(adjustment_factor = 1, quantification_depth = depth_of_estimate) %>% 
      mutate(method = "Mass Correction", 
             sample_depths = paste(depth_vals, collapse = ', '))
    
    mass_corr_unbiased_res <- to_depths %>% df_agg_to_depths(depth_vals) %>% 
      group_run_MC(adjustment_factor = 'standard_correction',  quantification_depth = depth_of_estimate) %>% 
      mutate(method = "Mass Correction, SOC_corrected?", 
             sample_depths = paste(depth_vals, collapse = ', '))
    
    all_res <- rbind(all_res, mass_correction_res, mass_corr_unbiased_res)
    
  }
  
  all_res <-all_res %>% ungroup() %>%select(ID, Rep, Upper_cm, Lower_cm, 
                                            Cum_SOC_g_cm2,
                                            method, sample_depths) %>% merge(comparisons)
  
  baseline_data <- all_res %>%ungroup() %>% filter( method == 'Single Depth' & Ref_ID == ID) %>% 
    rename( Cum_SOC_g_cm2_baseline = Cum_SOC_g_cm2) %>%
    select(Ref_ID, Rep, Cum_SOC_g_cm2_baseline)
  
  res_changes <- all_res %>% filter(ID != Ref_ID) %>% merge(baseline_data) %>%
    mutate(soc_change = Cum_SOC_g_cm2 - Cum_SOC_g_cm2_baseline) %>% 
    select(ID, Rep, Cum_SOC_g_cm2, method, sample_depths, Cum_SOC_g_cm2_baseline,
           soc_change)
  
  return (res_changes)
}


mut_Ref_ID <- function(id){
  return (paste('Y5', substring(id, 4), sep = '_'))
}










#Run comparisons

filtered_FD <- read.VanHaden.data() %>% preprocess.fd.data()
#filtered_FD <- filtered_FD %>% mutate(Ref_ID = sapply(ID, set_Ref_ID ))



vanHoden_res <- compare.methods(filtered_FD, 
                                ESM_depths = list(c(10, 30, 50, 100),  c(10, 30), c(30, 50),
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


#folwer_dt < - fowler_sim_data$fowler_dt
#sample_dt <- fowler_sim_data$sample_dt

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



#Data from the Wisconsin Integrated Cropping Systems Trial
sanford_fd <- read.csv('Source_Data/Sanford2012_data.csv') %>% mutate(Rep = 1)

sanford_res<- compare.methods(sanford_fd, ESM_depths = list(c(15, 30), c(15,30, 60), c(15, 30, 60, 90), c(30, 60)),
                              MC_depths = list(c(15, 30), c(30)))


MC_SSurgo <- run_SSurgo_Mass_Corr(lat = 43.29586, lon = -89.38068,
                     data = sanford_fd)

sanford_res <- sanford_res %>% rbind(MC_SSurgo) %>% arrange(ID) %>% 
  tidyr::fill(Cum_SOC_g_cm2_baseline, .direction = 'down') %>% 
  mutate(soc_change = Cum_SOC_g_cm2 - Cum_SOC_g_cm2_baseline )

#set the benchmark as the ESM method using all available depths
sanford_res <- sanford_res %>% merge(sanford_res %>% filter(sample_depths == '15, 30, 60, 90') 
                                     %>% select(ID, soc_change) %>% rename(soc_change_benchmark = soc_change)) %>% 
  mutate(dif_from_benchmark = soc_change - soc_change_benchmark)


summary<- sanford_res %>% group_by(method, sample_depths) %>% 
  summarize(RMSE_from_benchmark = sqrt(mean(dif_from_benchmark**2)), 
            Bias.from.benchmark = mean(dif_from_benchmark))



#summarize the between-treatment comparisons
#Compare the calculated differences to the benchmark (ESM on all available data)
calced_difs <- sanford_res %>% filter(ID != 'all_t1') %>% group_by(method, sample_depths) %>% 
  reframe(changes = dif_matrix(soc_change))

as.df <- as.data.frame(calced_difs$changes) %>% 
  mutate(method = calced_difs$method, sample_depths = calced_difs$sample_depths)
calced_difs <- as.df

benchmark_difs <- as.df %>% filter(sample_depths == '15, 30, 60, 90') %>% select(-c(method, sample_depths))


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

sanford_summary <- merge(summary, differences_summary )
