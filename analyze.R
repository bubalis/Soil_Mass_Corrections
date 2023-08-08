
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
    ESM_res <-  rbind(df_agg_to_depths(fd_data %>% filter(ID != Ref_ID), depth_vals),
                          df_agg_to_depths(fd_data %>% filter(ID == Ref_ID), c(depth_of_estimate))
                                 )%>% 
      calc_ESM_VanHaden(T, depth_vals) %>% 
      filter(Lower_cm == depth_of_estimate) %>% #only keep the estimate for the depth we are calculating
      mutate(sample_depths = paste(depth_vals, collapse = ', '),
             method = 'ESM')
    all_res <- rbind(all_res, ESM_res)
  }
  
  for (depth_vals in MC_depths){
    print(depth_vals)
    to_depths <-  rbind(df_agg_to_depths(fd_data %>% filter(ID != Ref_ID), depth_vals),
                        df_agg_to_depths(fd_data %>% filter(ID == Ref_ID), c(depth_of_estimate))
    )
    
    mass_correction_res <- to_depths %>% 
      group_run_MC(adjustment_factor = 1, quantification_depth = depth_of_estimate) %>% 
      mutate(method = "Mass Correction", 
             sample_depths = paste(depth_vals, collapse = ', '))
    
   
    
    all_res <- rbind(all_res, mass_correction_res)
    
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


#'Calculate the mass change to depth of quantification
mass_change <- function(FD_data, quant_depth){
  cum_FD <- calc.cumulative_masses(FD_data) 
  
  return (cum_FD %>% filter(ID != Ref_ID & Lower_cm == quant_depth) %>% 
    merge(cum_FD %>% filter(Lower_cm == quant_depth) %>% select(ID, Cum_Min_Soil_g_cm2), 
          by.x = 'Ref_ID', by.y = 'ID') %>% 
    mutate(mass.change = Cum_Min_Soil_g_cm2.x -  Cum_Min_Soil_g_cm2.y) %>% 
    select(ID, mass.change))
  
}


#'Compare all individual treatment/core data to every possible ref
space_for_time_comparison <- function(FD_data, lat, lon, muname, 
                                      ESM_depths, MC_depths, depth_of_estimate = 30, 
                                      study_name){
  all_out <- data.frame()
  
  for (.ID in unique(FD_data$ID)){
    FD_data <- FD_data %>% mutate(Ref_ID = .ID) 
    
    res<- compare.methods(FD_data, ESM_depths = ESM_depths,
                          MC_depths = MC_depths, 
                          depth_of_estimate = depth_of_estimate)
    
    
    MC_SSurgo <- run_SSurgo_Mass_Corr(lat = lat, lon = lon, muname = muname, 
                                      data = FD_data, 
                                      depth_of_estimate = depth_of_estimate)
    
    res <- res %>% rbind(MC_SSurgo) %>% arrange(ID) %>% 
      tidyr::fill(Cum_SOC_g_cm2_baseline, .direction = 'down') %>% 
      mutate(soc_change = Cum_SOC_g_cm2 - Cum_SOC_g_cm2_baseline, Ref = .ID )
    res <- res %>% merge( mass_change(FD_data, depth_of_estimate), by = 'ID')
    
    all_out <- rbind(all_out, res)
    
  }
  
  max_depths <- all_out %>% filter(nchar(sample_depths) == max(nchar(all_out$sample_depths))) %>% 
    pull(sample_depths) %>% first()
  all_out <- all_out %>% merge(all_out %>% filter(sample_depths == max_depths) 
                       %>% select(ID, soc_change, Ref) %>% 
                         rename(soc_change_benchmark = soc_change)) %>% 
    mutate(dif_from_benchmark = soc_change - soc_change_benchmark) %>% 
    left_join(mass_change(FD_data, depth_of_estimate)) %>% 
    mutate(study = study_name, weight = 1/ (nrow(distinct(FD_data, ID)) -1))
  
  
  summary<- all_out %>% group_by(method, sample_depths) %>% 
    summarize(RMSE_from_benchmark = sqrt(mean(dif_from_benchmark**2)), 
              Bias.from.benchmark = mean(dif_from_benchmark * ifelse(mass.change >0, 1, -1)))
  return (list(summary = summary, res = all_out))
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
  mutate(dif_from_benchmark = soc_change - soc_change_benchmark) %>% merge(mass_change(sanford_fd, 30))


summary<- sanford_res %>% group_by(method, sample_depths) %>% 
  summarize(RMSE_from_benchmark = sqrt(mean(dif_from_benchmark**2)), 
            Bias.from.benchmark = mean(dif_from_benchmark * ifelse(mass.change >0, 1, -1)))



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





franz_FD <- read.csv('Source_data/Franzluebbers_2013.csv')  %>% mutate(Rep = 1)
muname <- "Cecil fine sandy loam, 2 to 6 percent slopes"  
ESM_depths <- list(c(20, 40), c(20,40, 60), c(20, 40, 60, 90), c(40, 60), 
                  c(20, 40 ,60, 90, 120), c(20, 40 ,60, 90, 120, 150))
MC_depths <- list(c(20, 40), c(40))


franz_summary <- space_for_time_comparison(franz_FD, lat= NULL, lon = NULL,
                                           muname = muname, ESM_depths = ESM_depths,
                                           MC_depths = MC_depths, depth_of_estimate = 40)




sainju_FD <- read.csv('Source_data/Sainju_2013.csv') %>% mutate(Rep = 1)
muname <- 'Dooley sandy loam, 0 to 4 percent slopes'
ESM_depths = list(c(15, 30), c(15,30, 60), c(15, 30, 60, 90), c(30, 60), c(7.5, 15, 30, 60),
                  c(7.5, 15, 30, 60, 90, 120), c(7.5, 15,30))
MC_depths = list(c(15, 30), c(30), c(7.5, 15, 30))

sainju_summary <- space_for_time_comparison(sainju_FD, lat = NULL, lon = NULL, 
                                            muname = muname, ESM_depths = ESM_depths,
                                            MC_depths = MC_depths, study_name = 'Sainju_2013.csv'
                                            )



mishra_soil_types = list(NAEW = 'Westmoreland silt loam, 3 to 8 percent slopes', #2-6% slopes
                           WARS = 'Crosby silt loam, 0 to 2 percent slopes',
                         NWARS = 'Hoytville silty clay loam, 0 to 2 percent slopes',
                         Delaware = 'Glynwood silt loam',
                         Coshocton = 'Allegheny silt loam',
                         Hoytville = 'Hoytville clay loam')


mishra_fd <- read.csv('Source_data/mishra_2010.csv') %>% mutate(Rep = 1)
study_name = 'mishra_2010'
mishra_ESM <- list(c(10, 20, 30, 40), c(20, 30), c(10,20,30),
                    c(30, 40), c(20,30, 40))
mishra_MC <- list(c(20, 30), c(10, 20, 30), c(10, 30), c(30))

mishra_fd <- mishra_fd %>% mutate(ID = paste(ID, site, sep = '_'))

res_NAEW <- space_for_time_comparison(mishra_fd %>% filter(site == 'NAEW'),
                                      muname = 'Westmoreland silt loam, 3 to 8 percent slopes',
                                      ESM_depths = mishra_ESM,
                                      MC_depths = mishra_MC, 
                                      lat = NULL, lon = NULL, 
                                      study_name = study_name) 

res_WARS <- space_for_time_comparison(mishra_fd %>% filter(site == 'WARS'),
                                      muname = 'Crosby silt loam, 0 to 2 percent slopes',
                                      ESM_depths = mishra_ESM,
                                      MC_depths = mishra_MC, 
                                      lat = NULL, lon = NULL, 
                                      study_name = study_name) 

res_NWARS <- space_for_time_comparison(mishra_fd %>% filter(site == 'NWARS'),
                                       muname = 'Hoytville silty clay loam, 0 to 1 percent slopes',
                                       ESM_depths = mishra_ESM,
                                       MC_depths = mishra_MC, 
                                       lat = NULL, lon = NULL, 
                                       study_name = study_name) 

res_DEL <- space_for_time_comparison(mishra_fd %>% filter(site == 'Delaware'),
                                     #note: slope data not given in paper
                                     muname = 'Glynwood silt loam, 0 to 2 percent slopes',
                                     ESM_depths = mishra_ESM,
                                     MC_depths = mishra_MC, 
                                     lat = NULL, lon = NULL, 
                                     study_name = study_name) 

res_Cos<- space_for_time_comparison(mishra_fd %>% filter(site == 'Coshocton'),
                                     #note: slope data not given in paper
                                     muname = 'Allegheny silt loam, 0 to 3 percent slopes',
                                     ESM_depths = mishra_ESM,
                                     MC_depths = mishra_MC, 
                                    lat = NULL, lon = NULL, 
                                    study_name = study_name) 

res_Hoyt <- space_for_time_comparison(mishra_fd %>% filter(site == 'Hoytville'),
                                      #note: slope data not given in paper, but this is the only
                                      #muname for Hoytville clay
                                      muname = 'Hoytville clay loam, 0 to 1 percent slopes',
                                      ESM_depths = mishra_ESM,
                                      MC_depths = mishra_MC, 
                                      lat = NULL, lon = NULL,
                                      study_name = study_name)
  

#exclude data/comparisons from woodlots
ag_comparisons_summary <- rbind(res_NAEW$res, res_WARS$res, res_NWARS$res, 
                                res_DEL$res, res_Cos$res, res_Hoyt$res) %>% 
  filter(!grepl('WL', Ref) & !grepl('WL', ID)) %>% 
  group_by(method, sample_depths) %>% summarise(RMSE = sqrt(mean(dif_from_benchmark**2)), 
                                                bias = mean(dif_from_benchmark * ifelse(mass.change>0, 1, -1)))
 
#just the comparisons that include the woodlots                                    
rbind(res_NAEW$res, res_WARS$res, res_NWARS$res, res_DEL$res, res_Cos$res, res_Hoyt$res) %>% 
  filter(Ref == 'WL' | ID == 'WL') %>% 
  group_by(method, sample_depths) %>% summarise(RMSE = sqrt(mean(dif_from_benchmark**2)),
                                                bias = mean(dif_from_benchmark * ifelse(mass.change>0, 1, -1)))



vanDoren_fd <- read.csv('Source_data/Van_doren_1986.csv') %>% mutate(Rep = 1)

muname <- 'Wooster silt loam, 2 to 6 percent slopes'



ESM_depths <- list(c(5, 10, 15,20, 25, 30,35, 40, 45), c(30, 45), c(25, 30, 35), c(15,30, 45), 
                   c(20, 30, 40), c(30, 40),
                   c(10, 20, 30), c(15, 30), c(20, 30))
MC_depths <- list(c(30), c(15,30), c(20, 30), c(25, 30), c(10, 20, 30))



vanDoren_sp_4_t <- space_for_time_comparison(vanDoren_fd %>% filter(!grepl('t1', ID)), 
                                             lat = NULL, lon = NULL, muname = muname,
                                           ESM_depths = ESM_depths, 
                                           MC_depths = MC_depths, study_name = 'Van_doren_1986.csv')


devine_fd <- read.csv('Source_Data/Devine_2014.csv') %>% mutate(Rep = 1)
#Location of Horsehoe-bend research station from google maps
lat <- 33.933715
lon <- -83.356535

devine_res <- space_for_time_comparison(devine_fd, lat = lat, lon = lon, muname = NULL, 
                                        ESM_depths = list(c(5, 15,30, 50, 100), c(5, 15, 30), c(15, 30), c(15,30, 50),
                                        c(15, 30, 50), c(15, 30, 50, 100), c(30,50)),
                                        MC_depths = list(c(5, 15, 30), c(30), c(15, 30)),
                                        study_name = 'Devine_2014'
                                        )

blanco_fd <- read.csv('Source_Data/Blanco_Canqui(2008).csv') %>% mutate(Rep = 1)

#coordinates for the fields representing the various MLRAS
MLRA_coords <- list('122KY' = list(37, 0.1226, 85, 55.5832),
                    '125KY' = list(37, 25.8868, 83, 59.5938),
                    '99OH' = list(41, 21.5594, 83, 5.2101),
                    '124 OH' = list(38,58.379, 82, 47.3865),
                    '139A OH' = list(40, 52.7280, 81, 38.4214),
                    '139B PA' = list(41, 13.2813, 80, 4.5554),
                    '140 PA' = list(41,49.3611, 76,51.7202)
                    )

blanco_res <- data.frame()

for (name in names(MLRA_coords)){
  subset <- filter(blanco_fd, MLRA == name)
  res <- compare.methods(subset, ESM_depths = list(c(5, 10, 30), c(10, 30),  c(30, 50),
                                                   c(10, 30, 50), c(10, 30, 50, 60),
                                                   c(5, 10, 30, 50, 60)
                                                   ),
      MC_depths = list(c(30), c(10, 30))
                     )
  pt_dat <- MLRA_coords[name][[1]]
  
  lat <- pt_dat[[1]] + pt_dat[[2]]/60
  lon <- (pt_dat[[3]] + pt_dat[[4]] /60 ) * -1
  muname <- NULL
  if (name == "139A OH"){
    lat <- NULL
    lon <- NULL
    muname <- 'Chili loam, 2 to 6 percent slopes'
  }
  
  ssurgo_res <- run_SSurgo_Mass_Corr(lat= lat, lon = lon, muname = muname, 
                                     data = subset, depth_of_estimate = 30)
  res <- res %>% rbind(ssurgo_res) %>% arrange(ID) %>% 
    tidyr::fill(Cum_SOC_g_cm2_baseline, .direction = 'down') %>% 
    mutate(soc_change = Cum_SOC_g_cm2 - Cum_SOC_g_cm2_baseline )
  
  blanco_res <- rbind(blanco_res, res) 
}

blanco_res <- blanco_res %>% merge(blanco_res %>% filter(sample_depths == "5, 10, 30, 50, 60") 
                                       %>% select(ID, soc_change) %>% rename(soc_change_benchmark = soc_change)) %>% 
  mutate(dif_from_benchmark = soc_change - soc_change_benchmark) %>% 
  merge(mass_change(blanco_fd, 30)) %>% mutate(study = 'Blanco_Canqui_2008')

blanco_summary <- blanco_res %>% group_by(method, sample_depths) %>% 
  summarize(RMSE_from_benchmark = sqrt(mean(dif_from_benchmark**2)), 
            Bias.from.benchmark = mean(dif_from_benchmark * ifelse(mass.change >0, 1, -1)))


venterea_fd <- read.csv('Source_data/Venterea_2006.csv') %>% mutate(Rep = 1)
muname <- 'Waukegan silt loam, 0 to 2 percent slopes'


venterea_s4_t <- space_for_time_comparison(venterea_fd %>% filter(years == 2000),
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
                                           muname = muname, lat = NULL, lon = NULL, 
                                           study_name = 'Venterea_2006'
                                           
                                           )
chatterjee_munames <- list( '99' = 'Pewamo clay loam, 0 to 2 percent slopes',
                            '111' ='Hoytville silty clay loam, 0 to 1 percent slopes',
                            '124' = 'Genesee loam',
                            #'139' = 'Chili silt loam, 2 to 6 percent slopes',
                            '127' = 'Cookport loam, 3 to 8 percent slopes')

chatterjee_fd <- read.csv('Source_data/chatterjee_2009.csv') %>% mutate(Rep = 1)

chj_res <- data.frame()

for (name in names(chatterjee_munames)){
  subset <- filter(chatterjee_fd, MLRA == name)
  res <- space_for_time_comparison(subset, lat = NULL, lon = NULL,
                                     ESM_depths = list(c(5, 10, 30), c(10, 30),  c(30, 50),
                                                       c(10, 30, 50), c(10, 30, 50, 60),
                                                       c(5, 10, 30, 50, 60)
                                     ),
                                   MC_depths = list(c(30), c(10, 30), c(5,10,30)),
                                   muname = chatterjee_munames[[name]],
                                   study_name = 'chatterjee_2009'
                                   )
    
  chj_res <- rbind(chj_res, res$res) 
}
blanco_munames2 <- list('Hutchinson' = 'Ost loam, 0 to 1 percent slopes',
                        'Hays' = 'Harney silt loam, 0 to 1 percent slopes')

blanco_fd2 <- read.csv('Source_data/Blanco-canqui_2011.csv') %>% mutate(Rep = 1)
blanc_res_2 <- data.frame()

for (name in names(blanco_munames2)){
  subset <- filter(blanco_fd2, site == name)
  res <- compare.methods(subset, ESM_depths = list(c(30,40), c(20,30), c(20,30, 40), c(30,60),
                                                   c(20,30,60), c(20, 30,40, 60),
                                                   c(15, 20, 30,40, 60),
                                                   c(15,30), c(10, 30),
                                                   c(20, 30, 40, 60, 80),
                                                   c(15, 20, 30, 40, 60, 80),
                                                   c(2.5, 5, 10, 15, 20, 30, 40, 60, 80, 100)
  ),
  MC_depths = list(c(30), c(10, 30), c(20, 30), c(15,30))
  )
  
  
  
  ssurgo_res <- run_SSurgo_Mass_Corr(lat= NULL, lon = NULL, muname = blanco_munames2[[name]], 
                                     data = subset, depth_of_estimate = 30)
  res <- res %>% rbind(ssurgo_res) %>% arrange(ID) %>% 
    tidyr::fill(Cum_SOC_g_cm2_baseline, .direction = 'down') %>% 
    mutate(soc_change = Cum_SOC_g_cm2 - Cum_SOC_g_cm2_baseline )
  res <- res %>% merge(res %>% filter(sample_depths == "2.5, 5, 10, 15, 20, 30, 40, 60, 80, 100") 
                              %>% select(ID, soc_change) %>% rename(soc_change_benchmark = soc_change)) %>% 
    mutate(dif_from_benchmark = soc_change - soc_change_benchmark) %>% 
    merge(mass_change(blanco_fd2, 30)) %>% mutate(study = 'Blanco_Canqui_2011',
                                                  Ref = 0, weight = 1)
  blanc_res_2 <- rbind(blanc_res_2, res)
  
}




sanford_st_res <- space_for_time_comparison(sanford_fd %>% filter(ID != Ref_ID & ID != 'all_t0'),
                                            lat = 43.29586, lon = -89.38068, muname = NULL,
                                            ESM_depths = list(c(15, 30), c(15,30, 60), c(15, 30, 60, 90), c(30, 60)),
                                            MC_depths = list(c(15, 30), c(30)),
                                            study_name = 'Sanford_2012'
                                            )

out_s4t <- rbind(res_NAEW$res, res_WARS$res, res_NWARS$res, 
      res_DEL$res, res_Cos$res, res_Hoyt$res, blanco_res %>% mutate(Ref = 0, weight = 1), 
      devine_res$res, vanDoren_sp_4_t$res,
      sainju_summary$res, sanford_st_res$res, venterea_s4_t$res, chj_res, blanc_res_2)


summarize_res <- function(grouped_res){
  grouped_res %>% summarise(RMSE = sqrt(weighted.mean(dif_from_benchmark**2, weight)), 
                            RMPSE = sqrt(weighted.mean((dif_from_benchmark /soc_change_benchmark )**2, weight )), 
                            bias = weighted.mean(dif_from_benchmark * ifelse(mass.change >0, 1, -1), weight), 
                            MAPE = weighted.mean(abs(dif_from_benchmark/soc_change_benchmark), weight), ,
                            Max_err =  max(abs(dif_from_benchmark)),
                            mean_abs_err = mean(abs(dif_from_benchmark)))
}


cross_fold_weights_res <- data.frame()
for (study_name in unique(out_s4t$study)){
  training_data <- out_s4t %>% filter(study != study_name)
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
   testing_data %>% filter(method %in% c('Single Depth', 'Mass Correction') & sample_depths == 30) %>% 
     pivot_wider(id_cols = c(ID, Ref), names_from = method, values_from = soc_change) %>% 
     mutate(soc_change = `Mass Correction` * MC_weight + `Single Depth` * sd_weight ) %>% 
     merge(testing_data %>% select(ID, Ref, soc_change_benchmark, mass.change, weight) %>% distinct(ID, Ref, .keep_all = T)) %>%
     mutate(dif_from_benchmark = soc_change - soc_change_benchmark, study = study_name,
            method = 'Mass Correction, empirical weights') %>%
     rbind(cross_fold_weights_res ) -> cross_fold_weights_res
   
}





vanDoren_res<- compare.methods(vanDoren_fd, ESM_depths, MC_depths)

MC_SSurgo <- run_SSurgo_Mass_Corr(lat = NULL, lon = NULL,
                                  muname,
                                  data = vanDoren_fd )

vanDoren_res <- vanDoren_res %>% rbind(MC_SSurgo) %>% arrange(ID) %>% 
  tidyr::fill(Cum_SOC_g_cm2_baseline, .direction = 'down') %>% 
  mutate(soc_change = Cum_SOC_g_cm2 - Cum_SOC_g_cm2_baseline )


vanDoren_res <- vanDoren_res %>% merge(vanDoren_res %>% filter(sample_depths == "5, 10, 15, 20, 25, 30, 35, 40, 45") 
                                       %>% select(ID, soc_change) %>% rename(soc_change_benchmark = soc_change)) %>% 
  mutate(dif_from_benchmark = soc_change - soc_change_benchmark) %>% merge(mass_change(vanDoren_fd, 30))


summary<- vanDoren_res %>% group_by(method, sample_depths) %>% 
  summarize(RMSE_from_benchmark = sqrt(mean(dif_from_benchmark**2)), 
            Bias.from.benchmark = mean(dif_from_benchmark * ifelse(mass.change >0, 1, -1)))


calced_difs <- vanDoren_res  %>% group_by(method, sample_depths) %>% 
  reframe(changes = dif_matrix(soc_change))

as.df <- as.data.frame(calced_difs$changes) %>% 
  mutate(method = calced_difs$method, sample_depths = calced_difs$sample_depths)
calced_difs <- as.df

benchmark_difs <- as.df %>% filter(sample_depths == "5, 10, 15, 20, 25, 30, 35, 40, 45") %>% select(-c(method, sample_depths))


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
