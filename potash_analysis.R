library(sf)
library(parallel)
source('analysis_functions.R')

potash_data <- read.csv('Source_Data/Potash_et_al_data/measurements.csv') %>% 
      rename(Upper_cm = sample_depth_min, 
                           Lower_cm = sample_depth_max,
                           SOC_pct = SOCc, 
                           BD_g_cm3 = BD, 
                          ID = location_id) %>% 
       mutate(SOM_pct = SOC_pct * 1.9, 
                           Rep = 1, 
                           Ref_ID = '1') %>% 
  mutate(ID = paste(site, ID, sep = '_')) %>% 
  calc.cumulative_masses() %>%
 mutate(study_name = site)


locs <- read.csv('Source_Data/Potash_et_al_data/locations.csv')

locs <- st_as_sf(locs, coords = c('X', 'Y'), remove = F) %>% 
  mutate(ID = paste(site, location_id, sep = '_'))

d <- st_distance(locs)

closest <- apply(d, 2, function(arr){sort(arr, index.return = T)$ix[2:11]})



fun <- function(i){
  out <- data.frame()
  .ID <- locs[i,] %>% pull(ID)
  lat <- locs[i,] %>% pull(Y)
  lon <- locs[i,] %>% pull(X)
  sub_data <- filter(potash_data, ID == .ID)
  for (i2 in 1:3){
    ref_I <- closest[i2, i]
    .Ref_ID <- locs[ref_I,] %>% pull(ID)
    compare_data <-  filter(potash_data, ID == .Ref_ID)
    if (length(setdiff(compare_data$Lower_cm, sub_data$Lower_cm))==0  ){
      
      run_data <- bind_rows(sub_data, compare_data) %>%
        mutate(Ref_ID =.Ref_ID)
      
      res <- tryCatch({compare.methods(run_data, lat = lat, lon = lon)},
                      error = function(cond){data.frame()})
      out <- bind_rows(out, res)}
  }
  
  return(out)
  }

cl <- makeCluster(detectCores() -1)
clusterExport(cl, c('fun', 'locs', 'potash_data', 'closest'), envir = environment())
clusterEvalQ(cl, {library(dplyr)})
clusterEvalQ(cl, {library(tidyr)})
clusterEvalQ(cl, {source('analysis_functions.R')})
#clusterExport(cl, c('check_errors', 'spline_error_within', 'sum_weights_mass',  'sum_weights_mass_2',
#                    'hyman_extrapolating_splinefun', "make_zerodepth_row", "calc.linear_SOC_pct",
#                    'df_agg_to_depths', 'agg_to_depths', 'agg_to_depth',
#                    'calc.cumulative_masses'))

out <- parLapply(cl, 1:nrow(locs), fun) %>% bind_rows()


save_path <- file.path('Results', 'potash_data_comparisons.csv')
write.csv(out, save_path)


out <- mutate(out, study_name = sub('_.*', "", ID))


sub_potash <- filter(potash_data, ID %in% unique(out$ID) | ID %in% unique(out$Ref_ID)  )
get_fd_data <- function(id){filter(sub_potash, study_name == id)}

out <- bind_rows(lapply(unique(out$study_name), function(x){
  add_reference_changes(x, out, load_fd_data_func = get_fd_data, nsims = 50)}))
# 
# out <- data.frame()
# for (i in 1:nrow(locs)){
#   print(i)
#   .ID <- locs[i,] %>% pull(ID)
#   lat <- locs[i,] %>% pull(Y)
#   lon <- locs[i,] %>% pull(X)
#   sub_data <- filter(potash_data, ID == .ID)
#   for (i2 in 1:3){
#     ref_I <- closest[i2, i]
#     .Ref_ID <- locs[ref_I,] %>% pull(ID)
#     compare_data <-  filter(potash_data, ID == .Ref_ID)
#     if (length(setdiff(compare_data$Lower_cm, sub_data$Lower_cm))==0  ){
#     run_data <- bind_rows(sub_data, compare_data) %>%
#       mutate(Ref_ID =.Ref_ID)
#     
#   res <- compare.methods(run_data, lat = lat, lon = lon)
#   out <- bind_rows(out, res)
#   }}}

get_IDs_for_analysis <- function(res_data, depths1, depths2){
  IDs1 <- res_data %>% filter(sample_depths == depths1) %>% pull(ID) %>% unique()
  IDs2 <- res_data %>% filter(sample_depths == depths2) %>% pull(ID) %>% unique()

  intersect(IDs1, IDs2)}





summary_potash <- out %>% 
  group_by(method, sample_depths, study_name, simulation.N) %>% 
  mutate(weight = 1) %>% 
  summarize_res() %>% 
  mutate(data_type = 'Single Core Data')


summary_fieldtrials <- ag_res %>% 
  group_by(study_name, method, sample_depths, simulation.N) %>% 
  summarize_res() %>%
  mutate(data_type = 'Field-Level Aggregated Data')

onedepth_sum <-  rbind(summary_potash, summary_fieldtrials) %>% filter(sample_depths == '30') %>% 
  mutate(method = factor(method)) %>% 
  filter(method != 'ESM Spline') %>% 
  mutate(method_group = case_when(
    method %in% c('Fixed Depth', 'Linear') ~ "Raw Data Methods",
    method %in% c('Linear Average', 'Linear Average, Fowler Weights') ~ 'Averaging Methods',
    grepl('SSurgo', method) ~ 'SSurgo-based methods'
  )) %>% 
  mutate(method_group = factor(method_group, 
                               levels = c('Raw Data Methods', 'Averaging Methods', 'SSurgo-based methods'))) %>%
  mutate(method = factor(method, levels = c("Fixed Depth",   "Linear",   
                                             "Linear Average",                 
                                             "Linear Average, Fowler Weights",
                                             "Mass Correction, SSurgo avg", "Mass Correction, SSurgo EXP",     
                                          "Mass Correction, SSurgo_spline")))


  ggplot(onedepth_sum, aes(x = method_group, y= bias * 100, color = method#, fill = method#, 
                           #shape = data_type
                           )) + 
    geom_boxplot(position = position_dodge(.75) ) + 
    geom_hline(yintercept = 0, linetype = 'dashed') +
  scale_y_continuous(name = 'Bulk-Density Bias, T C/ha ', breaks = c(-2,0,2,4), 
                     label = c(-2,0,2,4), limits = c(-2, 4.5)) + xlab('') + 
  ggtitle('Bias of Different Methods of Estimating Soil Organic Carbon Change \nfrom Single-Depth Samples, by Site/Study and Simulation')
  
  ggplot(onedepth_sum, aes(x = method_group, y= bias_rate, color = method, fill = method)) + 
    geom_point(position = position_dodge(.5) ) + 
    geom_hline(yintercept = 0, linetype = 'dashed') +
    scale_y_continuous(name = 'Bias Rate: Tons C / Tons Mineral Soil', breaks = c(-.01, 0, .01, .02), 
                       label = c(-.1, 0, .1, .2), limits = c(-.01, .02)) + 
    xlab('') + 
    ggtitle('Bias, as a proportion of change in mass to 30 cm, \n of Different Methods of Estimating Soil Organic Carbon Change \nfrom Single-Depth Samples,by Site/Study and Simulation')
  

  ggplot(onedepth_sum, aes(x = method_group, y= RMSE * 100, color = method)) + 
    geom_boxplot(position = position_dodge(.75) ) + 
    scale_y_continuous(name = 'RMSE, T C/ha', breaks = 0:5, labels = 0:5) + xlab('') + 
    ggtitle('Error of Different Methods of Estimating Soil Organic Carbon Change \nfrom Single-Depth Samples, by Site/Study and Simulation')
  
  
  ggplot(onedepth_sum, aes(x = method_group, y = MAPE_change * 100, color = method, fill = method)) + 
    geom_boxplot(position = position_dodge(.5) ) + 
    #geom_hline(yintercept = 0, linetype = 'dashed') +
    scale_y_continuous(trans = 'log', breaks = c(10, 100, 1000, 2000)) +
    xlab('') + 
    ggtitle('Bias, as a proportion of change in mass to 30 cm, \n of Different Methods of Estimating Soil Organic Carbon Change \nfrom Single-Depth Samples, by Site/Study')
  
  
  plot_for_present_in_all <- function(sample_depths1, sample_depths2){
  
    ids_potash <- get_IDs_for_analysis(out, sample_depths1, sample_depths2)
    
    
    #out_s4t <- read.csv( file.path('Results', 'all_data_space_for_time_comparison.csv'))
    
    ids_other <- get_IDs_for_analysis(out_s4t, sample_depths1, sample_depths2)
    
  sum2.potash <- out %>% filter(ID %in% ids_potash) %>% mutate(weight =1) %>% 
    group_by(study_name, method, sample_depths, simulation.N) %>% 
    summarize_res()
  
  sum2.other <- out_s4t %>% filter(ID %in% ids_other) %>% 
    group_by(study_name, method, sample_depths, simulation.N) %>% 
    summarize_res()
  
  
  sum2_to_plot <- bind_rows(sum2.potash, sum2.other)
  
  
  
  
  p <- sum2_to_plot %>% filter(sample_depths %in% c(sample_depths1, '30', sample_depths2)) %>% 
    filter(sample_depths!= '30' | grepl('SSurgo', method) | method == 'Fixed Depth'| method == 'Linear') %>% 
    ggplot(aes(x = sample_depths, y = RMSE *100, color = method)) + 
    geom_boxplot(position = position_dodge(.75), alpha = .5) + 
    #geom_hline(yintercept = 0, linetype = 'dashed')
    ylab('RMSE T C/ha') + xlab('Depths sampled') +
    ggtitle('RMSE by sampling depth and quantification method \n aggregated by site/study and simulation')
  print(p)
  
  sum2_to_plot %>% filter(sample_depths %in% c(sample_depths1, '30', sample_depths2)) %>% 
    filter(sample_depths!= '30' | grepl('SSurgo', method) | method == 'Fixed Depth' | method == 'Linear') %>% 
    ggplot(aes(x = sample_depths, y = bias *100, color = method)) + 
    geom_boxplot(position = position_dodge(.75)) + 
    geom_hline(yintercept = 0, linetype = 'dashed') +
    ylab('Bulk Density bias, T C/ha') + xlab('Depths sampled') +
    ggtitle('Bias by sampling depth and quantification method \n aggregated by site/study and simulation')
  
  }
  
plot_for_present_in_all('15, 30', '30, 45')

plot_for_present_in_all('15, 30', '30, 60')
  